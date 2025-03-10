import numpy as np
from ufl import (
    TestFunction,
    TrialFunction,
    inner,
    dot,
    grad,
    div,
    dx,
    as_vector,
    cos,
    sin,
)
from dolfinx import default_scalar_type
from dolfinx.geometry import bb_tree, compute_collisions_points, compute_colliding_cells
from dolfinx.fem import (
    functionspace,
    Constant,
    Function,
    Expression,
    locate_dofs_geometrical,
    dirichletbc,
    form,
    assemble_scalar,
)
from dolfinx.fem.petsc import LinearProblem


def calcul_de_la_contrainte(domaine, modele, uu):
    Vsig = functionspace(domaine.maillage, ("DG", 0, (2, 2)))
    sig_expr = Expression(modele.sig(uu), Vsig.element.interpolation_points())
    sig_func = Function(Vsig)
    sig_func.interpolate(sig_expr)
    return sig_func


def extraction_contrainte_sur_ligne(sig_func, line):
    # Generate the bounding box tree
    mesh = sig_func.function_space.mesh
    tree = bb_tree(mesh, mesh.topology.dim)
    # Find cells whose bounding-box collide with the the points
    cell_candidates = compute_collisions_points(tree, line)
    # For each points, choose one of the cells that contains the point
    colliding_cells = compute_colliding_cells(mesh, cell_candidates, line)
    cells = [colliding_cells.links(i)[0] for i, x in enumerate(line)]
    # Evaluate the stress function along the line
    sig_along_line = sig_func.eval(line, cells)
    return sig_along_line


def calcul_KI_via_interp_uy(domaine, modele, uu, a0):
    # NOTE: Try this https://www.fracturemechanics.org/kfromdisps.html
    # Détermination de KI pour interpolation sur le champ de déplacement

    # TODO: Try what is done here: https://www-cast3m.cea.fr/index.php?page=procedures&procedure=sif
    xs = domaine.maillage.geometry.x
    # Création d'un masque pour extraire les points en y=0 et x<a0
    masque_theta_180 = np.isclose(xs[:, 1], 0) & (0.9 * a0 < xs[:, 0]) & (xs[:, 0] < a0)

    # Extraction de la composante selon y du champ de déplacement
    composante_y = [i for i in range(1, uu.x.array.shape[0], 2)]
    uy = uu.x.array[composante_y]

    # Extraction du déplacement vertical le long de la fissure
    r_en_y0 = abs(xs[masque_theta_180, 0] - a0)
    uy_en_y0 = uy[masque_theta_180]

    # DEBUG
    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.xlabel("Distance to crack tip $r$")
    # plt.ylabel(r"Displacement $u_y(r, \theta = \pi)$")
    # plt.scatter(r_en_y0, uy_en_y0**2)
    # plt.show()

    # Linear regression
    pol = np.polynomial.Polynomial.fit(r_en_y0, uy_en_y0**2, 2)
    A = pol.convert().coef[1]
    KI_interp = 2 * modele.mu / (modele.k + 1) * np.sqrt(2 * np.pi * A)

    return KI_interp


def calcul_du_champ_theta(domaine, pointe_fissure, R_int, R_ext):
    # Define the distance to the crack tip
    def r(x):
        return np.sqrt(
            (x[0] - pointe_fissure[0]) ** 2 + (x[1] - pointe_fissure[1]) ** 2
        )

    # Define the variational problem to define theta
    V_theta = functionspace(domaine.maillage, ("Lagrange", 1))
    theta, theta_ = TrialFunction(V_theta), TestFunction(V_theta)
    a = dot(grad(theta), grad(theta_)) * dx
    L = Constant(domaine.maillage, default_scalar_type(0.0)) * theta_ * dx
    # Set the boundary conditions
    # Imposing 1 in the inner circle and zero in the outer circle
    dofs_inner = locate_dofs_geometrical(V_theta, lambda x: r(x) <= R_int)
    dofs_out = locate_dofs_geometrical(V_theta, lambda x: r(x) >= R_ext)
    bc_inner = dirichletbc(default_scalar_type(1.0), dofs_inner, V_theta)
    bc_out = dirichletbc(default_scalar_type(0.0), dofs_out, V_theta)
    bcs = [bc_out, bc_inner]
    # Solve the problem
    problem = LinearProblem(a, L, bcs=bcs)
    return problem.solve()


def calcul_G_via_G_theta(domaine, modele, uu, pointe_fissure, R_int, R_ext, phi0=0.0):
    # Definition de l'amplitude du champ theta
    theta_amp = calcul_du_champ_theta(domaine, pointe_fissure, R_int, R_ext)

    # DEBUG
    # from affichage import affichage_fonction_scalaire_2D
    # affichage_fonction_scalaire_2D(domaine, theta_amp)

    # Definition de l'expression de l'intégrale de volume
    eps = modele.eps(uu)
    sig = modele.sig(uu)
    theta_vector = as_vector([cos(phi0), sin(phi0)]) * theta_amp
    G_expr = (
        inner(sig, grad(uu) * grad(theta_vector)) * dx
        - 1 / 2 * inner(sig, eps) * div(theta_vector) * dx
    )
    # Calcul de G
    return assemble_scalar(form(G_expr))
