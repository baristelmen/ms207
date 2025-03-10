from dolfinx.fem import Constant, Function
from dolfinx.fem.petsc import LinearProblem
from ufl import Measure, dot, derivative, replace, TestFunction, TrialFunction, lhs, rhs


def resolution(domaine, modele, bcs, tractions, Vu):
    # Definition des mesures
    dx = Measure("dx", domain=domaine.maillage, subdomain_data=domaine.cell_markers)
    ds = Measure("ds", domain=domaine.maillage, subdomain_data=domaine.facet_markers)
    # Definition de la fonction de déplacement
    u = Function(Vu, name="Déplacement")
    # Calcul de l'énergie élastique
    energie_elastique = modele.energie_elastique(u)
    # Calcul du travail des efforts exterieurs
    fv = Constant(Vu.mesh, [0.0, 0.0])
    travail_exterieur = dot(fv, u) * dx
    for traction in tractions:
        # Extraction des données
        T_vec = traction["T"]
        b_id = traction["bord"]
        # Définition du vecteur de traction
        T = Constant(Vu.mesh, T_vec)
        # Ajout de la contribution au travail des efforts extérieurs
        travail_exterieur += dot(T, u) * ds(b_id)
    # Calcul de l'énergie potentielle totale
    P = energie_elastique - travail_exterieur

    # Dérivation de l'énergie potentielle (pour écrire l'équation de stationnarité)
    dP_du_temp = derivative(P, u, TestFunction(Vu))
    dP_du = replace(dP_du_temp, {u: TrialFunction(Vu)})

    # Définition du problème élastique
    problem = LinearProblem(
        lhs(dP_du),
        rhs(dP_du),
        bcs=bcs,
        petsc_options={
            "ksp_type": "preonly",
            "pc_type": "cholesky",
            "pc_factor_mat_solver_type": "cholmod",
        },
        u=u,
    )
    # Résolution du problème
    problem.solve()
    # Retourne le champ de déplacement solution
    return u
