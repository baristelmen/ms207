import numpy as np
import matplotlib.pyplot as plt
import gmsh
from dolfinx.fem import (
    locate_dofs_topological,
    functionspace,
    dirichletbc,
    Constant,
    Function,
    Expression,
)

from domaine import Domaine
from modeles import ModeleMecaniqueLineaireElastiqueIsotrope
from resolution import resolution
from affichage import (
    affichage_fonction_vectorielle_2D,
    affichage_deformee_2D,
    affichage_fonction_tensorielle_2D,
)
from post_traitements import extraction_contrainte_sur_ligne

# ----------- DONNEES-------------#
# Materiau
YOU = 10_000.0
v = 0.3
# Chargement
sig_inf = 1.0
# Geometrie
a0 = 1.0
big = 10.0
# big = 20.0
l0 = big * a0
h0 = big * a0
# Numerique
ordre_elem = 1
# Maillage
# de1 = 0.01 * a0
# de2 = 0.5 * a0
de1 = 0.001 * a0
de2 = 0.1 * a0
# NOTE: Il faut aussi raffiner de2 sinon ça ne change pas beaucoup le maillage.
#   Sinon on peut rajouter un "size field" mais ça complexifie un peu le code

# Desactivation des tracés
afficher_les_traces = False


# ----------- GEOMETRIE ET MAILLAGE -------------#
# Parametre
dim = 2
# Initialize GMSH
gmsh.initialize()
# Choix de l'algorithme de maillage (ici, 5 = Delaunay)
gmsh.option.setNumber("Mesh.Algorithm", 5)
# Define points
p0 = gmsh.model.geo.addPoint(0, 0, 0, meshSize=de2)
p1 = gmsh.model.geo.addPoint(a0, 0, 0, meshSize=de1)  # Crack tip
p2 = gmsh.model.geo.addPoint(l0, 0, 0, meshSize=de2)
p3 = gmsh.model.geo.addPoint(l0, h0, 0, meshSize=de2)
p4 = gmsh.model.geo.addPoint(0, h0, 0, meshSize=de2)
# Define lines
d01 = gmsh.model.geo.addLine(p0, p1)
d12 = gmsh.model.geo.addLine(p1, p2)
d23 = gmsh.model.geo.addLine(p2, p3)
d34 = gmsh.model.geo.addLine(p3, p4)
d40 = gmsh.model.geo.addLine(p4, p0)
# Create curve loop and surface
cl = gmsh.model.geo.addCurveLoop([d01, d12, d23, d34, d40])
sur = gmsh.model.geo.addPlaneSurface([cl])
# Add physical groups (for boundary conditions
gmsh.model.geo.addPhysicalGroup(dim, [sur])
gmsh.model.geo.addPhysicalGroup(dim - 1, [d12], tag=d12)
gmsh.model.geo.addPhysicalGroup(dim - 1, [d34], tag=d34)
gmsh.model.geo.addPhysicalGroup(dim - 1, [d40], tag=d40)
# Synchronize to apply the changes
gmsh.model.geo.synchronize()
# Generate mesh
gmsh.model.mesh.generate(dim)

# Show the mesh (comment to hide)
# gmsh.fltk.run()

# Creation du maillage à partir à partir du modèle GMSH
domaine = Domaine(gmsh.model)
if afficher_les_traces:
    domaine.affichage()
    # NOTE: Contrôles de la fenêtre interactive: https://docs.pyvista.org/api/plotting/plotting

# Stop gmsh
gmsh.finalize()

# ----------- DEFINITION DE MODELE ------------- #
modele = ModeleMecaniqueLineaireElastiqueIsotrope(E=YOU, nu=v)

# ----------- DEFINITION DU CHAMPS DE DÉPLACEMENT ------------- #
# Définition de l'espace de fonction pour le champ de déplacement
Vu = functionspace(domaine.maillage, ("Lagrange", ordre_elem, (dim,)))
Vux = Vu.sub(0)
Vuy = Vu.sub(1)
# NOTE:  Cette section n'est pas dans CAST3M (on peut la cacher)
#   Par exemple, mettre ça dans le modèle (avec la définition de u)


# ----------- DEFINITION DU CHARGEMENT ------------- #
# Définition du déplacement imposé d12
d12_fm = domaine.facet_markers.find(d12)
d12_dofs = locate_dofs_topological(Vuy, dim - 1, d12_fm)
cl1 = dirichletbc(Constant(domaine.maillage, 0.0), d12_dofs, Vuy)
# Définition du déplacement imposé d40
d40_fm = domaine.facet_markers.find(d40)
d40_dofs = locate_dofs_topological(Vux, dim - 1, d40_fm)
cl2 = dirichletbc(Constant(domaine.maillage, 0.0), d40_dofs, Vux)
# Définition de l'effort imposée sur d34
cha = {"T": [0.0, sig_inf], "bord": d34}
# NOTE: On peut faire des fonctions pour rendre plus clair ce qui précède

# ----------- CALCUL DE LA SOLUTION ---------------- #
# Calcul du deplacement solution
uu = resolution(domaine, modele, [cl1, cl2], [cha], Vu)

# ----------- POST-TRAITEMENT ---------------- #
# Affichage du champ de déplacement
if afficher_les_traces:
    affichage_fonction_vectorielle_2D(domaine, uu)
    affichage_deformee_2D(domaine, uu)
# Calcul des contraintes
# NOTE: On ne calcule pas la composante sig_zz mais elle est non-nulle.
Vsig = functionspace(domaine.maillage, ("DG", 0, (2, 2)))
sig_expr = Expression(modele.sig(uu), Vsig.element.interpolation_points())
sig_func = Function(Vsig)
sig_func.interpolate(sig_expr)
# TODO: Peut-être cacher un peu ça ?

# Affichage du champs de contrainte
if afficher_les_traces:
    affichage_fonction_tensorielle_2D(domaine, sig_func)

# Extraction du champ de contrainte en y=0
N = 1_000
ligne = np.array([(0.0, l0 * i / (N - 1), 0.0) for i in range(N)])
sig_sur_ligne = extraction_contrainte_sur_ligne(sig_func, ligne)

# Affichage de la contrainte le long de la ligne y=0
if afficher_les_traces:
    plt.figure()
    plt.xlabel(r"Coordonnée $x$ (m)")
    plt.ylabel(r"Contrainte $\sigma_{yy}$ en $y=0$ (Pa)")
    plt.scatter(ligne[:, 1], sig_sur_ligne[:, 3], marker="x")
    plt.grid()
    plt.show()

# Détermination de KI pour interpolation sur le champ de déplacement
xs = domaine.maillage.geometry.x
# Création d'un masque pour extraire les points en y=0
masque_y0 = np.isclose(xs[:, 1], 0)

# Extraction de la composante selon y du champ de déplacement
composante_y = [i for i in range(1, uu.x.array.shape[0], 2)]
uy = uu.x.array[composante_y]

# Extraction de la
x_en_y0 = xs[masque_y0, 0]
uy_en_y0 = uy[masque_y0]

# Tri des données selon x croissant
masque_de_tri = np.argsort(x_en_y0)
x_en_y0 = x_en_y0[masque_de_tri]
uy_en_y0 = uy_en_y0[masque_de_tri]


# Determination de KI par interpolation
N_point_calcul_KI = 3
masque_calcul_KI = x_en_y0 < a0

if False:
    uy_calcul_KI = uy_en_y0[masque_calcul_KI][-N_point_calcul_KI:]
    r_calcul_KI = abs(x_en_y0[masque_calcul_KI] - a0)[-N_point_calcul_KI:]
    Cs = uy_calcul_KI / np.sqrt(r_calcul_KI)
    C = np.average(Cs)
    KI_interp = C * 2 * modele.mu * np.sqrt(2 * np.pi) / (modele.k + 1)
else:
    uy_calcul_KI_sqrd = uy_en_y0[masque_calcul_KI][-N_point_calcul_KI:] ** 2
    r_calcul_KI = abs(x_en_y0[masque_calcul_KI] - a0)[-N_point_calcul_KI:]
    res = np.polynomial.Polynomial.fit(r_calcul_KI, uy_calcul_KI_sqrd, deg=1)
    C1 = res.convert().coef[1]
    KI_interp = np.sqrt(C1) * 2 * modele.mu * np.sqrt(2 * np.pi) / (modele.k + 1)

# Détermination de K1 théorique
KI_th = sig_inf * np.sqrt(np.pi * a0)


plt.figure()
plt.scatter(x_en_y0, uy_en_y0, label="Éléments finis")
plt.scatter(
    x_en_y0,
    KI_interp
    / (2 * modele.mu)
    * np.sqrt((a0 - x_en_y0) / (2 * np.pi))
    * (modele.k + 1),
    color="r",
    label="Reconstruit via $K_I$ interpolé",
)
plt.scatter(
    x_en_y0,
    KI_th / (2 * modele.mu) * np.sqrt((a0 - x_en_y0) / (2 * np.pi)) * (modele.k + 1),
    color="g",
    label="Reconstruit via $K_I$ théorique",
)
plt.xlabel("Coordonnée $x$ (m)")
plt.ylabel("Déplacement selon e_y en x=0 (m)")
plt.xlim([0.0, 1.5 * a0])
plt.legend()
plt.grid()
plt.show()


print(f"KI interpolé : {KI_interp} Pa.m^1/2")
print(f"KI théorique : {KI_th} Pa.m^1/2")
