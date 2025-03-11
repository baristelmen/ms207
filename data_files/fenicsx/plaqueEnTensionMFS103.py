# ----------- IMPORT DES MODULES EXTERNES -------------#
import numpy as np
import gmsh
from ufl import Measure, FacetNormal, dot
from dolfinx.fem import (
    locate_dofs_topological,
    functionspace,
    dirichletbc,
    Constant,
    form,
    assemble_scalar,
)

# ----------- IMPORT DES MODULES CRÉÉS POUR LE COURS -------------#
from domaine import Domaine
from modeles import ModeleMecaniqueLineaireElastiqueIsotrope
from resolution import resolution
from post_traitements import (
    calcul_de_la_contrainte,
    calcul_KI_via_interp_uy,
    calcul_G_via_G_theta,
)
from affichage import (
    affichage_fonction_vectorielle_2D,
    affichage_deformee_2D,
    affichage_contour_deformee_2D,
    affichage_fonction_tensorielle_2D,
)

# ----------- DONNEES -------------#
# Matériau
YOU = 1.0
v = 0.4
# Chargement
udon = 1.0
# Géometrie (mm)
hypothese_2D = "contraintes_planes"
dim = 2
a0 = 20 / 100
l0 = 1.0
h0 = 4.5 / 2
# Maillage
f0 = 0.25
# f0 = 0.5
de1 = f0 * 0.01 * a0
de2 = f0 * 0.5 * a0
# G-theta
R_int = a0 / 8
R_ext = a0 / 4
# Numerique
ordre_elem = 1  # NOTE: Ne pas changer
# Desactivation des tracés
# NOTE: Contrôles des fenêtres interactives: https://docs.pyvista.org/api/plotting/plotting
afficher_les_traces = True

# ----------- GEOMETRIE ET MAILLAGE -------------#
# Initialisation de GMSH
gmsh.initialize()
# Définition des points
p0 = gmsh.model.geo.addPoint(0, 0, 0, meshSize=de2)
p1 = gmsh.model.geo.addPoint(a0, 0, 0, meshSize=de2)  # Pointe de la fissure
p2 = gmsh.model.geo.addPoint(l0, 0, 0, meshSize=de2)
p3 = gmsh.model.geo.addPoint(l0, h0, 0, meshSize=de2)
p4 = gmsh.model.geo.addPoint(0, h0, 0, meshSize=de2)
# Définition des lignes
d01 = gmsh.model.geo.addLine(p0, p1)
d12 = gmsh.model.geo.addLine(p1, p2)
d23 = gmsh.model.geo.addLine(p2, p3)
d34 = gmsh.model.geo.addLine(p3, p4)
d40 = gmsh.model.geo.addLine(p4, p0)
# Définition des surfaces
cl = gmsh.model.geo.addCurveLoop([d01, d12, d23, d34, d40])
sur = gmsh.model.geo.addPlaneSurface([cl])
# Création des groupes physiques (pour les conditions aux limites)
gmsh.model.geo.addPhysicalGroup(dim, [sur])
gmsh.model.geo.addPhysicalGroup(dim - 1, [d12], tag=d12)
gmsh.model.geo.addPhysicalGroup(dim - 1, [d34], tag=d34)
gmsh.model.geo.addPhysicalGroup(dim - 1, [d40], tag=d40)
# Synchronisation du modèle GMSH
gmsh.model.geo.synchronize()
# Raffinement du maillage autour de la pointe de fissure
gmsh.model.mesh.field.add("Distance", 1)
gmsh.model.mesh.field.setNumbers(1, "PointsList", [p1])
gmsh.model.mesh.field.add("Threshold", 2)
gmsh.model.mesh.field.setNumber(2, "IField", 1)
gmsh.model.mesh.field.setNumber(2, "LcMin", de1)
gmsh.model.mesh.field.setNumber(2, "LcMax", de2)
gmsh.model.mesh.field.setNumber(
    2, "DistMin", R_ext
)  # Distance sous laquelle appliquer le maillage fin
gmsh.model.mesh.field.setNumber(
    2, "DistMax", 2 * R_ext
)  # Distance au-delà de laquelle appliquer le maillage grossier
gmsh.model.mesh.field.setAsBackgroundMesh(2)

# # NOTE: Pour debugger lors de la création de maillage
# # Display and exit for debug purposes
# # Synchronize the model
# gmsh.model.geo.synchronize()
# # Display the GMSH window
# gmsh.fltk.run()
# exit()

# Génération du maillage
gmsh.model.mesh.generate(dim)

# Création du maillage à partir à partir du modèle GMSH
domaine = Domaine(gmsh.model)
if afficher_les_traces:
    # NOTE: Contrôles de la fenêtre interactive: https://docs.pyvista.org/api/plotting/plotting
    domaine.affichage()

# Arrêt de gmsh
gmsh.finalize()

# ----------- DEFINITION DE MODELE ------------- #
modele = ModeleMecaniqueLineaireElastiqueIsotrope(
    E=YOU, nu=v, hypothese_2D=hypothese_2D
)

# ----------- DISCRETISATION DU CHAMP DE DÉPLACEMENT ------------- #
# Définition de l'espace (discret) de fonctions pour le champ de déplacement
Vu = functionspace(domaine.maillage, ("Lagrange", ordre_elem, (dim,)))
Vux = Vu.sub(0)
Vuy = Vu.sub(1)

# ----------- DEFINITION DU CHARGEMENT ------------- #
# Définition du déplacement imposé d12 (plan de symétrie)
d12_fm = domaine.facet_markers.find(d12)
d12_dofs = locate_dofs_topological(Vuy, dim - 1, d12_fm)
cl1 = dirichletbc(Constant(domaine.maillage, 0.0), d12_dofs, Vuy)
# Définition du déplacement imposé d34
d34_fm = domaine.facet_markers.find(d34)
d34_dofs = locate_dofs_topological(Vu, dim - 1, d34_fm)
cl2 = dirichletbc(Constant(domaine.maillage, [0.0, udon / 2]), d34_dofs, Vu)
# Groupements de conditions aux limites en déplacement
us_imp = [cl1, cl2]
# Pas d'effort imposé
Ts_imp = []

# ----------- CALCUL DE LA SOLUTION ---------------- #
# Calcul du deplacement solution
uu = resolution(domaine, modele, us_imp, Ts_imp, Vu)

# ----------- POST-TRAITEMENT ---------------- #
# Affichage du champ de déplacement
if afficher_les_traces:
    affichage_fonction_vectorielle_2D(domaine, uu)
    affichage_deformee_2D(domaine, uu, amplitude=0.47)
    affichage_contour_deformee_2D(domaine, uu, amplitude=0.47)
# Calcul des contraintes
sig = calcul_de_la_contrainte(domaine, modele, uu)

# Affichage du champs de contrainte
if afficher_les_traces:
    affichage_fonction_tensorielle_2D(domaine, sig)

# Calcul de l'effort de réaction
ds = Measure("ds", domain=domaine.maillage, subdomain_data=domaine.facet_markers)
n = FacetNormal(domaine.maillage)
ey = Constant(domaine.maillage, [0.0, 1.0])
T = dot(modele.sig(uu), n)
R0_ufl = dot(T, ey) * ds(d34)
R0 = assemble_scalar(form(R0_ufl))

# Determination de K1 par interpolation du champ de déplacement
KI_interp = calcul_KI_via_interp_uy(domaine, modele, uu, a0)
# Determination de K1 par interpolation du champ de déplacement
pointe_fissure = [a0, 0]
G = 2 * calcul_G_via_G_theta(domaine, modele, uu, pointe_fissure, R_int, R_ext)
KI_G_theta = np.sqrt(G * modele.Ep)

# ----------- AFFICHAGE DES RÉSULTATS ---------------- #
print("=== Effort de reaction")
print(f"R0            : {R0:.3g}")
print("=== Facteur d'intensité des contraintes (mode I)")
print(f"KI interpolé  : {KI_interp:.3g}")
print(f"KI G et Irwin : {KI_G_theta:.3g}")
