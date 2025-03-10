from mpi4py import MPI
from dolfinx.io.gmshio import model_to_mesh
import dolfinx.plot as plot
import pyvista as pv


class Domaine:
    def __init__(self, gmsh_model):
        dim = gmsh_model.get_dimension()
        self.maillage, self.cell_markers, self.facet_markers = model_to_mesh(
            gmsh_model, MPI.COMM_WORLD, 0, gdim=dim
        )

    def affichage(self):
        cells, types, x = plot.vtk_mesh(self.maillage)
        grid = pv.UnstructuredGrid(cells, types, x)

        plotter = pv.Plotter()
        plotter.add_mesh(grid, show_edges=True)
        plotter.view_xy()
        plotter.show()
