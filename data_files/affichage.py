import numpy as np
import pyvista as pv

import dolfinx.plot as plot


def affichage_fonction_vectorielle_2D(domaine, uu):
    # Preparation des données
    cells, types, x = plot.vtk_mesh(domaine.maillage)
    grid = pv.UnstructuredGrid(cells, types, x)
    uu2D = uu.x.array.reshape(x.shape[0], uu.function_space.dofmap.index_map_bs)
    grid.point_data["u"] = uu2D

    # Affichage de ux
    pl = pv.Plotter()
    pl.add_mesh(grid, component=0, cmap="magma")
    pl.add_text("Déplacement selon x", color="k")
    pl.view_xy()
    pl.show()

    # Affichage de uy
    pl = pv.Plotter()
    pl.add_mesh(grid, component=1, cmap="magma")
    pl.add_text("Déplacement selon y", color="k")
    pl.view_xy()
    pl.show()


def affichage_deformee_2D(domaine, uu):
    # Preparation des données
    cells, types, x = plot.vtk_mesh(domaine.maillage)
    grid = pv.UnstructuredGrid(cells, types, x)
    uu2D = uu.x.array.reshape(x.shape[0], uu.function_space.dofmap.index_map_bs)
    uu3D = np.hstack((uu2D, np.zeros((uu2D.shape[0], 1))))
    grid.point_data["u"] = uu3D

    # Affichage de la déformée
    deformed_grid = grid.copy()
    deformed_grid.warp_by_vector(vectors="u", factor=200, inplace=True)
    pl = pv.Plotter()
    pl.add_mesh(grid, style="wireframe", color="blue", opacity=0.25, line_width=2)
    pl.add_mesh(deformed_grid, style="wireframe", color="red", line_width=2)
    pl.add_text("Déformée", color="k")
    pl.view_xy()
    pl.show()


def affichage_fonction_tensorielle_2D(domaine, T):
    # Preparation des données
    cells, types, x = plot.vtk_mesh(domaine.maillage)
    grid = pv.UnstructuredGrid(cells, types, x)
    # Create mask to recover the components
    xx_component = [i for i in range(0, cells.shape[0], 4)]
    xy_component = [i for i in range(1, cells.shape[0] + 1, 4)]
    yx_component = [i for i in range(2, cells.shape[0] + 1, 4)]
    yy_component = [i for i in range(3, cells.shape[0] + 1, 4)]
    # Stockage des composantes dans la grille
    grid.cell_data["sigxx"] = T.x.array[xx_component]
    grid.cell_data["sigxy"] = T.x.array[xy_component]
    grid.cell_data["sigyx"] = T.x.array[yx_component]
    grid.cell_data["sigyy"] = T.x.array[yy_component]

    # Affichage du tenseur
    pl = pv.Plotter(shape=(2, 2))
    pl.disable_shadows()

    pl.subplot(0, 0)
    pl.add_mesh(
        grid,
        scalars="sigxx",
        scalar_bar_args={"title": "Contrainte sig_xx [Pa]"},
        cmap="magma",
    )
    pl.view_xy()

    pl.subplot(0, 1)
    pl.add_mesh(
        grid,
        scalars="sigxy",
        scalar_bar_args={"title": "Contrainte sig_xy [Pa]"},
        cmap="magma",
    )
    pl.view_xy()

    pl.subplot(1, 0)
    pl.add_mesh(
        grid,
        scalars="sigyx",
        scalar_bar_args={"title": "Contrainte sig_yx [Pa]"},
        cmap="magma",
    )
    pl.view_xy()

    pl.subplot(1, 1)
    pl.add_mesh(
        grid,
        scalars="sigyy",
        scalar_bar_args={"title": "Contrainte sig_yy [Pa]"},
        cmap="magma",
    )
    pl.view_xy()

    pl.show()
