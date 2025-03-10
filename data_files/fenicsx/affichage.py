import numpy as np
import pyvista as pv

import dolfinx.plot as plot


def affichage_fonction_scalaire_2D(domaine, ff):
    # Preparation des données
    cells, types, x = plot.vtk_mesh(domaine.maillage)
    grid = pv.UnstructuredGrid(cells, types, x)
    grid.point_data["u"] = ff.x.array

    # Affichage de ux
    pl = pv.Plotter()
    pl.add_mesh(grid, component=0, cmap="magma")
    pl.add_text(ff.name, color="k")
    pl.view_xy()
    pl.show()


def affichage_fonction_vectorielle_2D(domaine, uu):
    # Preparation des données
    cells, types, x = plot.vtk_mesh(domaine.maillage)
    grid = pv.UnstructuredGrid(cells, types, x)

    # Creation d'un masque pour récupérer les composantes
    x_component = [i for i in range(0, uu.x.array.shape[0], 2)]
    y_component = [i for i in range(1, uu.x.array.shape[0] + 1, 2)]

    # Stockage des composantes dans la grille
    grid.point_data["ux"] = uu.x.array[x_component]
    grid.point_data["uy"] = uu.x.array[y_component]

    # Liste des composantes
    composantes = ["ux", "uy"]
    indice_composante = 0

    # Creation du plot pyvista
    pl = pv.Plotter()
    pl.disable_shadows()

    # Fonction pour mettre à jour le tracé changement de composante
    def update_plot(component):
        composante = composantes[indice_composante]
        pl.clear()
        pl.add_mesh(
            grid,
            scalars=composante,
            scalar_bar_args={"title": f"Déplacement {composante}"},
            cmap="magma",
        )
        pl.view_xy()
        pl.reset_camera()
        pl.add_text(
            "Appuyer sur Espace pour changer de composante.",
            position="upper_right",
            font_size=12,
            color="black",
            name="instruction_text",
        )
        pl.render()

    # Fonction pour changer la composante (déclenchée par appui sur Espace)
    def cycle_components():
        nonlocal indice_composante
        indice_composante = (indice_composante + 1) % len(composantes)
        update_plot(indice_composante)

    # Ajouter un événement d'appui sur la touche Espace au tracé
    pl.add_key_event("space", cycle_components)
    # Initialisation du tracé
    update_plot(composantes[indice_composante])
    pl.show()


def affichage_deformee_2D(domaine, uu, amplitude):
    # Preparation des données
    cells, types, x = plot.vtk_mesh(domaine.maillage)
    grid = pv.UnstructuredGrid(cells, types, x)
    uu2D = uu.x.array.reshape(x.shape[0], uu.function_space.dofmap.index_map_bs)
    uu3D = np.hstack((uu2D, np.zeros((uu2D.shape[0], 1))))
    grid.point_data["u"] = uu3D

    # Affichage de la déformée
    deformed_grid = grid.copy()
    deformed_grid = grid.warp_by_vector(vectors="u", factor=amplitude)
    pl = pv.Plotter()
    pl.add_mesh(grid, style="wireframe", color="blue", opacity=0.25, line_width=2)
    pl.add_mesh(deformed_grid, style="wireframe", color="red", line_width=2)
    pl.add_text("Déformée", color="k")
    pl.view_xy()
    pl.show()


def affichage_contour_deformee_2D(domaine, uu, amplitude):
    # Preparation des données
    cells, types, x = plot.vtk_mesh(domaine.maillage)
    grid = pv.UnstructuredGrid(cells, types, x)
    uu2D = uu.x.array.reshape(x.shape[0], uu.function_space.dofmap.index_map_bs)
    uu3D = np.hstack((uu2D, np.zeros((uu2D.shape[0], 1))))
    grid.point_data["u"] = uu3D
    grid_edges = grid.extract_feature_edges()

    # Affichage de la déformée
    deformed_grid = grid.warp_by_vector(vectors="u", factor=amplitude)
    deformed_grid_edges = deformed_grid.extract_feature_edges()

    pl = pv.Plotter()
    pl.add_mesh(grid_edges, style="wireframe", color="blue", line_width=2)
    pl.add_mesh(deformed_grid_edges, style="wireframe", color="red", line_width=2)
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

    # Create a PyVista plotter
    pl = pv.Plotter()
    pl.disable_shadows()

    # List of stress components to cycle through
    components = ["sigxx", "sigxy", "sigyx", "sigyy"]
    component_index = 0

    # Function to update the plot based on the selected stress component
    def update_plot():
        nonlocal component_index
        component = components[component_index]
        pl.clear()
        pl.add_mesh(
            grid,
            scalars=component,
            scalar_bar_args={"title": f"Contrainte {component}"},
            cmap="magma",
        )
        pl.view_xy()
        pl.reset_camera()
        pl.add_text(
            "Appuyer sur Espace pour changer de composante.",
            position="upper_right",
            font_size=12,
            color="black",
            name="instruction_text",
        )
        pl.add_text(
            "Clic droit sur une cellule pour voir le tenseur de contrainte.",
            position="lower_left",
            font_size=12,
            color="black",
            name="value_text",
        )
        pl.render()

    # Define a callback function for mouse clicks
    def callback(point):
        # Find the closest cell to the clicked point
        cell_id = grid.find_closest_cell(point)
        cell_sigxx = grid["sigxx"][cell_id]
        cell_sigxy = grid["sigxy"][cell_id]
        cell_sigyy = grid["sigyy"][cell_id]

        # Format the probe text with consistent number formatting
        probe_text = (
            f"Point ({point[0]:.3g}, {point[1]:.3g}):\n"
            f"sig_xx = {cell_sigxx:.3g}\n"
            f"sig_xy = {cell_sigxy:.3g}\n"
            f"sig_yy = {cell_sigyy:.3g}"
        )

        # Update the text on the plot
        pl.remove_actor("value_text")  # Remove previous text actor
        pl.add_text(
            probe_text,
            position="lower_left",
            font_size=12,
            color="black",
            name="value_text",
        )
        pl.render()

    # Add the callback to the plotter
    pl.enable_point_picking(callback=callback, show_message=False)

    # Function to handle space bar press events
    def cycle_components():
        nonlocal component_index
        component_index = (component_index + 1) % len(components)
        update_plot()

    # Add space bar press event to the plotter
    pl.add_key_event("space", cycle_components)

    # Initialize the plot with the default component
    update_plot()

    pl.show()
