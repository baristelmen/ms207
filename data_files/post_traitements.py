from dolfinx.geometry import bb_tree, compute_collisions_points, compute_colliding_cells


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
