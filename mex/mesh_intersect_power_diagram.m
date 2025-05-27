%% [output_vertices, output_triangulation, power_cell_indices, mesh_indices] = mesh_intersect_power_diagram(mesh_vertices, mesh_triangulation, circle_centers, circle_radii, bbox_bottom_left, bbox_top_right)
% /// This function takes a triangular mesh, a collection of circles, and a bounding box and computes the intersections of the triangles in the mesh with 
% /// the cells in the power diagram. It returns triangles formed via triangulating the resulting intersections as well as which mesh triangle and which power
% /// cell each of the output triangles belongs to. 
% /// Reference: Altschuler, J.M. and Boix-Adserà, E., 2021. Wasserstein barycenters can be computed in polynomial time in fixed dimension. Journal of Machine Learning Research, 22(44), pp.1-19. See https://github.com/eboix/high_precision_barycenters/tree/master.
% /// Inputs:
% ///     mesh_vertices: two-column matrix containing the vertices of the triangular mesh
% ///     mesh_triangulation: three-column matrix representing the mesh triangulation
% ///     circle_centers: two-column matrix containing the centers of the circles
% ///     circle_sqradii: vector containing the squared radii of the circles
% ///     bbox_bottom_left: two-element vector representing the coordinate of the bottom-left corner of the bounding box; the bounding box must contain the entire triangular mesh
% ///     bbox_top_right: two-element vector representing the coordinate of the top-right corner of the bounding box; the bounding box must contain the entire triangular mesh
% /// Outputs:
% ///     output_vertices: two-column matrix containing the vertices in the output triangular mesh
% ///     output_triangulation: three-column matrix representing the output triangulation
% ///     power_cell_indices: vector containing the power cell indices of the output triangles
% ///     mesh_indices: vector containing the mesh triangle indices of the output triangles