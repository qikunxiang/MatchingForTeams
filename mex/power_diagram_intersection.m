%% cell_indices = power_diagram_intersection(circles_cell, bbox_bottom_left, bbox_top_right)
% /// This function takes multiple collections of circles as well as a bounding box and computes the intersections of the cells in the corresponding 
% /// power diagrams formed by these collections of circles. It returns a list of indices where the intersections of the corresponding power cells
% /// have a non-empty intersection with the bounding box.
% /// Reference: Altschuler, J.M. and Boix-Adsera, E., 2021. Wasserstein barycenters can be computed in polynomial time in fixed dimension. 
% /// Journal of Machine Learning Research, 22(44), pp.1-19. See https://github.com/eboix/high_precision_barycenters/tree/master.
% /// Inputs:
% ///     circles_cell: cell array with two columns where the first column contains two-column matrices representing the centers of the circles and 
% ///     the second column contains vectors representing the squared radii of the circles
% ///     bbox_bottom_left: two-element vector representing the coordinate of the bottom-left corner of the bounding box
% ///     bbox_top_right: two-element vector representing the coordinate of the top-right corner of the bounding box
% /// Output:
% ///     cell_indices: a matrix where each row corresponds to a combination of cell indices which results in a non-empty intersection and each column
% ///     corresponds to a collection of circles