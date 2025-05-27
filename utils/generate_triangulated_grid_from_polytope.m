function [points, ...
    triangulation] ...
    = generate_triangulated_grid_from_polytope( ...
    vertices,  ...
    grid_x, ...
    grid_y)
% Function that takes a polytope and an existing grid that covers the polytope as inputs and returns a set of new grid points and a
% triangulation of them to cover the polytope exactly with triangles. Each triangle will be a subset of a rectangle in the input grid.
% Inputs: 
%   vertices: two-column matrix containing the vertices of the polytope
%   grid_x: vector containing the x-coordinates of the grid points
%   grid_y: vector containing the y-coordinates of the grid points
% Output:
%   points: two-column matrix containing the vertices in the triangulation
%   triangulation: three-column matrix containing the vertex indices of the triangles

assert(min(grid_x) <= min(vertices(:, 1)) && max(grid_x) >= max(vertices(:, 1)) ...
    && min(grid_y) <= min(vertices(:, 2)) && max(grid_y) >= max(vertices(:, 2)), 'the grid needs to cover the polytope');

grid_x_num = length(grid_x);
grid_y_num = length(grid_y);

ps = polyshape(vertices(:, 1), vertices(:, 2), 'Simplify', false, 'KeepCollinearPoints', true);

[grid_x_mat, grid_y_mat] = meshgrid(grid_x, grid_y);

% identify all the grid points contained inside the polytope
inside_mat = reshape(isinterior(ps, grid_x_mat(:), grid_y_mat(:)), grid_y_num, grid_x_num)';
[inside_x_list, inside_y_list] = find(inside_mat);

% add the grid points inside the polytope to the output triangulation
% reserve enough space for the points in the output triangulation
points = zeros(length(inside_x_list) * 4 + 100, 2);
points(1:length(inside_x_list), :) = [grid_x(inside_x_list), grid_y(inside_y_list)];
point_counter = length(inside_x_list);

inside_mat_bl = inside_mat(1:end - 1, 1:end - 1);
inside_mat_br = inside_mat(2:end, 1:end - 1);
inside_mat_tl = inside_mat(1:end - 1, 2:end);
inside_mat_tr = inside_mat(2:end, 2:end);

% identify all the squares that are partially overlapping with the polytope (have non-empty intersection but not fully contained)
[partial_x_list, partial_y_list] = find((inside_mat_bl | inside_mat_br | inside_mat_tl | inside_mat_tr) ...
    & ~(inside_mat_bl & inside_mat_br & inside_mat_tl & inside_mat_tr));

for rect_id = 1:length(partial_x_list)
    x_id = partial_x_list(rect_id);
    y_id = partial_y_list(rect_id);
    rect = polyshape([grid_x(x_id); grid_x(x_id); grid_x(x_id + 1); grid_x(x_id + 1)], ...
        [grid_y(y_id); grid_y(y_id + 1); grid_y(y_id + 1); grid_y(y_id)]);

    % if the intersection is non-empty, add the vertices of the polytope formed by the intersection to the output triangulation
    intersection = intersect(ps, rect);

    vert_num = size(intersection.Vertices, 1);
    points(point_counter + (1:vert_num), :) = intersection.Vertices;
    point_counter = point_counter + vert_num;
end

% add the vertices in the input polytope to the top and remove (approximate) duplicates
points = [vertices; points(1:point_counter, :)];
[~, uind] = unique(round(points, 4), 'rows', 'stable');
points = points(uind, :);

triangulation = delaunay(points(:, 1), points(:, 2));

end