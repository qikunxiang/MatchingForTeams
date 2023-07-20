function [empty, is_simpcover, condition_mat] ...
    = check_triangulation_simplicial_cover(vertices, triangles)
% Function that determines whether all triangles in a triangulation are
% non-empty and satisfy the condition of being a simplicial cover. To test
% for illegal intersections (i.e., the intersection of two triangles is not
% a face of one of the triangles), trianges will be expanded by a factor of
% 1e-6. Thus, two interior-disjoint triangles that are too close to each 
% other are considered illegal as well. 
% Inputs: 
%   vertices: two-clumn matrix containing the vertices of all triangles
%   triangles: three-column matrix containing the vertex indices in the
%   triangulation
% Output:
%   empty: boolean vector indicating whether the triangles are empty
%   is_simpcover: boolean indicating whether the triangles form a
%   simplicial cover, that is, any non-empty intersection of two triangles
%   is a face of both triangles
%   condition_mat: matrix indicating whether each pair of triangles
%   satisfies the condition

tri_num = size(triangles, 1);

poly_list = repmat(polyshape, tri_num, 1);

for tri_id = 1:tri_num
    poly_list(tri_id) = polyshape(vertices(triangles(tri_id, :), 1), ...
        vertices(triangles(tri_id, :), 2));
end

empty = area(poly_list) < 1e-10;

condition_mat = true(tri_num, tri_num);

expand = 1e-6;
expand_mat = [1 + 2 * expand, -expand, -expand;
              -expand, 1 + 2 * expand, -expand;
              -expand, -expand, 1 + 2 * expand];
edge_expand_mat = [1, 0, 0;
                   0, 1 + expand, -expand;
                   0, -expand, 1 + expand];

% iterate through pairs of triangles
for tri_id1 = 1:tri_num
    for tri_id2 = tri_id1 + 1:tri_num
        % count the number of distince vertices in the two triangles
        unique_vert_num = length(unique([triangles(tri_id1, :), ...
            triangles(tri_id2, :)]));
        
        if unique_vert_num == 6
            % case 1: no shared vertex, thus the condition is satisfied if
            % and only if the interior of the two triangles are still
            % disjoint after expanding both triangles slightly
            expanded_tri1_v = expand_mat ...
                * vertices(triangles(tri_id1, :), :);
            expanded_poly1 = polyshape(expanded_tri1_v(:, 1), ...
                expanded_tri1_v(:, 2));
            expanded_tri2_v = expand_mat ...
                * vertices(triangles(tri_id2, :), :);
            expanded_poly2 = polyshape(expanded_tri2_v(:, 1), ...
                expanded_tri2_v(:, 2));
            condition_mat(tri_id1, tri_id2) = ~overlaps(expanded_poly1, ...
                expanded_poly2);

        elseif unique_vert_num == 5
            % case 2: the two triangles share one vertex, thus the
            % condition is satisfied if and only if the interior of the two
            % triangles are still disjoint after expanding the opposite
            % edges of the two trianges slightly

            % find the shared vertex
            shared_v = intersect(triangles(tri_id1, :), ...
                triangles(tri_id2, :));
            
            % reorder the vertices such that the shared vertex comes first
            tri1_reordered = [shared_v, setdiff(triangles(tri_id1, :), ...
                shared_v)];
            tri2_reordered = [shared_v, setdiff(triangles(tri_id2, :), ...
                shared_v)];
            expanded_tri1_v = edge_expand_mat ...
                * vertices(tri1_reordered, :);
            expanded_poly1 = polyshape(expanded_tri1_v(:, 1), ...
                expanded_tri1_v(:, 2));
            expanded_tri2_v = edge_expand_mat ...
                * vertices(tri2_reordered, :);
            expanded_poly2 = polyshape(expanded_tri2_v(:, 1), ...
                expanded_tri2_v(:, 2));
            condition_mat(tri_id1, tri_id2) = ~overlaps(expanded_poly1, ...
                expanded_poly2);

        elseif unique_vert_num == 4
            % case 3: the two triangles share an edge, thus the condition
            % is satisfied if and only if the interior of the two triangles
            % are disjoint
            condition_mat(tri_id1, tri_id2) = ~overlaps( ...
                poly_list(tri_id1), poly_list(tri_id2));
        else
            % case 4: two triangles are identical, hence illegal
            condition_mat(tri_id1, tri_id2) = false;
        end
    end
end

condition_mat = condition_mat & condition_mat';

is_simpcover = all(all(condition_mat));

end