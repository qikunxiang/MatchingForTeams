function [T_clean, P_keep_list] = cleanup_triangles(P, T, thres, ...
    remove_redundant_vertices)
% Given a set of points and a triangulation, clean up triangles that are close to being degenerate
% Inputs:
%   P: matrix where each row represents a point
%   T: matrix where each row contains point indices of a triangle
%   thres: threshold for the condition number (default is 1e12)
%   remove_redundant_vertices: logical variable indicating whether to remove redundant vertices from the original triangulation 
%   (default is false)
% Output:
%   T_clean: matrix where the almost degenerate triangles are removed
%   P_keep_list: matrix where the redundant vertices after cleaning up are removed

if ~exist('thres', 'var') || isempty(thres)
    thres = 1e12;
end

if ~exist('remove_redundant_vertices', 'var') || isempty(remove_redundant_vertices)
    remove_redundant_vertices = false;
end

tri_num = size(T, 1);
dim = size(P, 2);

keep_list = true(tri_num, 1);

for tri_id = 1:tri_num
    tri_cond = cond([P(T(tri_id, :), :)'; ones(1, dim + 1)]);

    keep_list(tri_id) = tri_cond <= thres;
end

T_clean = T(keep_list, :);
T_clean_num = size(T_clean, 1);

if remove_redundant_vertices
    % remove redundant vertices
    [P_keep_list, ~, T_clean_unique_mapping] = unique(T_clean(:), 'sorted');
    T_clean = reshape(T_clean_unique_mapping, T_clean_num, 3);
else
    P_keep_list = true(size(P, 1), 1);
end

end