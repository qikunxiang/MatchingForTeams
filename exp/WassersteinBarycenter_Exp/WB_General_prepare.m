% Example computing the Wasserstein barycenter of probability measures supported on the union of two triangles with continuous 
% piece-wise affine density functions

CONFIG = WB_General_config();

test_num = 6;
marg_granularity_list = [1; 2; 4; 8; 16; 32] + 1;
quality_granularity_list = [4; 8; 16; 32; 64; 128] + 1;

% we will generate a total of 20 input measures

rng(2000, 'combRecursive');

marg_num = 20;

marg_weights = ones(marg_num, 1) / marg_num;

marg_vertices_pretransform_cell = cell(marg_num, 1);
marg_vertices_cell = cell(marg_num, 1);
marg_triangles_cell = cell(marg_num, 1);
marg_densities_cell = cell(marg_num, 1);

rand_trans_cell = cell(marg_num, 1);
rand_shift_cell = cell(marg_num, 1);

for marg_id = 1:marg_num
    % first randomly generate two triangles before the affine transformation
    angle_upper = rand * pi * 2/3 + pi * 1/6;
    angle_lower = rand * pi * 2/3 + pi * 1/6 + pi;

    vertex_center_left = [-1/2; 0];
    vertex_center_right = [1/2; 0];
    vertex_center_center = [0; 0];
    vertex_upper = [cos(angle_upper); sin(angle_upper)];
    vertex_lower = [cos(angle_lower); sin(angle_lower)];
    vertex_mid_upper_left = (vertex_upper + vertex_center_left) / 2;
    vertex_mid_upper_right = (vertex_upper + vertex_center_right) / 2;
    vertex_mid_lower_left = (vertex_lower + vertex_center_left) / 2;
    vertex_mid_lower_right = (vertex_lower + vertex_center_right) / 2;

    marg_vertices_pretransform = [vertex_center_left'; vertex_center_right'; vertex_upper'; vertex_lower'; vertex_center_center'; ...
        vertex_mid_upper_left'; vertex_mid_upper_right'; vertex_mid_lower_left'; vertex_mid_lower_right'];
    marg_triangles = [ ...
        1, 5, 6; ...
        2, 5, 7; ...
        5, 6, 7; ...
        3, 6, 7; ...
        1, 5, 8; ...
        2, 5, 9; ...
        5, 8, 9; ...
        4, 8, 9
    ];

    % randomly generate an affine transformation
    rand_shift = randn(2, 1) * 0.5;
    eigvals = gamrnd(10, 0.15, 2, 1);
    eigvec_angle = rand(1, 1) * 2 * pi;
    eigvec = [cos(eigvec_angle), -sin(eigvec_angle); sin(eigvec_angle), cos(eigvec_angle)];
    rand_trans = eigvec * diag(eigvals) * eigvec';

    rand_shift_cell{marg_id} = rand_shift;
    rand_trans_cell{marg_id} = rand_trans;

    % randomly generate the densities (without normalization)
    while true
        marg_densities = randi([0, 2], 9, 1);

        if any(marg_densities > 0)
            break;
        end
    end

    % remove triangles with 0 probabilities
    nonzero_list = any(marg_densities(marg_triangles')' > 0, 2);
    marg_triangles = marg_triangles(nonzero_list, :);
    marg_triangles_unwrapped = marg_triangles(:);
    [uverts, ~, umap] = unique(marg_triangles_unwrapped);
    marg_vertices_pretransform = marg_vertices_pretransform(uverts, :);
    marg_densities = marg_densities(uverts);
    marg_triangles = reshape(umap, size(marg_triangles, 1), 3);

    marg_vertices_pretransform_cell{marg_id} = marg_vertices_pretransform;
    marg_vertices_cell{marg_id} = marg_vertices_pretransform * rand_trans' + rand_shift';
    marg_triangles_cell{marg_id} = marg_triangles;
    marg_densities_cell{marg_id} = marg_densities;
end

% create a grid for the test functions
marg_testfuncs_cell = cell(test_num, 1);

for test_id = 1:test_num
    testfuncs_x_num = marg_granularity_list(test_id);
    testfuncs_y_num = marg_granularity_list(test_id);

    % build the test functions in a triangle: [0, 0; 1, 0; 0, 1]
    testfuncs_x_grid = linspace(0, 1, testfuncs_x_num);
    testfuncs_y_grid = linspace(0, 1, testfuncs_y_num);

    [G_x, G_y] = meshgrid(testfuncs_x_grid, testfuncs_y_grid);
    canonical_grid_points = [G_x(:), G_y(:)];
    canonical_grid_points = canonical_grid_points(sum(canonical_grid_points, 2) < 1 + 1e-6, :);
    canonical_grid_point_num = size(canonical_grid_points, 1);
    canonical_grid_triangulation = delaunay(canonical_grid_points(:, 1), canonical_grid_points(:, 2));
    canonical_grid_triangulation = cleanup_triangles(canonical_grid_points, canonical_grid_triangulation, 1e12, false);

    marg_testfuncs_cell{test_id} = cell(marg_num, 1);

    for marg_id = 1:marg_num
        tri_num = size(marg_triangles_cell{marg_id}, 1);
        marg_testfuncs_vertices_cell = cell(tri_num, 1);

        for tri_id = 1:tri_num
            % compute the affine transformation from the triangle [0, 0; 1, 0; 0, 1] to the target triangle in the density mesh
            tri_vertices = marg_vertices_pretransform_cell{marg_id}(marg_triangles_cell{marg_id}(tri_id, :), :);
            trans_mat = [tri_vertices'; ones(1, 3)] / [0, 1, 0; 0, 0, 1; 1, 1, 1];

            marg_testfuncs_vertices_cell{tri_id} = [canonical_grid_points, ones(canonical_grid_point_num, 1)] * trans_mat(1:2, :)';
        end

        marg_testfuncs_vertices_agg = [marg_vertices_pretransform_cell{marg_id}; vertcat(marg_testfuncs_vertices_cell{:})];
        [~, uind, umap] = unique(round(marg_testfuncs_vertices_agg, 6), 'rows', 'stable');
        marg_testfuncs_vertices_unique = marg_testfuncs_vertices_agg(uind, :);

        marg_testfuncs_triangles_cell = cell(tri_num, 1);
        counter = size(marg_vertices_pretransform_cell{marg_id}, 1);

        for tri_id = 1:tri_num
            umap_tri = umap(counter + (1:canonical_grid_point_num), :);
    
            marg_testfuncs_triangles_cell{tri_id} = umap_tri(canonical_grid_triangulation')';
    
            counter = counter + canonical_grid_point_num;
        end

        marg_testfuncs_vertices = marg_testfuncs_vertices_unique * rand_trans_cell{marg_id}' + rand_shift_cell{marg_id}';
        marg_testfuncs_cell{test_id}{marg_id} = {marg_testfuncs_vertices, marg_testfuncs_triangles_cell, false};
    end
end

% construct the test functions on the quality space
weighted_marg_vertices_cell = cell(marg_num, 1);

for marg_id = 1:marg_num
    weighted_marg_vertices_cell{marg_id} = marg_vertices_cell{marg_id} * marg_weights(marg_id);
end

quality_vertices = minkowski_sum_2d(weighted_marg_vertices_cell);

quality_x_min = min(quality_vertices(:, 1));
quality_x_max = max(quality_vertices(:, 1));
quality_y_min = min(quality_vertices(:, 2));
quality_y_max = max(quality_vertices(:, 2));

quality_testfuncs_cell = cell(test_num, 1);

for test_id = 1:test_num
    quality_testfuncs_x_num = quality_granularity_list(test_id);
    quality_testfuncs_y_num = quality_granularity_list(test_id);

    if (quality_x_max - quality_x_min) > (quality_y_max - quality_y_min)
        quality_testfuncs_y_num = ceil(quality_testfuncs_x_num ...
            / (quality_x_max - quality_x_min) ...
            * (quality_y_max - quality_y_min));
    else
        quality_testfuncs_x_num = ceil(quality_testfuncs_y_num ...
            / (quality_y_max - quality_y_min) ...
            * (quality_x_max - quality_x_min));
    end

    quality_testfuncs_x = linspace(quality_x_min - 1e-3, quality_x_max + 1e-3, ...
        quality_testfuncs_x_num)';
    quality_testfuncs_y = linspace(quality_y_min - 1e-3, quality_y_max + 1e-3, ...
        quality_testfuncs_y_num)';
    
    [quality_testfuncs_vertices, quality_testfuncs_triangles] ...
        = generate_triangulated_grid_from_polytope(quality_vertices, quality_testfuncs_x, quality_testfuncs_y);

    [quality_testfuncs_triangles, vertex_keep_list] = cleanup_triangles( ...
        quality_testfuncs_vertices, quality_testfuncs_triangles, 1e12, true);
    
    quality_testfuncs_cell{test_id} = ...
        {quality_testfuncs_vertices(vertex_keep_list, :), quality_testfuncs_triangles, true, false};
end

% grid used to plot the histogram 
quality_hist_x_num = 300;
quality_hist_y_num = 300;

quality_hist_edge_x = linspace(quality_x_min, quality_x_max, quality_hist_x_num + 1)';
quality_hist_edge_y = linspace(quality_y_min, quality_y_max, quality_hist_y_num + 1)';
quality_plot_hist_x = (quality_hist_edge_x(1:end - 1) + quality_hist_edge_x(2:end)) / 2;
quality_plot_hist_y = (quality_hist_edge_y(1:end - 1) + quality_hist_edge_y(2:end)) / 2;
[quality_plot_hist_grid_x, quality_plot_hist_grid_y] = meshgrid(quality_plot_hist_x, quality_plot_hist_y);
quality_plot_hist_grid = [quality_plot_hist_grid_x(:), quality_plot_hist_grid_y(:)];

save(CONFIG.SAVEPATH_INPUTS, ...
    'test_num', ...
    'marg_num', ...
    'marg_weights', ...
    'marg_vertices_cell', ...
    'marg_triangles_cell', ...
    'marg_densities_cell', ...
    'marg_testfuncs_cell', ...
    'quality_vertices', ...
    'quality_testfuncs_cell', ...
    'quality_x_min', ...
    'quality_x_max', ...
    'quality_y_min', ...
    'quality_y_max', ...
    'quality_hist_x_num', ...
    'quality_hist_y_num', ...
    'quality_hist_edge_x', ...
    'quality_hist_edge_y', ...
    'quality_plot_hist_grid_x', ...
    'quality_plot_hist_grid_y', ...
    'quality_plot_hist_grid', ...
    '-v7.3');
 