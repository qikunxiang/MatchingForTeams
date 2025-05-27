% Example illustrating the Wasserstein barycenter problem when all measures belong to the same elliptical family

CONFIG = WB_Elliptical_config();

test_num = 4;
marg_granularity_list = [1; 2; 4; 8] + 1;
quality_granularity_list = [16; 32; 64; 128] + 1;

rng(1000, 'combRecursive');

% the initial measure used to generate the measures from an elliptical family via affine transformations
% the support is a regular 36-gon which approximates a circle, since this implementation restricts the support to the union of
% triangles
meas_init_vertices = [ ...
    0, 0;
    1/3 * [cos((0:5)' / 6 * 2 * pi), sin((0:5)' / 6 * 2 * pi)];
    2/3 * [cos((0:17)' / 18 * 2 * pi), sin((0:17)' / 18 * 2 * pi)];
    3/3 * [cos((0:35)' / 36 * 2 * pi), sin((0:35)' / 36 * 2 * pi)];
];
meas_init_triangles = delaunay(meas_init_vertices);

% the measure has three equally weighted mixture components
comp_num = 2;
meas_init_mixnorm = struct;
meas_init_mixnorm.weights = [1/4; 3/4];
meas_init_mixnorm.components = cell(comp_num, 1);

meas_init_mixnorm.components{1} = struct;
meas_init_mixnorm.components{1}.mean_vec = [0; 0];
meas_init_mixnorm.components{1}.cov_mat = [0.04, 0; 0, 0.04];

meas_init_mixnorm.components{2} = struct;
meas_init_mixnorm.components{2}.mean_vec = [0; 0];
meas_init_mixnorm.components{2}.cov_mat = [1, 0; 0, 1];

% create a grid for the test functions
meas_init_testfuncs_cell = cell(test_num, 1);

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

    tri_num = size(meas_init_triangles, 1);
    marg_init_testfuncs_vertices_cell = cell(tri_num, 1);

    for tri_id = 1:tri_num
        % compute the affine transformation from the triangle [0, 0; 1, 0; 0, 1] to the target triangle in the density mesh
        tri_vertices = meas_init_vertices(meas_init_triangles(tri_id, :), :);
        trans_mat = [tri_vertices'; ones(1, 3)] / [0, 1, 0; 0, 0, 1; 1, 1, 1];

        marg_init_testfuncs_vertices_cell{tri_id} = [canonical_grid_points, ones(canonical_grid_point_num, 1)] * trans_mat(1:2, :)';
    end

    marg_init_testfuncs_vertices_agg = [meas_init_vertices; vertcat(marg_init_testfuncs_vertices_cell{:})];
    [~, uind, umap] = unique(round(marg_init_testfuncs_vertices_agg * 1024, 3), 'rows', 'stable');
    marg_init_testfuncs_vertices_unique = marg_init_testfuncs_vertices_agg(uind, :);

    marg_init_testfuncs_triangles_cell = cell(tri_num, 1);
    counter = size(meas_init_vertices, 1);

    for tri_id = 1:tri_num
        umap_tri = umap(counter + (1:canonical_grid_point_num), :);

        marg_init_testfuncs_triangles_cell{tri_id} = umap_tri(canonical_grid_triangulation')';

        counter = counter + canonical_grid_point_num;
    end

    meas_init_testfuncs_cell{test_id} = {marg_init_testfuncs_vertices_unique, marg_init_testfuncs_triangles_cell, false};
end


Meas_init = ProbMeas2DMixNorm(meas_init_vertices, meas_init_triangles, meas_init_mixnorm);

meas_init_meanvec = Meas_init.meanVector();
meas_init_covmat = Meas_init.covarianceMatrix();
meas_init_trans = inv(chol(meas_init_covmat))';
meas_init_shift = -meas_init_trans * meas_init_meanvec;

% create the reference measure which is the affine transformation of the initial measure which has zero mean and unit covariance
meas_ref_vertices = meas_init_vertices * meas_init_trans' + meas_init_shift';
meas_ref_triangles = meas_init_triangles;

meas_ref_mixnorm = struct;
meas_ref_mixnorm.weights = meas_init_mixnorm.weights;
meas_ref_mixnorm.components = cell(comp_num, 1);

for comp_id = 1:comp_num
    meas_ref_mixnorm.components{comp_id} = struct;
    meas_ref_mixnorm.components{comp_id}.mean_vec = meas_init_trans * meas_init_mixnorm.components{comp_id}.mean_vec + meas_init_shift;
    meas_ref_mixnorm.components{comp_id}.cov_mat = meas_init_trans * meas_init_mixnorm.components{comp_id}.cov_mat * meas_init_trans';
end

% grid for the test functions
meas_ref_testfuncs_cell = cell(test_id, 1);

for test_id = 1:test_num
    meas_ref_testfuncs_cell{test_id} = {meas_init_testfuncs_cell{test_id}{1} * meas_init_trans' + meas_init_shift', ...
        meas_init_testfuncs_cell{test_id}{2}, false};
end


% generate the marginals via shifting and shearing
marg_num = 5;
marg_weights = ones(marg_num, 1) / marg_num;

% randomly generate the mean vectors and covariance matrices
elliptical_cell = cell(marg_num, 2);

for marg_id = 1:marg_num
    elliptical_cell{marg_id, 1} = randn(2, 1) * 0.5;
    eigvals = gamrnd(10, 1/10, 2, 1);
    eigvec_angle = rand(1, 1) * 2 * pi;
    eigvec = [cos(eigvec_angle), -sin(eigvec_angle); sin(eigvec_angle), cos(eigvec_angle)];
    elliptical_cell{marg_id, 2} = eigvec * diag(eigvals) * eigvec';
end

marg_vertices_cell = cell(marg_num, 1);
marg_triangles_cell = cell(marg_num, 1);
marg_mixnorm_cell = cell(marg_num, 1);

marg_testfuncs_cell = cell(test_num, 1);

for test_id = 1:test_num
    marg_testfuncs_cell{test_id} = cell(marg_num, 1);
end

marg_cell = cell(marg_num, 1);

for marg_id = 1:marg_num
    trans_mat = elliptical_cell{marg_id, 2};
    shift_vec = elliptical_cell{marg_id, 1};
    marg_verts = meas_ref_vertices * trans_mat' + shift_vec';
    marg_triangles = meas_ref_triangles;

    marg_mixnorm = struct;
    marg_mixnorm.weights = meas_ref_mixnorm.weights;
    marg_mixnorm.components = cell(comp_num, 1);

    for comp_id = 1:comp_num
        marg_mixnorm.components{comp_id} = struct;
        marg_mixnorm.components{comp_id}.mean_vec = trans_mat * meas_ref_mixnorm.components{comp_id}.mean_vec + shift_vec;
        marg_mixnorm.components{comp_id}.cov_mat = trans_mat * meas_ref_mixnorm.components{comp_id}.cov_mat * trans_mat';
    end

    marg_vertices_cell{marg_id} = marg_verts;
    marg_triangles_cell{marg_id} = marg_triangles;
    marg_mixnorm_cell{marg_id} = marg_mixnorm;

    for test_id = 1:test_num
        marg_testfuncs_cell{test_id}{marg_id} = {meas_ref_testfuncs_cell{test_id}{1} * trans_mat' + shift_vec', ...
            meas_ref_testfuncs_cell{test_id}{2}, false};
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
        {quality_testfuncs_vertices(vertex_keep_list, :), quality_testfuncs_triangles, false, false};
end

% grid used to plot the histogram and density
quality_hist_x_num = 500;
quality_hist_y_num = 500;

quality_hist_edge_x = linspace(quality_x_min, quality_x_max, quality_hist_x_num + 1)';
quality_hist_edge_y = linspace(quality_y_min, quality_y_max, quality_hist_y_num + 1)';
quality_plot_hist_x = (quality_hist_edge_x(1:end - 1) + quality_hist_edge_x(2:end)) / 2;
quality_plot_hist_y = (quality_hist_edge_y(1:end - 1) + quality_hist_edge_y(2:end)) / 2;
[quality_plot_hist_grid_x, quality_plot_hist_grid_y] = meshgrid(quality_plot_hist_x, quality_plot_hist_y);
quality_plot_hist_grid = [quality_plot_hist_grid_x(:), quality_plot_hist_grid_y(:)];

save(CONFIG.SAVEPATH_INPUTS, ...
    'test_num', ...
    'marg_num', ...
    'meas_ref_vertices', ...
    'meas_ref_triangles', ...
    'meas_ref_mixnorm', ...
    'meas_ref_testfuncs_cell', ...
    'marg_weights', ...
    'marg_vertices_cell', ...
    'marg_triangles_cell', ...
    'marg_mixnorm_cell', ...
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