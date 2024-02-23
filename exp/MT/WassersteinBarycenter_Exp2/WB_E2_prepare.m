% Example illustrating the Wasserstein barycenter problem with arbitrary
% mixture of bivariate Gaussians supported on a polytope

CONFIG = WB_E2_config();

test_num = 6;
marg_granularity_list = [5; 10; 15; 20; 25; 30] + 1;
quality_granularity_list = [10; 20; 30; 40; 50; 60] + 1;

% we will generate a total of 20 input measures:
% 10 input measures have 2 mixture components,
% 10 input measures have 3 mixture components.
% we will first generate the input measures on a fixed support, and then we
% will apply a random affine transformation

rng(2000, 'combRecursive');

marg_num = 20;

marg_weights = ones(marg_num, 1) / marg_num;

% the initial support
meas_supp_vertices = [1, 0; 0, -1; -1, 0; 0, 1];
meas_supp_triangles = [1, 2, 3; 3, 4, 1];

comp_num_list = [2; 2; 2; 2; 2; 2; 2; 2; 2; 2; ...
    3; 3; 3; 3; 3; 3; 3; 3; 3; 3];

marg_vertices_cell = cell(marg_num, 1);
marg_triangles_cell = cell(marg_num, 1);
marg_mixnorm_cell = cell(marg_num, 1);

rand_trans_cell = cell(marg_num, 1);
rand_shift_cell = cell(marg_num, 1);

for marg_id = 1:marg_num
    comp_num = comp_num_list(marg_id);
    mixnorm = struct;

    % randomly generate an affine transformation
    rand_shift = randn(2, 1) * 0.5;
    eigvals = gamrnd(10, 0.15, 2, 1);
    eigvec_angle = rand(1, 1) * 2 * pi;
    eigvec = [cos(eigvec_angle), -sin(eigvec_angle); ...
        sin(eigvec_angle), cos(eigvec_angle)];
    rand_trans = eigvec * diag(eigvals) * eigvec';

    rand_shift_cell{marg_id} = rand_shift;
    rand_trans_cell{marg_id} = rand_trans;

    marg_vertices_cell{marg_id} = meas_supp_vertices * rand_trans' ...
        + rand_shift';
    marg_triangles_cell{marg_id} = meas_supp_triangles;
    
    % the weights are Dirichlet(3, ..., 3) distributed
    mixnorm.weights = gamrnd(3, 1, comp_num, 1);
    mixnorm.weights = mixnorm.weights / sum(mixnorm.weights);
    
    mixnorm.components = cell(comp_num, 1);

    for comp_id = 1:comp_num
        mixnorm.components{comp_id} = struct;
        mixnorm.components{comp_id}.mean_vec = randn(2, 1) * 0.6;
        eigvals = 1 ./ gamrnd(3, 1, 2, 1);
        eigvec_angle = rand(1, 1) * 2 * pi;
        eigvec = [cos(eigvec_angle), -sin(eigvec_angle); ...
            sin(eigvec_angle), cos(eigvec_angle)];
        mixnorm.components{comp_id}.cov_mat = ...
            eigvec * diag(eigvals) * eigvec';

        % apply the transformation
        mixnorm.components{comp_id}.mean_vec = rand_trans ...
            * mixnorm.components{comp_id}.mean_vec + rand_shift;
        mixnorm.components{comp_id}.cov_mat = rand_trans ...
            * mixnorm.components{comp_id}.cov_mat * rand_trans';
    end

    marg_mixnorm_cell{marg_id} = mixnorm;
end

% create a grid for the test functions
marg_testfuncs_cell = cell(test_num, 1);

for test_id = 1:test_num
    testfuncs_x_num = marg_granularity_list(test_id);
    testfuncs_y_num = marg_granularity_list(test_id);

    meas_supp_testfuncs_x = linspace(-1, 1, testfuncs_x_num)';
    meas_supp_testfuncs_y = linspace(-1, 1, testfuncs_y_num)';
    [meas_supp_testfuncs_grid_x, meas_supp_testfuncs_grid_y] = ...
        meshgrid(meas_supp_testfuncs_x, meas_supp_testfuncs_y);
    meas_supp_testfuncs_vertices = [meas_supp_testfuncs_grid_x(:), ...
        meas_supp_testfuncs_grid_y(:)];
    meas_supp_testfuncs_vertices = meas_supp_testfuncs_vertices ...
        * [1/2, 1/2; -1/2, 1/2];
    meas_supp_testfuncs_vertices = unique([meas_supp_vertices; ...
        meas_supp_testfuncs_vertices], 'rows', 'stable');
    testfuncs_vertices_lower_list = ...
        find(meas_supp_testfuncs_vertices(:, 2) <= 1e-10);
    testfuncs_vertices_upper_list = ...
        find(meas_supp_testfuncs_vertices(:, 2) >= -1e-10);
    testfuncs_triangles_lower = delaunay(meas_supp_testfuncs_vertices( ...
        testfuncs_vertices_lower_list, :));
    testfuncs_triangles_lower = cleanup_triangles( ...
        meas_supp_testfuncs_vertices(testfuncs_vertices_lower_list, :), ...
        testfuncs_triangles_lower, [], false);
    testfuncs_triangles_upper = delaunay(meas_supp_testfuncs_vertices( ...
        testfuncs_vertices_upper_list, :));
    testfuncs_triangles_upper = cleanup_triangles( ...
        meas_supp_testfuncs_vertices(testfuncs_vertices_upper_list, :), ...
        testfuncs_triangles_upper, [], false);

    marg_testfuncs_cell{test_id} = cell(marg_num, 1);

    for marg_id = 1:marg_num
        marg_testfuncs_vertices = meas_supp_testfuncs_vertices ...
            * rand_trans_cell{marg_id}' + rand_shift_cell{marg_id}';
        marg_testfuncs_triangles = ...
            {testfuncs_vertices_lower_list(testfuncs_triangles_lower); ...
        testfuncs_vertices_upper_list(testfuncs_triangles_upper)};

        marg_testfuncs_cell{test_id}{marg_id} = ...
            {marg_testfuncs_vertices, marg_testfuncs_triangles, false};
    end
end

marg_cell = cell(marg_num, 1);

for marg_id = 1:marg_num
    marg_cell{marg_id} = ProbMeas2D_MixNorm( ...
        marg_vertices_cell{marg_id}, ...
        marg_triangles_cell{marg_id}, ...
        marg_mixnorm_cell{marg_id});
end

% construct the quality space from the convex hull of the union of the
% support of the input measures
union_support_vertices = vertcat(marg_vertices_cell{:});
convhull_indices = convhull(union_support_vertices);
quality_vertices = union_support_vertices(convhull_indices(1:end - 1), :);

WB_problem = MT2DQuad_ParTrans(marg_cell, marg_weights, ...
    {quality_vertices});

% construct the test functions on the quality space
quality_x_min = min(quality_vertices(:, 1));
quality_x_max = max(quality_vertices(:, 1));
quality_y_min = min(quality_vertices(:, 2));
quality_y_max = max(quality_vertices(:, 2));

quality_hp_w = WB_problem.Quality.Hyperplane.w;
quality_hp_steep_list = abs(quality_hp_w(:, 1)) > abs(quality_hp_w(:, 2));
quality_hp_b = WB_problem.Quality.Hyperplane.b;
quality_hp_num = length(quality_hp_b);

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
    
    quality_testfuncs_x = linspace(quality_x_min, quality_x_max, ...
        quality_testfuncs_x_num)';
    quality_testfuncs_y = linspace(quality_y_min, quality_y_max, ...
        quality_testfuncs_y_num)';
    [quality_testfuncs_grid_x, quality_testfuncs_grid_y] = ...
        meshgrid(quality_testfuncs_x, quality_testfuncs_y);
    quality_testfuncs_grid = [quality_testfuncs_grid_x(:), ...
        quality_testfuncs_grid_y(:)];

    % compute the intersections of the grid lines with the edges of the
    % polytope; if an edge is steeper than 45 degrees, then we compute the
    % intersections with horizontal grid lines; if an edge is no steeper
    % than 45 degrees, then we compute the intersections with vertical grid
    % lines
    quality_intersection_cell = cell(quality_hp_num, 1);

    for hp_id = 1:quality_hp_num
        if quality_hp_steep_list(hp_id)
            quality_intersection_cell{hp_id} = [(quality_hp_b(hp_id) ...
                - quality_hp_w(hp_id, 2) * quality_testfuncs_y) ...
                / quality_hp_w(hp_id, 1), quality_testfuncs_y];
        else
            quality_intersection_cell{hp_id} = [quality_testfuncs_x, ...
                (quality_hp_b(hp_id) - quality_hp_w(hp_id, 1) ...
                * quality_testfuncs_x) / quality_hp_w(hp_id, 2)];
        end
    end

    quality_testfuncs_vertices = unique([quality_vertices; ...
        vertcat(quality_intersection_cell{:}); ...
        quality_testfuncs_grid], 'rows', 'stable');
    quality_testfuncs_vertices_inside = ...
        WB_problem.checkIfInsideQualitySpace(quality_testfuncs_vertices);
    quality_testfuncs_vertices = ...
        quality_testfuncs_vertices(quality_testfuncs_vertices_inside, :);
    quality_testfuncs_triangles = delaunay(quality_testfuncs_vertices);
    [quality_testfuncs_triangles, quality_testfuncs_vertices_keep_list] ...
        = cleanup_triangles(quality_testfuncs_vertices, ...
        quality_testfuncs_triangles, [], true);
    quality_testfuncs_vertices = quality_testfuncs_vertices( ...
        quality_testfuncs_vertices_keep_list, :);

    quality_testfuncs_cell{test_id} = ...
        {quality_testfuncs_vertices, quality_testfuncs_triangles, false};
end

% grid used to plot the histogram 
quality_hist_x_num = 500;
quality_hist_y_num = 500;

quality_hist_edge_x = linspace(quality_x_min, quality_x_max, ...
    quality_hist_x_num + 1)';
quality_hist_edge_y = linspace(quality_y_min, quality_y_max, ...
    quality_hist_y_num + 1)';
quality_plot_hist_x = (quality_hist_edge_x(1:end - 1) ...
    + quality_hist_edge_x(2:end)) / 2;
quality_plot_hist_y = (quality_hist_edge_y(1:end - 1) ...
    + quality_hist_edge_y(2:end)) / 2;
[quality_plot_hist_grid_x, quality_plot_hist_grid_y] = ...
    meshgrid(quality_plot_hist_x, quality_plot_hist_y);
quality_plot_hist_grid = [quality_plot_hist_grid_x(:), ...
    quality_plot_hist_grid_y(:)];

save(CONFIG.SAVEPATH_INPUTS, ...
    'test_num', ...
    'marg_num', ...
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
 