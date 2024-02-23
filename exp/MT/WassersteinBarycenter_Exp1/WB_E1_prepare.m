% Example illustrating the Wasserstein barycenter problem when all measures
% belong to the same location-scatter family

CONFIG = WB_E1_config();

test_num = 6;
marg_granularity_list = [4; 12; 20; 28; 36; 44] + 1;
quality_granularity_list = [10; 30; 50; 70; 90; 110] + 1;

rng(1000, 'combRecursive');

% the initial measure used to generate the measures from a
% location-scatter family via affine transformations
meas_init_vertices = [1, 0; 0, -1; -1, 0; 0, 1];
meas_init_triangles = [1, 2, 3; 3, 4, 1];

% the measure has three equally weighted mixture components
comp_num = 3;
meas_init_mixnorm = struct;
meas_init_mixnorm.weights = [1/3; 1/3; 1/3];
meas_init_mixnorm.components = cell(comp_num, 1);

meas_init_mixnorm.components{1} = struct;
meas_init_mixnorm.components{1}.mean_vec = [0; 0];
meas_init_mixnorm.components{1}.cov_mat = [9, 0; 0, 9];

meas_init_mixnorm.components{2} = struct;
meas_init_mixnorm.components{2}.mean_vec = [0; 0];
meas_init_mixnorm.components{2}.cov_mat = [16, 0; 0, 0.04];

meas_init_mixnorm.components{3} = struct;
meas_init_mixnorm.components{3}.mean_vec = [0; 0];
meas_init_mixnorm.components{3}.cov_mat = [0.04, 0; 0, 16];

% create a grid for the test functions
meas_init_testfuncs_vertices_cell = cell(test_num, 1);
meas_init_testfuncs_triangles_cell = cell(test_num, 1);

for test_id = 1:test_num
    testfuncs_x_num = marg_granularity_list(test_id);
    testfuncs_y_num = marg_granularity_list(test_id);

    meas_init_testfuncs_x = linspace(-1, 1, testfuncs_x_num)';
    meas_init_testfuncs_y = linspace(-1, 1, testfuncs_y_num)';
    [meas_init_testfuncs_grid_x, meas_init_testfuncs_grid_y] = ...
        meshgrid(meas_init_testfuncs_x, meas_init_testfuncs_y);
    meas_init_testfuncs_vertices = [meas_init_testfuncs_grid_x(:), ...
        meas_init_testfuncs_grid_y(:)];
    meas_init_testfuncs_vertices = meas_init_testfuncs_vertices ...
        * [1/2, 1/2; -1/2, 1/2];
    meas_init_testfuncs_vertices = unique([meas_init_vertices; ...
        meas_init_testfuncs_vertices], 'rows', 'stable');
    testfuncs_vertices_lower_list = ...
        find(meas_init_testfuncs_vertices(:, 2) <= 1e-10);
    testfuncs_vertices_upper_list = ...
        find(meas_init_testfuncs_vertices(:, 2) >= -1e-10);
    testfuncs_triangles_lower = delaunay(meas_init_testfuncs_vertices( ...
        testfuncs_vertices_lower_list, :));
    testfuncs_triangles_lower = cleanup_triangles( ...
        meas_init_testfuncs_vertices(testfuncs_vertices_lower_list, :), ...
        testfuncs_triangles_lower, [], false);
    testfuncs_triangles_upper = delaunay(meas_init_testfuncs_vertices( ...
        testfuncs_vertices_upper_list, :));
    testfuncs_triangles_upper = cleanup_triangles( ...
        meas_init_testfuncs_vertices(testfuncs_vertices_upper_list, :), ...
        testfuncs_triangles_upper, [], false);

    meas_init_testfuncs_vertices_cell{test_id} = ...
        meas_init_testfuncs_vertices;
    meas_init_testfuncs_triangles_cell{test_id} = ...
        {testfuncs_vertices_lower_list(testfuncs_triangles_lower); ...
        testfuncs_vertices_upper_list(testfuncs_triangles_upper)};
end


Meas_init = ProbMeas2D_MixNorm(meas_init_vertices, ...
    meas_init_triangles, meas_init_mixnorm);

meas_init_meanvec = Meas_init.meanVector();
meas_init_covmat = Meas_init.covarianceMatrix();
meas_init_trans = inv(chol(meas_init_covmat))';
meas_init_shift = -meas_init_trans * meas_init_meanvec;

% create the reference measure which is the affine transformation of the
% initial measure which has zero mean and unit covariance
meas_ref_vertices = meas_init_vertices * meas_init_trans' ...
    + meas_init_shift';
meas_ref_triangles = meas_init_triangles;

meas_ref_mixnorm = struct;
meas_ref_mixnorm.weights = meas_init_mixnorm.weights;
meas_ref_mixnorm.components = cell(comp_num, 1);

for comp_id = 1:comp_num
    meas_ref_mixnorm.components{comp_id} = struct;
    meas_ref_mixnorm.components{comp_id}.mean_vec = meas_init_trans ...
        * meas_init_mixnorm.components{comp_id}.mean_vec + meas_init_shift;
    meas_ref_mixnorm.components{comp_id}.cov_mat = meas_init_trans ...
        * meas_init_mixnorm.components{comp_id}.cov_mat * meas_init_trans';
end

% grid for the test functions
meas_ref_testfuncs_vertices_cell = cell(test_id, 1);
meas_ref_testfuncs_triangles_cell = meas_init_testfuncs_triangles_cell;

for test_id = 1:test_num
    meas_ref_testfuncs_vertices_cell{test_id} = ...
        meas_init_testfuncs_vertices_cell{test_id} * meas_init_trans' ...
        + meas_init_shift';
end


% generate the marginals via shifting and shearing
marg_num = 5;

marg_weights = ones(marg_num, 1) / marg_num;

% randomly generate the mean vectors and covariance matrices
locationscatter_cell = cell(marg_num, 3);

for marg_id = 1:marg_num
    locationscatter_cell{marg_id, 1} = randn(2, 1) * 0.5;
    eigvals = gamrnd(3, 0.5, 2, 1);
    eigvec_angle = rand(1, 1) * 2 * pi;
    eigvec = [cos(eigvec_angle), -sin(eigvec_angle); ...
        sin(eigvec_angle), cos(eigvec_angle)];
    locationscatter_cell{marg_id, 2} = eigvec * diag(eigvals) * eigvec';
    locationscatter_cell{marg_id, 3} = ...
        locationscatter_cell{marg_id, 2} ^ 2;
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
    trans_mat = locationscatter_cell{marg_id, 2};
    shift_vec = locationscatter_cell{marg_id, 1};
    marg_verts = meas_ref_vertices * trans_mat' + shift_vec';
    marg_triangles = meas_ref_triangles;

    marg_mixnorm = struct;
    marg_mixnorm.weights = meas_ref_mixnorm.weights;
    marg_mixnorm.components = cell(comp_num, 1);

    for comp_id = 1:comp_num
        marg_mixnorm.components{comp_id} = struct;
        marg_mixnorm.components{comp_id}.mean_vec = trans_mat ...
            * meas_ref_mixnorm.components{comp_id}.mean_vec + shift_vec;
        marg_mixnorm.components{comp_id}.cov_mat = trans_mat ...
            * meas_ref_mixnorm.components{comp_id}.cov_mat * trans_mat';
    end

    marg_vertices_cell{marg_id} = marg_verts;
    marg_triangles_cell{marg_id} = marg_triangles;
    marg_mixnorm_cell{marg_id} = marg_mixnorm;

    for test_id = 1:test_num
        % grid for the test functions
        marg_testfuncs_vertices = ...
            meas_ref_testfuncs_vertices_cell{test_id} * trans_mat' ...
            + shift_vec';
        marg_testfuncs_triangles = ...
            meas_ref_testfuncs_triangles_cell{test_id};
    
        marg_testfuncs_cell{test_id}{marg_id} = ...
            {marg_testfuncs_vertices, marg_testfuncs_triangles, false};
    end

    % instantiate the marginal
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

% grid used to plot the histogram and density
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
    'meas_ref_vertices', ...
    'meas_ref_triangles', ...
    'meas_ref_mixnorm', ...
    'meas_ref_testfuncs_vertices_cell', ...
    'meas_ref_testfuncs_triangles_cell', ...
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


