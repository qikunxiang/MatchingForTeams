% Example benchmarking the running time of the proposed algorithm with 1D
% marginals and continuous piece-wise affine cost functions

CONFIG = BM1D_config();

% randomly generate 10 problem settings
randtrial_num = 10;

% maximum number of marginals to test
marg_max_num = 100;

% the support of all marginals are [0, 1], and the density of every
% marginal is continuous piece-wise affine on [0, 0.25], [0.25, 0.5],
% [0.5, 0.75], [0.75, 1]
marg_bounds = [0, 1];
marg_knots = linspace(0, 1, 5)';

% the quality space is the simplex formed by [0; 0], [0; 1]; [1; 0]
quality = struct;
quality.dim = 2;
quality.aux_num = 0;
quality.ineq_A = sparse( ...
    [1, 1]);
quality.ineq_rhs = 1;
quality.lb = [0; 0];

quality_vertices = [0, 0; 0, 1; 1, 0];

% the number of test functions
testfunc_num = 50;

testfuncs_cell = cell(marg_max_num, 1);

for marg_id = 1:marg_max_num
    marg_grid_pts = linspace(marg_bounds(1), marg_bounds(2), ...
        testfunc_num)';
    testfuncs_cell{marg_id} = {marg_grid_pts};
end

q_grid_pt_x = (0:1/32:1)';
q_grid_pt_y = (0:1/32:1)';
[q_grid_x, q_grid_y] = meshgrid(q_grid_pt_x, q_grid_pt_y);
q_grid_pts = [q_grid_x(:), q_grid_y(:)];
q_testfunc_vertices = q_grid_pts(sum(q_grid_pts, 2) <= 1, :);
q_testfunc_triangles = delaunay(q_testfunc_vertices);
q_tri_num = size(q_testfunc_triangles, 1);
q_degen_list = true(q_tri_num, 1);

for tri_id = 1:q_tri_num
    cond_num = cond([q_testfunc_vertices( ...
        q_testfunc_triangles(tri_id, :), :)'; ones(1, 3)]);
    if cond_num <= 1e12
        q_degen_list(tri_id) = false;
    end
end

q_testfunc_triangles = q_testfunc_triangles(~q_degen_list, :);

quality_testfuncs = {q_testfunc_vertices, q_testfunc_triangles, false};

for trial_id = 1:randtrial_num

    rng(1000 + trial_id * 100, 'combRecursive');
    
    marg_knots_cell = cell(marg_max_num, 1);
    marg_dens_cell = cell(marg_max_num, 1);
    
    for marg_id = 1:marg_max_num
        marg_knots_cell{marg_id} = marg_knots;
        marg_dens_cell{marg_id} = gamrnd(2, 1, length(marg_knots), 1);
    end
    
    costfunc_cell = cell(marg_max_num, 1);
    
    for marg_id = 1:marg_max_num
        weights = randn(2, 1);
        weights = weights / norm(weights);
    
        thres_list = rand(2, 1);
        thres1 = min(thres_list);
        thres2 = max(thres_list);

        cf_knots = [-2; -thres2; -thres1; thres1; thres2; 2];
        cf_vals = [thres2 - thres1; thres2 - thres1; ...
            0; 0; thres2 - thres1; thres2 - thres1];
        costfunc_cell{marg_id} = struct('weights', weights, ...
            'knots', cf_knots, 'values', cf_vals);
    end
    
    save(sprintf(CONFIG.SAVEPATH_INPUTS, trial_id), ...
        'marg_knots_cell', 'marg_dens_cell', ...
        'costfunc_cell', 'testfuncs_cell', ...
        'quality', 'quality_vertices', 'quality_testfuncs');
end