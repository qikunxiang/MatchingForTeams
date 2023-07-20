save_file_path = ['experiments/experiment2/quality10d/' ...
    'exp2_q10d_inputs_T%02d.mat'];

% randomly generate 10 problem settings
randtrial_num = 10;

% maximum number of marginals to test
marg_max_num = 18;

% the support of all marginals are [0, 1], and the density of every
% marginal is continuous piece-wise affine on [0, 0.25], [0.25, 0.5],
% [0.5, 0.75], [0.75, 1]
marg_bounds = [0, 1];
marg_knots = linspace(0, 1, 5)';

% the quality space is the simplex formed by [0; 0], [0; 1]; [1; 0]
quality = struct;
quality.dim = 10;
quality.aux_num = 0;
quality.ineq_A = sparse( ...
    ones(1, 10));
quality.ineq_rhs = 1;
quality.lb = zeros(10, 1);

% the number of test functions
testfunc_num = 50;

testfuncs_cell = cell(marg_max_num, 1);

for marg_id = 1:marg_max_num
    marg_grid_pts = linspace(marg_bounds(1), marg_bounds(2), ...
        testfunc_num)';
    testfuncs_cell{marg_id} = {marg_grid_pts};
end

for trial_id = 1:randtrial_num

    rng(20000 + trial_id * 100, 'combRecursive');
    
    marg_knots_cell = cell(marg_max_num, 1);
    marg_dens_cell = cell(marg_max_num, 1);
    
    for marg_id = 1:marg_max_num
        marg_knots_cell{marg_id} = marg_knots;
        marg_dens_cell{marg_id} = gamrnd(2, 1, length(marg_knots), 1);
    end
    
    costfunc_cell = cell(marg_max_num, 1);
    
    for marg_id = 1:marg_max_num
        weights = randn(10, 1);
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
    
    save(sprintf(save_file_path, trial_id), ...
        'marg_knots_cell', 'marg_dens_cell', ...
        'costfunc_cell', 'testfuncs_cell', 'quality');
end