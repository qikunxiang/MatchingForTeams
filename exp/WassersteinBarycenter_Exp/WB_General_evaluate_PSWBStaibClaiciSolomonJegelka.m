% Evaluate the performance of the Parallel Streaming Wasserstein Barycenter (PSWB) method of Staib, Claici, Solomon, Jegelka via 
% semi-discrete Wasserstein-2 OT

CONFIG = WB_General_config();

load(CONFIG.SAVEPATH_INPUTS);

grid = quality_testfuncs_cell{end}{1};

optimization_options = struct;
optimization_options.OptimalityTolerance = 1e-5;
optimization_options.Display = 'iter-detailed';

marg_cell = cell(marg_num, 1);

for marg_id = 1:marg_num
    marg_cell{marg_id} = ProbMeas2DCPWADens( ...
        marg_vertices_cell{marg_id}, ...
        marg_triangles_cell{marg_id}, ...
        marg_densities_cell{marg_id});
end

test_epoch = 10000000;


MT = MatchTeam2DWassersteinBarycenter( ...
    marg_cell, ...
    marg_weights, ...
    {quality_vertices}, ...
    [], []);

barycenter_weights = readmatrix(sprintf(CONFIG.PATHFORMAT_RESULTS_PSWBSTAIBCLAICISOLOMONJEGELKA, test_epoch));

grid_positive = barycenter_weights >= 1e-6;
subgrid = grid(grid_positive, :);
barycenter_weights = barycenter_weights(grid_positive);
barycenter_weights = barycenter_weights / sum(barycenter_weights);

[objective, W2_dist_list] = MT.evaluateObjectiveViaW2( ...
    subgrid, barycenter_weights, optimization_options);

fprintf('epoch %d: objective = %10.4f\n', test_epoch, objective);

save(CONFIG.SAVEPATH_EVALUATION_PSWBSTAIBCLAICISOLOMONJEGELKA, ...
    'objective', ...
    'W2_dist_list', ...
    '-v7.3');