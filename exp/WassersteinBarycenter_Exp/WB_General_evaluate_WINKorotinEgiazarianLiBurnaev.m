% Evaluate the performance of the Wasserstein Iterative Network method of Korotin, Egiazarian, Li, Burnaev via semi-discrete 
% Wasserstein-2 OT

CONFIG = WB_General_config();

load(CONFIG.SAVEPATH_INPUTS);

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

test_epoch = 300;


MT = MatchTeam2DWassersteinBarycenter( ...
    marg_cell, ...
    marg_weights, ...
    {quality_vertices}, ...
    [], []);

barycenter_samples = readmatrix(sprintf(CONFIG.PATHFORMAT_RESULTS_WINKOROTINEGIAZARIANLIBURNAEV, test_epoch));
samp_num = size(barycenter_samples, 1);

[objective, W2_dist_list] = MT.evaluateObjectiveViaW2( ...
    barycenter_samples, ones(samp_num, 1) / samp_num, optimization_options);

fprintf('epoch %d: objective = %10.4f\n', test_epoch, objective);

save(CONFIG.SAVEPATH_EVALUATION_WINKOROTINEGIAZARIANLIBURNAEV, ...
    'objective', ...
    'W2_dist_list', ...
    '-v7.3');