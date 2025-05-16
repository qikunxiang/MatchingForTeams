% Evaluate the performance of our algorithm via semi-discrete Wasserstein-2 OT

CONFIG = WB_General_config();

load(CONFIG.SAVEPATH_INPUTS);
load(CONFIG.SAVEPATH_OUTPUTS_LB);

options.W2OT = struct;
options.W2OT.optimization_options = struct;
options.W2OT.optimization_options.OptimalityTolerance = 1e-5;
options.W2OT.optimization_options.Display = 'iter-detailed';

test_id = 6;

marg_cell = cell(marg_num, 1);

for marg_id = 1:marg_num
    marg_cell{marg_id} = ProbMeas2DCPWADens( ...
        marg_vertices_cell{marg_id}, ...
        marg_triangles_cell{marg_id}, ...
        marg_densities_cell{marg_id});
    marg_cell{marg_id}.setSimplicialTestFuncs(marg_testfuncs_cell{test_id}{marg_id}{:});
    
end


MT = MatchTeam2DWassersteinBarycenter( ...
    marg_cell, ...
    marg_weights, ...
    {quality_vertices}, ...
    [], options);

MT.setLSIPSolutions( ...
    LSIP_primal_cell{test_id}, ...
    LSIP_dual_cell{test_id}, ...
    LSIP_UB_list(test_id), ...
    LSIP_LB_list(test_id));

[~, W2OT_info] = MT.computeW2OptimalCouplings();

W2_dist_list = zeros(marg_num, 1);

for marg_id = 1:marg_num
    W2_dist_list(marg_id) = W2OT_info{marg_id}.cost;
end

objective = MT.getMTUpperBoundViaW2DiscreteCoupling();

fprintf('objective = %10.4f\n', objective);

save(CONFIG.SAVEPATH_EVALUATION_OURALGO, ...
    'objective', ...
    'W2_dist_list', ...
    '-v7.3');