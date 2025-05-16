% Compute Wasserstein-2 semi-discrete optimal transport of each marginal

CONFIG = WB_General_config();

load(CONFIG.SAVEPATH_INPUTS);
load(CONFIG.SAVEPATH_OUTPUTS_LB);

options = struct;
options.log_file = CONFIG.LOGPATH_W2OT;
options.W2OT = struct;
options.W2OT.optimization_options = struct;
options.W2OT.optimization_options.OptimalityTolerance = 1e-7;
options.W2OT.optimization_options.Display = 'iter-detailed';

marg_cell = cell(marg_num, 1);

for marg_id = 1:marg_num
    marg_cell{marg_id} = ProbMeas2DCPWADens( ...
        marg_vertices_cell{marg_id}, ...
        marg_triangles_cell{marg_id}, ...
        marg_densities_cell{marg_id});
end

W2OT_info_cell = cell(test_num, 2);

for test_id = 1:test_num
    for marg_id = 1:marg_num
        marg_cell{marg_id}.setSimplicialTestFuncs(marg_testfuncs_cell{test_id}{marg_id}{:});
    end

    MT = MatchTeam2DWassersteinBarycenter( ...
        marg_cell, ...
        marg_weights, ...
        {quality_vertices}, ...
        [], ...
        options);

    MT.setLSIPSolutions( ...
        LSIP_primal_cell{test_id}, ...
        LSIP_dual_cell{test_id}, ...
        LSIP_UB_list(test_id), ...
        LSIP_LB_list(test_id));

    [W2OT_info_cell{test_id, 1}, W2OT_info_cell{test_id, 2}] = MT.computeW2OptimalCouplings();

    save(CONFIG.SAVEPATH_OUTPUTS_W2OT, ...
        'W2OT_info_cell', ...
        '-v7.3');
end