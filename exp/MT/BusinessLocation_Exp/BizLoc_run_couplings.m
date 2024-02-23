% Plot the computed couplings

CONFIG = BizLoc_config();

load(CONFIG.SAVEPATH_INPUTS);
load(CONFIG.SAVEPATH_OUTPUTS);
load(CONFIG.SAVEPATH_OT);

MCsamp_num = 1e3;

marg_cell = cell(marg_num, 1);

for marg_id = 1:marg_num
    marg_cell{marg_id} = ProbMeas2D_CPWADens( ...
        marg_vertices_cell{marg_id}, ...
        marg_triangles_cell{marg_id}, ...
        marg_density_cell{marg_id});
end

RS = RandStream('mrg32k3a', 'Seed', 9999);

test_id = 5;

for marg_id = 1:marg_num
    marg_cell{marg_id}.setSimplicialTestFuncs( ...
        marg_testfuncs_cell{test_id}{marg_id}{:});
end

MT = MT2DBizLoc_ParTrans(marg_cell, costfuncs, ...
    quality_vertices, ...
    quality_testfuncs_cell{test_id}, ...
    [], [], []);
MT.setLSIPSolutions(LSIP_primal_cell{test_id}, ...
    LSIP_dual_cell{test_id}, ...
    LSIP_UB_list(test_id), LSIP_LB_list(test_id));

MT.loadOptimalTransportInfo(OT_info_cell{test_id});

[~, ~, MCsamp] = MT.getMTUpperBounds(MCsamp_num, RS);

save(CONFIG.SAVEPATH_COUPLINGS, ...
    'MCsamp', ...
    '-v7.3');