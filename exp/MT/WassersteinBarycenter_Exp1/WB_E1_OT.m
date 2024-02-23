% Compute semi-discrete optimal transport of each marginal

CONFIG = WB_E1_config();

load(CONFIG.SAVEPATH_INPUTS);

options = struct;
options.log_file = CONFIG.LOGPATH_OT;
options.OT = struct;
options.OT.optimization_options = struct;
options.OT.optimization_options.Display = 'iter-detailed';

% larger numbers of angles are used when the discrete measure has fewer
% atoms
angle_num_list = [1e5; 1e4; 5e3; 2e3; 1e3; 1e3];
optimality_tolerance_list = [1e-5; 1e-5; 1e-5; 1e-6; 1e-6; 1e-6];

marg_cell = cell(marg_num, 1);

for marg_id = 1:marg_num
    marg_cell{marg_id} = ProbMeas2D_MixNorm( ...
        marg_vertices_cell{marg_id}, ...
        marg_triangles_cell{marg_id}, ...
        marg_mixnorm_cell{marg_id});
end

OT_info_cell = cell(test_num, 1);

for test_id = 1:test_num
    options.OT.angle_num = angle_num_list(test_id);
    options.OT.optimization_options.OptimalityTolerance = ...
        optimality_tolerance_list(test_id);

    for marg_id = 1:marg_num
        marg_cell{marg_id}.setSimplicialTestFuncs( ...
            marg_testfuncs_cell{test_id}{marg_id}{:});
    end

    MT = MT2DQuad_ParTrans(marg_cell, marg_weights, ...
        {quality_vertices}, [], options);
    
    OT_info_cell{test_id} = MT.performReassembly();

    save(CONFIG.SAVEPATH_OT, ...
        'OT_info_cell', ...
        '-v7.3');
end