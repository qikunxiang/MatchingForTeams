% Compute semi-discrete optimal transport of each marginal

CONFIG = BizLoc_config();

load(CONFIG.SAVEPATH_INPUTS);
load(CONFIG.SAVEPATH_OUTPUTS);

options = struct;
options.log_file = CONFIG.LOGPATH_OT;
options.OT = struct;
options.OT.optimization_options = struct;
options.OT.optimization_options.Display = 'iter-detailed';

% larger numbers of angles are used when the discrete measure has fewer
% atoms
angle_num_list = [1e5; 5e4; 1e4; 1e4; 5e3];
optimality_tolerance_list = ...
    [1e-5; 1e-5; 3e-6; 1e-6; 1e-6];
marginal_magnetic_grids = [2; 4; 8; 16; 32] * 100;

marg_cell = cell(marg_num, 1);

for marg_id = 1:marg_num
    marg_cell{marg_id} = ProbMeas2D_CPWADens( ...
        marg_vertices_cell{marg_id}, ...
        marg_triangles_cell{marg_id}, ...
        marg_density_cell{marg_id});
end

OT_info_cell = cell(test_num, 1);

for test_id = 1:test_num
    options.marginal_magnetic_grids = ones(marg_num, 2) ...
        * marginal_magnetic_grids(test_id);

    options.OT.angle_num = angle_num_list(test_id);
    options.OT.optimization_options.OptimalityTolerance = ...
        optimality_tolerance_list(test_id);

    for marg_id = 1:marg_num
        marg_cell{marg_id}.setSimplicialTestFuncs( ...
            marg_testfuncs_cell{test_id}{marg_id}{:});
    end

    MT = MT2DBizLoc_ParTrans(marg_cell, costfuncs, ...
            quality_vertices, ...
            quality_testfuncs_cell{test_id}, ...
            options, [], []);
    MT.setLSIPSolutions(LSIP_primal_cell{test_id}, ...
        LSIP_dual_cell{test_id}, ...
        LSIP_UB_list(test_id), LSIP_LB_list(test_id));

    MT.prepareDiscreteOT();
    OT_info_cell{test_id} = MT.prepareSemiDiscreteOT();

    save(CONFIG.SAVEPATH_OT, ...
        'OT_info_cell', ...
        '-v7.3');
end