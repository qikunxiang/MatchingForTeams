% Compute semi-discrete optimal transport of each marginal

CONFIG = BizLoc_config();

load(CONFIG.SAVEPATH_INPUTS);
load(CONFIG.SAVEPATH_OUTPUTS);


marg_cell = cell(marg_num, 1);

for marg_id = 1:marg_num
    marg_cell{marg_id} = ProbMeas2DCPWADens( ...
        marg_vertices_cell{marg_id}, ...
        marg_triangles_cell{marg_id}, ...
        marg_density_cell{marg_id});
end

MT_transfuncs_cell = cell(test_num, 1);


for test_id = 1:test_num
    for marg_id = 1:marg_num
        marg_cell{marg_id}.setSimplicialTestFuncs(marg_testfuncs_cell{test_id}{marg_id}{:});
    end

    MT = MatchTeam2DBusinessLocation( ...
            marg_cell, costfuncs, ...
            quality_vertices, ...
            quality_testfuncs_cell{test_id}, ...
            [], [], []);
    MT.setLSIPSolutions(LSIP_primal_cell{test_id}, ...
        LSIP_dual_cell{test_id}, ...
        LSIP_UB_list(test_id), LSIP_LB_list(test_id));


    MIP_largest_C = 0;
    MIP_largest_B = 0;
    for marg_id = 1:marg_num
        gurobimodel = MT.Storage.GlobalMin{marg_id}.GurobiModel;

        MIP_C = sum(gurobimodel.vtype == 'C');
        MIP_B = sum(gurobimodel.vtype == 'B');

        if MIP_C > MIP_largest_C
            MIP_largest_C = MIP_C;
            MIP_largest_B = MIP_B;
        end
    end

    fprintf('Test %d: largest MIP model contains %d continuous variables and %d binary variables.\n', ...
        test_id, MIP_largest_C, MIP_largest_B);

    MT_transfuncs_cell{test_id} = MT.evaluateOptTransferFuncs(quality_transfuncs_grid);

    log_text = sprintf('Test %d: transfer functions computed\n', test_id);

    fprintf(log_text);

    save(CONFIG.SAVEPATH_OUTPUTS_TRANSFUNCS, ...
        'MT_transfuncs_cell', ...
        '-v7.3');
end