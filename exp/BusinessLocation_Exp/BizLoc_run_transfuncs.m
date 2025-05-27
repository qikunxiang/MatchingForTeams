% Compute semi-discrete optimal transport of each marginal

CONFIG = BizLoc_config();

load(CONFIG.SAVEPATH_INPUTS);
load(CONFIG.SAVEPATH_OUTPUTS);

ref_pt = [0; 0];

options = struct;
options.log_file = CONFIG.LOGPATH_TRANSFUNCS;
options.display = true;

main_log_file = fopen(CONFIG.LOGPATH_MAIN, 'a');

if main_log_file < 0
    error('cannot open log file');
end

fprintf(main_log_file, '--- experiment starts ---\n');

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
            options, [], []);
    MT.setLSIPSolutions(LSIP_primal_cell{test_id}, ...
        LSIP_dual_cell{test_id}, ...
        LSIP_UB_list(test_id), LSIP_LB_list(test_id));

    MT_transfuncs_cell{test_id} = MT.evaluateOptTransferFuncs(quality_transfuncs_grid, ref_pt);

    log_text = sprintf('Test %d: transfer functions computed\n', test_id);

    fprintf(main_log_file, log_text);
    fprintf(log_text);

    save(CONFIG.SAVEPATH_OUTPUTS_TRANSFUNCS, ...
        'MT_transfuncs_cell', ...
        '-v7.3');
end

fprintf(main_log_file, '--- experiment ends ---\n\n');
fclose(main_log_file);