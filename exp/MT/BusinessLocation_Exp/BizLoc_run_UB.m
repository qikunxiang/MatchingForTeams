% Compute semi-discrete optimal transport of each marginal

CONFIG = BizLoc_config();

load(CONFIG.SAVEPATH_INPUTS);
load(CONFIG.SAVEPATH_OUTPUTS);
load(CONFIG.SAVEPATH_OT);

options = struct;
options.log_file = CONFIG.LOGPATH_UB;

MCsamp_num = 1e4;
MCrep_num = 100;

tolerance = 2e-4;

main_log_file = fopen(CONFIG.LOGPATH_MAIN, 'a');

if main_log_file < 0
    error('cannot open log file');
end

fprintf(main_log_file, '--- experiment starts ---\n');

marg_cell = cell(marg_num, 1);

for marg_id = 1:marg_num
    marg_cell{marg_id} = ProbMeas2D_CPWADens( ...
        marg_vertices_cell{marg_id}, ...
        marg_triangles_cell{marg_id}, ...
        marg_density_cell{marg_id});
end

MT_LB_list = zeros(test_num, 1);
MT_UB_cell = cell(test_num, 2);
MT_UB_mean_list = zeros(test_num, 2);
MT_diff_cell = cell(test_num, 2);
MT_diff_mean_list = zeros(test_num, 2);
MT_OTEB_list = zeros(test_num, 1);
MT_THEB_list = zeros(test_num, 1);
MT_histpdf_cell = cell(test_num, 1);

for test_id = 1:test_num
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

    MT.loadOptimalTransportInfo(OT_info_cell{test_id});

    MT_LB = MT.getMTLowerBound();

    RS = RandStream('mrg32k3a', 'Seed', 1000);
    RS.Substream = test_id;
    [MT_UB_list_disc, MT_UB_list_cont, MCsamp, MT_histpdf] ...
        = MT.getMTUpperBoundsWRepetition( ...
        MCsamp_num, MCrep_num, RS, [], ...
        quality_hist_edge_x, quality_hist_edge_y);
    MT_UB_mean = zeros(2, 1);
    MT_UB_mean(1) = mean(MT_UB_list_disc);
    MT_UB_mean(2) = mean(MT_UB_list_cont);
    MT_diff_list1 = MT_UB_list_disc - MT_LB;
    MT_diff_list2 = MT_UB_list_cont - MT_LB;
    MT_diff_mean = zeros(2, 1);
    MT_diff_mean(1) = mean(MT_diff_list1);
    MT_diff_mean(2) = mean(MT_diff_list2);

    MT_OTEB = MT.getMTErrorBoundBasedOnOT();
    MT_THEB = MT.getMTTheoreticalErrorBound(tolerance);

    log_text = sprintf(['Test %d: LB = %10.4f, ' ...
        'UB1 = %10.4f, UB2 = %10.4f, ' ...
        'diff1 = %10.6f, diff2 = %10.6f, ' ...
        'OTEB = %10.6f, THEB = %10.6f\n'], ...
        test_id, MT_LB, ...
        MT_UB_mean(1), MT_UB_mean(2), ...
        MT_diff_mean(1), MT_diff_mean(2), ...
        MT_OTEB, MT_THEB);

    fprintf(main_log_file, log_text);
    fprintf(log_text);

    MT_LB_list(test_id) = MT_LB;
    MT_UB_cell{test_id, 1} = MT_UB_list_disc;
    MT_UB_cell{test_id, 2} = MT_UB_list_cont;
    MT_UB_mean_list(test_id, :) = MT_UB_mean';
    MT_diff_cell{test_id, 1} = MT_diff_list1;
    MT_diff_cell{test_id, 2} = MT_diff_list2;
    MT_diff_mean_list(test_id, :) = MT_diff_mean';
    MT_OTEB_list(test_id) = MT_OTEB;
    MT_THEB_list(test_id) = MT_THEB;
    MT_histpdf_cell{test_id} = MT_histpdf;

    save(CONFIG.SAVEPATH_OUTPUTS_UB, ...
        'MT_LB_list', ...
        'MT_UB_cell', ...
        'MT_UB_mean_list', ...
        'MT_diff_cell', ...
        'MT_diff_mean_list', ...
        'MT_OTEB_list', ...
        'MT_THEB_list', ...
        'MT_histpdf_cell', ...
        '-v7.3');
end

fprintf(main_log_file, '--- experiment ends ---\n\n');
fclose(main_log_file);