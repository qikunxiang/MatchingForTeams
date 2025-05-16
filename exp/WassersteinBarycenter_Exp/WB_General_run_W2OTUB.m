% Compute the lower and upper bounds for the 2-Wasserstein barycenter problem

CONFIG = WB_General_config();

load(CONFIG.SAVEPATH_INPUTS);
load(CONFIG.SAVEPATH_OUTPUTS_LB);
load(CONFIG.SAVEPATH_OUTPUTS_W2OT);

options = struct;
options.log_file = CONFIG.LOGPATH_W2OTUB;
options.display = true;

enable_quadratic_testfuncs = false;
tolerance = 2e-4;

MCsamp_num = 1e6;
MCrep_num = 100;

MT_LB_list = zeros(test_num, 1);
MT_W2OTUB_disc_list = zeros(test_num, 1);
MT_W2OTdiff_disc_list = zeros(test_num, 1);
MT_W2OTUB_cell = cell(test_num, 1);
MT_W2OTUB_mean_list = zeros(test_num, 1);
MT_W2OTdiff_cell = cell(test_num, 1);
MT_W2OTdiff_mean_list = zeros(test_num, 1);
MT_THEB_list = zeros(test_num, 1);
WB_W2OT_histpdf_cell = cell(test_num, 1);


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
        marg_densities_cell{marg_id});
end


for test_id = 1:test_num
    quality_testfuncs_cell{test_id}{3} = enable_quadratic_testfuncs; %#ok<SAGROW>
    
    for marg_id = 1:marg_num
        marg_cell{marg_id}.setSimplicialTestFuncs(marg_testfuncs_cell{test_id}{marg_id}{:});
    end

    MT = MatchTeam2DWassersteinBarycenter( ...
        marg_cell, ...
        marg_weights, ...
        {quality_vertices}, ...
        quality_testfuncs_cell{test_id}, ...
        options, [], []);
    MT.setLSIPSolutions( ...
        LSIP_primal_cell{test_id}, ...
        LSIP_dual_cell{test_id}, ...
        LSIP_UB_list(test_id), ...
        LSIP_LB_list(test_id));
    MT.loadW2OptimalTransportInfo(W2OT_info_cell{test_id, 1}, W2OT_info_cell{test_id, 2});
    
    MT_LB = MT.getMTLowerBound();
    MT_W2OTUB_disc = MT.getMTUpperBoundViaW2DiscreteCoupling();
    MT_W2OTdiff_disc = MT_W2OTUB_disc - MT_LB;

    RS = RandStream('mrg32k3a', 'Seed', 3000);
    RS.Substream = test_id;
    [MT_W2OTUB_list, MCsamp, WB_histpdf] ...
        = MT.getMTUpperBoundViaW2CouplingWRepetition( ...
        MCsamp_num, ...
        MCrep_num, ...
        RS, ...
        quality_hist_edge_x, quality_hist_edge_y);
    MT_W2OTUB_mean = mean(MT_W2OTUB_list);
    MT_W2OTdiff_list = MT_W2OTUB_list - MT_LB;
    MT_W2OTdiff_mean = mean(MT_W2OTdiff_list);
    MT_THEB = MT.getMTTheoreticalErrorBound(tolerance);

    log_text = sprintf('test %2d: LB = %10.4f, UB1 = %10.4f, UB2 = %10.4f, diff1 = %10.6f, diff2 = %10.6f, THEB = %10.6f\n', ...
        test_id, ...
        MT_LB, ...
        MT_W2OTUB_disc, ...
        MT_W2OTUB_mean, ...
        MT_W2OTdiff_disc, ...
        MT_W2OTdiff_mean, ...
        MT_THEB);

    fprintf(main_log_file, log_text);
    fprintf(log_text);

    MT_LB_list(test_id) = MT_LB;
    MT_W2OTUB_disc_list(test_id) = MT_W2OTUB_disc;
    MT_W2OTdiff_disc_list(test_id) = MT_W2OTdiff_disc;
    MT_W2OTUB_cell{test_id} = MT_W2OTUB_list;
    MT_W2OTUB_mean_list(test_id) = MT_W2OTUB_mean;
    MT_W2OTdiff_cell{test_id} = MT_W2OTdiff_list;
    MT_W2OTdiff_mean_list(test_id) = MT_W2OTdiff_mean;
    MT_THEB_list(test_id) = MT_THEB;
    WB_W2OT_histpdf_cell{test_id} = WB_histpdf;

    save(CONFIG.SAVEPATH_OUTPUTS_W2OTUB, ...
        'MT_LB_list', ...
        'MT_W2OTUB_disc_list', ...
        'MT_W2OTdiff_disc_list', ...
        'MT_W2OTUB_cell', ...
        'MT_W2OTUB_mean_list', ...
        'MT_W2OTdiff_cell', ...
        'MT_W2OTdiff_mean_list', ...
        'MT_THEB_list', ...
        'WB_W2OT_histpdf_cell', ...
        '-v7.3');
end

fprintf(main_log_file, '--- experiment ends ---\n\n');
fclose(main_log_file);