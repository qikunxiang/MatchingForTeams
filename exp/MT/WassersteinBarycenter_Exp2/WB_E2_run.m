% Compute the lower and upper bounds for the 2-Wasserstein barycenter
% problem

CONFIG = WB_E2_config();

load(CONFIG.SAVEPATH_INPUTS);
load(CONFIG.SAVEPATH_OT);

options = struct;
options.log_file = CONFIG.LOGPATH_LSIP_MAIN;
options.display = true;
options.sanitation_threshold = 1e-8;
options.reduce = struct;
options.reduce.thres = 0.5;
options.reduce.max_iter = 2000;
options.reduce.freq = 25;

global_options = struct;
global_options.pool_size = 5;
global_options.low_memory_mode = true;
global_options.display = true;
global_options.log_file = CONFIG.LOGPATH_LSIP_GLOBAL;

LP_options = struct;
LP_options.OutputFlag = 1;
LP_options.LogToConsole = 1;
LP_options.LogFile = CONFIG.LOGPATH_LSIP_LP;

tolerance = 2e-4;
MCsamp_num = 1e6;
MCrep_num = 100;

output_cell = cell(test_num, 1);
LSIP_primal_cell = cell(test_num, 1);
LSIP_dual_cell = cell(test_num, 1);
LSIP_LB_list = zeros(test_num, 1);
LSIP_UB_list = zeros(test_num, 1);
MT_LB_list = zeros(test_num, 1);
MT_UB_cell = cell(test_num, 2);
MT_UB_mean_list = zeros(test_num, 2);
MT_diff_cell = cell(test_num, 2);
MT_diff_mean_list = zeros(test_num, 2);
MT_OTEB_list = zeros(test_num, 1);
MT_THEB_list = zeros(test_num, 1);
WB_histpdf_cell = cell(test_num, 1);


main_log_file = fopen(CONFIG.LOGPATH_MAIN, 'a');

if main_log_file < 0
    error('cannot open log file');
end

fprintf(main_log_file, '--- experiment starts ---\n');

marg_cell = cell(marg_num, 1);

for marg_id = 1:marg_num
    marg_cell{marg_id} = ProbMeas2D_MixNorm( ...
        marg_vertices_cell{marg_id}, ...
        marg_triangles_cell{marg_id}, ...
        marg_mixnorm_cell{marg_id});
end


for test_id = 1:test_num

    if test_id == 1
        for marg_id = 1:marg_num
            marg_cell{marg_id}.setSimplicialTestFuncs( ...
                marg_testfuncs_cell{test_id}{marg_id}{:});
        end
    
        MT = MT2DQuad_ParTrans(marg_cell, marg_weights, ...
            {quality_vertices}, ...
            quality_testfuncs_cell{test_id}, ...
            options, LP_options, global_options);
    
        coup_cell = MT.generateDiscreteCoupling();
    else
        coup_cell = MT.updateSimplicialTestFuncs( ...
            quality_testfuncs_cell{test_id}, ...
            marg_testfuncs_cell{test_id});
    end

    initial_constr = cell(marg_num, 1);

    for marg_id = 1:marg_num
        initial_constr{marg_id} = struct( ...
            'vertex_indices', coup_cell{marg_id}.vertex_indices, ...
            'points', coup_cell{marg_id}.points, ...
            'testfunc_vals', coup_cell{marg_id}.testfunc_vals, ...
            'costfunc_vals', coup_cell{marg_id}.costfunc_vals, ...
            'min_val', -inf);
    end

    output = MT.run(initial_constr, tolerance);

    MT.loadOptimalTransportInfo(OT_info_cell{test_id});

    LSIP_primal = MT.Runtime.PrimalSolution;
    LSIP_dual = MT.Runtime.DualSolution;
    LSIP_LB = MT.Runtime.LSIP_LB;
    LSIP_UB = MT.Runtime.LSIP_UB;

    MT_LB = MT.getMTLowerBound();

    RS = RandStream('mrg32k3a', 'Seed', 2000);
    RS.Substream = test_id;
    [MT_UB_list_disc, MT_UB_list_cont, MCsamp, WB_histpdf] ...
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

    MT_OTEB = MT.getMTErrorBoundBasedOnOT(tolerance);
    MT_THEB = MT.getMTTheoreticalErrorBound(tolerance);

    log_text = sprintf(['test %d: ' ...
        'LB = %10.4f, UB1 = %10.4f, UB2 = %10.4f, ' ...
        'diff1 = %10.6f, diff2 = %10.6f, ' ...
        'OTEB = %10.6f, THEB = %10.6f\n'], ...
        test_id, MT_LB, ...
        MT_UB_mean(1), MT_UB_mean(2), ...
        MT_diff_mean(1), MT_diff_mean(2), ...
        MT_OTEB, MT_THEB);

    fprintf(main_log_file, log_text);
    fprintf(log_text);

    output_cell{test_id} = output;
    LSIP_primal_cell{test_id} = LSIP_primal;
    LSIP_dual_cell{test_id} = LSIP_dual;
    LSIP_LB_list(test_id) = LSIP_LB;
    LSIP_UB_list(test_id) = LSIP_UB;
    MT_LB_list(test_id) = MT_LB;
    MT_UB_cell{test_id, 1} = MT_UB_list_disc;
    MT_UB_cell{test_id, 2} = MT_UB_list_cont;
    MT_UB_mean_list(test_id, :) = MT_UB_mean';
    MT_diff_cell{test_id, 1} = MT_diff_list1;
    MT_diff_cell{test_id, 2} = MT_diff_list2;
    MT_diff_mean_list(test_id, :) = MT_diff_mean';
    MT_OTEB_list(test_id) = MT_OTEB;
    MT_THEB_list(test_id) = MT_THEB;
    WB_histpdf_cell{test_id} = WB_histpdf;

    save(CONFIG.SAVEPATH_OUTPUTS, ...
        'output_cell', ...
        'LSIP_primal_cell', ...
        'LSIP_dual_cell', ...
        'LSIP_LB_list', ...
        'LSIP_UB_list', ...
        'MT_LB_list', ...
        'MT_UB_cell', ...
        'MT_UB_mean_list', ...
        'MT_diff_cell', ...
        'MT_diff_mean_list', ...
        'MT_OTEB_list', ...
        'MT_THEB_list', ...
        'WB_histpdf_cell', ...
        '-v7.3');
end

fprintf(main_log_file, '--- experiment ends ---\n\n');
fclose(main_log_file);