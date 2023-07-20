input_file_path = 'experiments/experiment1/exp1_inputs.mat';
OT_file_path = 'experiments/experiment1/exp1_OT.mat';

save_file_path = 'experiments/experiment1/exp1_MMOT_rst.mat';

rng(1000, 'combRecursive');

load(input_file_path);
load(OT_file_path);

options = struct;
options.log_file = 'exp1_MMOT_algo.log';
options.display = true;
options.reduce = struct;
options.reduce.thres = 5;
options.reduce.max_iter = 2000;
options.reduce.freq = 100;
options.OT = struct;
options.OT.optimization_options = struct;
options.OT.optimization_options.Display = 'iter-detailed';

global_options = struct;
global_options.tolerance = 2e-3;
global_options.branching = [100; 100];
global_options.termination_thres = -0.1;
global_options.display = true;
global_options.log_file = 'exp1_MMOT_globalmin_bnb.log';

LP_options = struct;
LP_options.OutputFlag = 1;
LP_options.LogToConsole = 0;
LP_options.LogFile = 'exp1_MMOT_gurobi_LP.log';

tolerance = 5e-3;
MCsamp_num = 1e6;
MCrep_num = 100;

exp_log_file_path = 'exp1_MMOT_outputs.log';

test_num = length(testfuncs_cell);
output_cell = cell(test_num, 1);
MCsamp_cell = cell(test_num, 1);
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

parfunc_cell = cell(test_num, marg_num);
transfunc_cell = cell(test_num, 1);
quality_cont_hist_cell = cell(test_num, 1);

exp_log_file = fopen(exp_log_file_path, 'a');

if exp_log_file < 0
    error('cannot open log file');
end

fprintf(exp_log_file, '--- experiment starts ---\n');

for test_id = 1:test_num
    if test_id == 1
        marg_cell = cell(marg_num, 1);

        for marg_id = 1:marg_num
            marg_cell{marg_id} = ProbMeas2D_CPWADens( ...
                marg_vert_cell{marg_id}, ...
                marg_tri_cell{marg_id}, ...
                marg_dens_cell{marg_id});
        end

        MT = MT2DQuad_MMOT(marg_cell, marg_weights, quality_poly_cell, ...
            options, LP_options, global_options);
        MT.setSimplicialTestFuncs(testfuncs_cell{test_id});
        [coup_indices, ~, opt_quality] = MT.generateHeuristicCoupling();
    else
        [coup_indices, ~, opt_quality] ...
            = MT.updateSimplicialTestFuncs(testfuncs_cell{test_id});
    end

    initial_constr = struct('vertex_indices', coup_indices, ...
        'points', opt_quality);

    output_cell{test_id} = MT.run(initial_constr, tolerance);

    MT.setSemiDiscreteOTWeights(OT_cell{test_id});

    LSIP_primal_cell{test_id} = MT.Runtime.PrimalSolution;
    LSIP_dual_cell{test_id} = MT.Runtime.DualSolution;
    LSIP_LB_list(test_id) = MT.Runtime.LSIP_LB;
    LSIP_UB_list(test_id) = MT.Runtime.LSIP_UB;

    MT_LB_list(test_id) = MT.getMTLowerBound();
    [MT_UB_cell{test_id, 1}, MT_UB_cell{test_id, 2}, ...
        MCsamp_cell{test_id}] = MT.getMTUpperBoundsWRepetition( ...
        MCsamp_num, MCrep_num);
    MT_UB_mean_list(test_id, 1) = mean(MT_UB_cell{test_id, 1});
    MT_UB_mean_list(test_id, 2) = mean(MT_UB_cell{test_id, 2});
    MT_diff_cell{test_id, 1} = MT_UB_cell{test_id, 1} ...
        - MT_LB_list(test_id);
    MT_diff_cell{test_id, 2} = MT_UB_cell{test_id, 2} ...
        - MT_LB_list(test_id);
    MT_diff_mean_list(test_id, 1) = mean(MT_diff_cell{test_id, 1});
    MT_diff_mean_list(test_id, 2) = mean(MT_diff_cell{test_id, 2});

    MT_OTEB_list(test_id) = MT.getMTErrorBoundBasedOnOT(tolerance);
    MT_THEB_list(test_id) = MT.getMTTheoreticalErrorBound(tolerance);

    quality_cont_hist_cell{test_id} = histcounts2( ...
        MCsamp_cell{test_id}.ContinuousQualities(:, 1), ...
        MCsamp_cell{test_id}.ContinuousQualities(:, 2), ...
        quality_hist_edge_x, quality_hist_edge_y)';

    for marg_id = 1:marg_num
        parfunc_cell{test_id, marg_id} = ...
            nan(size(plot_pts_cell{marg_id}, 1), 1);
        parfunc_cell{test_id, marg_id}(plot_inside_cell{marg_id}, :) ...
            = MT.evaluateOptParametricFunc(marg_id, ...
            plot_pts_cell{marg_id}(plot_inside_cell{marg_id}, :));
    end

    transfunc_cell{test_id} = nan(size(plot_quality_pts, 1), marg_num);
    transfunc_cell{test_id}(plot_quality_inside, :) ...
        = MT.evaluateOptTransferFuncs( ...
        plot_quality_pts(plot_quality_inside, :));
    
    fprintf(exp_log_file, ...
        ['test %d: LB = %10.4f, UB1 = %10.4f, UB2 = %10.4f, ' ...
        'diff1 = %10.6f, diff2 = %10.6f, ' ...
        'OTEB = %10.6f, THEB = %10.6f\n'], ...
        test_id, MT_LB_list(test_id), ...
        MT_UB_mean_list(test_id, 1), MT_UB_mean_list(test_id, 2), ...
        MT_diff_mean_list(test_id, 1), MT_diff_mean_list(test_id, 2), ...
        MT_OTEB_list(test_id), MT_THEB_list(test_id));

    fprintf(['test %d: LB = %10.4f, UB1 = %10.4f, UB2 = %10.4f, ' ...
        'diff1 = %10.6f, diff2 = %10.6f, ' ...
        'OTEB = %10.6f, THEB = %10.6f\n'], ...
        test_id, MT_LB_list(test_id), ...
        MT_UB_mean_list(test_id, 1), MT_UB_mean_list(test_id, 2), ...
        MT_diff_mean_list(test_id, 1), MT_diff_mean_list(test_id, 2), ...
        MT_OTEB_list(test_id), MT_THEB_list(test_id));

    save(save_file_path, 'output_cell', ...
        'MCsamp_cell', 'LSIP_primal_cell', 'LSIP_dual_cell', ...
        'LSIP_UB_list', 'LSIP_LB_list', ...
        'MT_UB_cell', 'MT_UB_mean_list', 'MT_LB_list', ...
        'MT_diff_cell', 'MT_diff_mean_list', ...
        'MT_OTEB_list', 'MT_THEB_list', ...
        'quality_cont_hist_cell', 'parfunc_cell', 'transfunc_cell', ...
        '-v7.3');

end

fprintf(exp_log_file, '--- experiment ends ---\n\n');
fclose(exp_log_file);