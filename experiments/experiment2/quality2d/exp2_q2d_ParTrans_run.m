input_file_path = ['experiments/experiment2/quality2d/' ...
    'exp2_q2d_inputs_T%02d.mat'];
result_file_path = 'experiments/experiment2/quality2d/%s_T%02d_M%03d.mat';
file_prefix = 'ParTrans_results/RST_LOG_DLOG';

options = struct;
options.log_file = 'exp2_q2d_ParTrans_algo.log';
options.global_formulation = 'LOG_DLOG';
options.sanitation_threshold = 1e-8;
options.time_limit = 86400;
options.display = true;
options.reduce = struct;
options.reduce.thres = 0.01;
options.reduce.max_iter = 3000;
options.reduce.freq = 50;

global_options = struct;
global_options.TimeLimit = 86400;
global_options.BestObjStop = -0.2;
global_options.PoolSearchMode = 1;
global_options.PoolSolutions = 30;
global_options.OutputFlag = 1;
global_options.LogFile = 'exp2_q2d_gurobi_MIP.log';

LP_options = struct;
LP_options.OutputFlag = 1;
LP_options.LogToConsole = 0;
LP_options.LogFile = 'exp2_q2d_gurobi_LP.log';

tolerance = 5e-5;
MCsamp_num = 1e6;
MCrep_num = 10;

exp_log_file_path = 'exp2_q2d_ParTrans_outputs.log';

test_settings = ...
    [ ...
      1,   4; ...
      2,   4; ...
      3,   4; ...
      4,   4; ...
      5,   4; ...
      6,   4; ...
      7,   4; ...
      8,   4; ...
      9,   4; ...
     10,   4; ...
      1,   6; ...
      2,   6; ...
      3,   6; ...
      4,   6; ...
      5,   6; ...
      6,   6; ...
      7,   6; ...
      8,   6; ...
      9,   6; ...
     10,   6; ...
      1,   8; ...
      2,   8; ...
      3,   8; ...
      4,   8; ...
      5,   8; ...
      6,   8; ...
      7,   8; ...
      8,   8; ...
      9,   8; ...
     10,   8; ...
      1,  10; ...
      2,  10; ...
      3,  10; ...
      4,  10; ...
      5,  10; ...
      6,  10; ...
      7,  10; ...
      8,  10; ...
      9,  10; ...
     10,  10; ...
      1,  12; ...
      2,  12; ...
      3,  12; ...
      4,  12; ...
      5,  12; ...
      6,  12; ...
      7,  12; ...
      8,  12; ...
      9,  12; ...
     10,  12; ...
      1,  14; ...
      2,  14; ...
      3,  14; ...
      4,  14; ...
      5,  14; ...
      6,  14; ...
      7,  14; ...
      8,  14; ...
      9,  14; ...
     10,  14; ...
      1,  16; ...
      2,  16; ...
      3,  16; ...
      4,  16; ...
      5,  16; ...
      6,  16; ...
      7,  16; ...
      8,  16; ...
      9,  16; ...
     10,  16; ...
      1,  18; ...
      2,  18; ...
      3,  18; ...
      4,  18; ...
      5,  18; ...
      6,  18; ...
      7,  18; ...
      8,  18; ...
      9,  18; ...
     10,  18; ...
      1,  20; ...
      2,  20; ...
      3,  20; ...
      4,  20; ...
      5,  20; ...
      6,  20; ...
      7,  20; ...
      8,  20; ...
      9,  20; ...
     10,  20; ...
      1,  50; ...
      2,  50; ...
      3,  50; ...
      4,  50; ...
      5,  50; ...
      6,  50; ...
      7,  50; ...
      8,  50; ...
      9,  50; ...
     10,  50; ...
      1,  80; ...
      2,  80; ...
      3,  80; ...
      4,  80; ...
      5,  80; ...
      6,  80; ...
      7,  80; ...
      8,  80; ...
      9,  80; ...
     10,  80; ...
      1, 100; ...
      2, 100; ...
      3, 100; ...
      4, 100; ...
      5, 100; ...
      6, 100; ...
      7, 100; ...
      8, 100; ...
      9, 100; ...
     10, 100; ...
     ];

for st_id = 1:size(test_settings, 1)

    trial_id = test_settings(st_id, 1);
    marg_num = test_settings(st_id, 2);

    rng(50000 + trial_id * 1000 + marg_num, 'combRecursive');

    load(sprintf(input_file_path, trial_id));

    exp_log_file = fopen(exp_log_file_path, 'a');

    if exp_log_file < 0
        error('cannot open log file');
    end

    fprintf(exp_log_file, '--- experiment starts ---\n');

    marg_cell = cell(marg_num, 1);
    scaled_costfunc_cell = cell(marg_num, 1);

    for marg_id = 1:marg_num
        marg_cell{marg_id} = ProbMeas1D_CPWADens( ...
            marg_knots_cell{marg_id}, ...
            marg_dens_cell{marg_id});

        % scale the cost functions such that the theoretical error
        % bounds stay the same
        scaled_costfunc_cell{marg_id} = costfunc_cell{marg_id};
        scaled_costfunc_cell{marg_id}.values = ...
            scaled_costfunc_cell{marg_id}.values / marg_num;
    end

    MT = MT1DCPWA_ParTrans(marg_cell, scaled_costfunc_cell, ...
        quality_vertices, [], options, LP_options, global_options);
    MT.setSimplicialTestFuncs(quality_testfuncs, ...
        testfuncs_cell(1:marg_num));
    coup_cell = MT.generateDiscreteCoupling();

    initial_constr = cell(marg_num, 1);

    for marg_id = 1:marg_num
        initial_constr{marg_id} = struct( ...
            'inputs', coup_cell{marg_id}.inputs, ...
            'qualities', coup_cell{marg_id}.qualities, ...
            'marg_testfunc_vals', ...
            coup_cell{marg_id}.marg_testfunc_vals, ...
            'quality_testfunc_vals', ...
            coup_cell{marg_id}.quality_testfunc_vals, ...
            'costfunc_vals', coup_cell{marg_id}.costfunc_vals, ...
            'min_val', -inf);
    end

    output = MT.run(initial_constr, tolerance);
    total_time = output.total_time;

    if output.time_limit_exceeded
        fprintf(exp_log_file, ['trial %2d, marg num = %3d: ' ...
            'time limit exceeded\n'], trial_id, marg_num);

        fprintf(['trial %2d, marg num = %3d: ' ...
            'time limit exceeded\n'], trial_id, marg_num);
        fprintf(exp_log_file, '--- experiment ends ---\n\n');
        fclose(exp_log_file);
        continue;
    end

    LSIP_primal = MT.Runtime.PrimalSolution;
    LSIP_dual = MT.Runtime.DualSolution;
    LSIP_LB = MT.Runtime.LSIP_LB;
    LSIP_UB = MT.Runtime.LSIP_UB;

    MT_LB = MT.getMTLowerBound();

    [MT_UB_list, MCsamps] ...
        = MT.getMTUpperBoundDiscreteWRepetition(MCsamp_num, MCrep_num);
    MT_UB_mean = mean(MT_UB_list);
    MT_diff_list = MT_UB_list - MT_LB;
    MT_diff_mean = mean(MT_diff_list);

    MT_OTEB = MT.getMTErrorBoundBasedOnOT();
    MT_THEB = MT.getMTTheoreticalErrorBound(tolerance);

    fprintf(exp_log_file, ...
        ['trial %2d, marg num = %3d: time = %10.4e, ' ...
        'LB = %6.4f, UB = %6.4f, ' ...
        'diff = %8.6f, OTEB = %8.6f, THEB = %8.6f\n'], trial_id, ...
        marg_num, total_time, MT_LB, MT_UB_mean, MT_diff_mean, ...
        MT_OTEB, MT_THEB);

    fprintf(['trial %2d, marg num = %3d: time = %10.4e, ' ...
        'LB = %6.4f, UB = %6.4f, ' ...
        'diff = %8.6f, OTEB = %8.6f, THEB = %8.6f\n'], trial_id, ...
        marg_num, total_time, MT_LB, MT_UB_mean, MT_diff_mean, ...
        MT_OTEB, MT_THEB);

    save(sprintf(result_file_path, ...
        file_prefix, trial_id, marg_num), ...
        'output', 'total_time', ...
        'MCsamps', 'LSIP_primal', 'LSIP_dual', ...
        'LSIP_UB', 'LSIP_LB', ...
        'MT_UB_list', 'MT_UB_mean', 'MT_LB', ...
        'MT_diff_list', 'MT_diff_mean', ...
        'MT_OTEB', 'MT_THEB', '-v7.3');

    fprintf(exp_log_file, '--- experiment ends ---\n\n');
    fclose(exp_log_file);
end