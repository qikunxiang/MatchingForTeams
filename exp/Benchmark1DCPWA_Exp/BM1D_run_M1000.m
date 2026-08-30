CONFIG = BM1D_config();

options = struct;
options.log_file = CONFIG.LOGPATH_LSIP_MAIN;
options.global_parfor = true;
options.global_formulation = 'LOG_DLOG';
options.dual_prob_thres = 1e-7;
options.new_cut_thres = 1e-8;
options.sanitation_threshold = 1e-8;
options.time_limit = inf;
options.display = true;
options.reduce = struct;
options.reduce.thres = inf;
options.reduce.max_iter = 0;
options.reduce.freq = inf;

global_options = struct;
global_options.TimeLimit = 60000;
global_options.BestObjStop = -inf;
global_options.PoolSearchMode = 1;
global_options.PoolSolutions = 1;
global_options.OutputFlag = 0;
global_options.LogToConsole = 0;

LP_options = struct;
LP_options.Method = 2;
LP_options.Presolve = 0;
LP_options.Crossover = 0;
LP_options.OutputFlag = 0;
LP_options.LogToConsole = 0;

tolerance = 5e-5;
MCsamp_num = 1e4;
MCrep_num = 1000;

exp_log_file_path = CONFIG.LOGPATH_MAIN;
randtrial_num = 10;
marg_num = 1000;

for trial_id = 1:randtrial_num

    rng(10000 + trial_id * 1000 + marg_num, 'combRecursive');

    load(sprintf(CONFIG.SAVEPATH_INPUTS, trial_id));

    exp_log_file = fopen(exp_log_file_path, 'a');

    if exp_log_file < 0
        error('cannot open log file');
    end

    fprintf(exp_log_file, '--- experiment starts ---\n');

    marg_cell = cell(marg_num, 1);
    scaled_costfunc_cell = cell(marg_num, 1);

    for marg_id = 1:marg_num
        marg_cell{marg_id} = ProbMeas1DCPWADens(marg_knots_cell{marg_id}, marg_dens_cell{marg_id});

        % scale the cost functions such that the theoretical error bounds stay the same
        scaled_costfunc_cell{marg_id} = costfunc_cell{marg_id};
        scaled_costfunc_cell{marg_id}.values = scaled_costfunc_cell{marg_id}.values / marg_num;
    end

    total_timer = tic;
    MT = MatchTeam1DCPWA(marg_cell, scaled_costfunc_cell, quality_vertices, [], options, LP_options, global_options);
    MT.setSimplicialTestFuncs(quality_testfuncs, testfuncs_cell(1:marg_num));
    coup_cell = MT.generateDiscreteCoupling();

    initial_constr = cell(marg_num, 1);

    for marg_id = 1:marg_num
        initial_constr{marg_id} = struct( ...
            'inputs', coup_cell{marg_id}.inputs, ...
            'qualities', coup_cell{marg_id}.qualities, ...
            'marg_testfunc_vals', coup_cell{marg_id}.marg_testfunc_vals, ...
            'quality_testfunc_vals', coup_cell{marg_id}.quality_testfunc_vals, ...
            'costfunc_vals', coup_cell{marg_id}.costfunc_vals, ...
            'min_val', -inf);
    end

    parpool(10);

    output = MT.run(initial_constr, tolerance);
    
    delete(gcp('nocreate'));

    if output.time_limit_exceeded
        fprintf(exp_log_file, 'trial %2d, marg num = %3d: time limit exceeded\n', trial_id, marg_num);

        fprintf('trial %2d, marg num = %3d: time limit exceeded\n', trial_id, marg_num);
        fprintf(exp_log_file, '--- experiment ends ---\n\n');
        fclose(exp_log_file);
        continue;
    end

    LSIP_primal = MT.Runtime.PrimalSolution;
    LSIP_dual = MT.Runtime.DualSolution;
    LSIP_LB = MT.Runtime.LSIP_LB;
    LSIP_UB = MT.Runtime.LSIP_UB;

    MT_LB = MT.getMTLowerBound();
    MT_sparsity = length(LSIP_dual{1}.Probabilities);

    [MT_UB_list, MCsamps] = MT.getMTUpperBoundDiscreteWRepetition(MCsamp_num, MCrep_num, ...
        RandStream('mrg32k3a', 'Seed', 10001 + trial_id * 1000 + marg_num), 1e3);
    MT_UB_mean = mean(MT_UB_list);
    MT_diff_list = MT_UB_list - MT_LB;
    MT_diff_mean = mean(MT_diff_list);

    MT_OTEB = MT.getMTErrorBoundBasedOnOT();
    MT_THEB = MT.getMTTheoreticalErrorBound(tolerance);

    total_time = toc(total_timer);

    fprintf(exp_log_file, ...
        ['trial %2d, marg num = %3d: time = %10.4e, ' ...
        'LB = %6.4f, UB = %6.4f, ' ...
        'diff = %8.6f, OTEB = %8.6f, THEB = %8.6f, sparsity = %5d\n'], trial_id, ...
        marg_num, total_time, MT_LB, MT_UB_mean, MT_diff_mean, ...
        MT_OTEB, MT_THEB, MT_sparsity);

    fprintf(['trial %2d, marg num = %3d: time = %10.4e, ' ...
        'LB = %6.4f, UB = %6.4f, ' ...
        'diff = %8.6f, OTEB = %8.6f, THEB = %8.6f, sparsity = %5d\n'], trial_id, ...
        marg_num, total_time, MT_LB, MT_UB_mean, MT_diff_mean, ...
        MT_OTEB, MT_THEB, MT_sparsity);

    save(sprintf(CONFIG.SAVEPATH_OUTPUTS, trial_id, marg_num), ...
        'output', 'total_time', ...
        'LSIP_UB', 'LSIP_LB', ...
        'MT_UB_list', 'MT_UB_mean', 'MT_LB', ...
        'MT_diff_list', 'MT_diff_mean', ...
        'MT_OTEB', 'MT_THEB', 'MT_sparsity', '-v7.3');

    fprintf(exp_log_file, '--- experiment ends ---\n\n');
    fclose(exp_log_file);
end