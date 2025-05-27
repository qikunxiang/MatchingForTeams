% Compute the lower and upper bounds for the 2-Wasserstein barycenter problem

CONFIG = WB_General_config();

load(CONFIG.SAVEPATH_INPUTS);

enable_quadratic_testfuncs = false;
tolerance = 2e-4;

options = struct;
options.log_file = CONFIG.LOGPATH_LSIP_MAIN;
options.display = true;
options.dual_prob_thres = 1e-8;
options.reduce = struct;
options.reduce.preserve_init_constr = true;
options.reduce.min_slack = 1e-7;
options.reduce.thres = 2e-3;
options.reduce.thres_quantile = 0.8;
options.reduce.max_iter = 5000;
options.reduce.freq = 300;


global_options = struct;
global_options.pool_size = 5;
global_options.boundary_buffer = 1e-3;
global_options.objective_threshold = -1e-6;
global_options.display = true;
global_options.log_file = CONFIG.LOGPATH_LSIP_GLOBAL;

LP_options = struct;
LP_options.Method = 2;
LP_options.Presolve = 0;
LP_options.Crossover = 0;
LP_options.OutputFlag = 1;
LP_options.LogToConsole = 1;
LP_options.LogFile = CONFIG.LOGPATH_LSIP_LP;

output_cell = cell(test_num, 1);
LSIP_primal_cell = cell(test_num, 1);
LSIP_dual_cell = cell(test_num, 1);
LSIP_LB_list = zeros(test_num, 1);
LSIP_UB_list = zeros(test_num, 1);
MT_LB_list = zeros(test_num, 1);

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

    if test_id == 1
        for marg_id = 1:marg_num
            marg_cell{marg_id}.setSimplicialTestFuncs(marg_testfuncs_cell{test_id}{marg_id}{:});
        end
    
        MT = MatchTeam2DWassersteinBarycenter( ...
            marg_cell, ...
            marg_weights, ...
            {quality_vertices}, ...
            quality_testfuncs_cell{test_id}, ...
            options, LP_options, global_options);
    
        coup_cell = MT.generateDiscreteCoupling();
    else
        coup_cell = MT.updateSimplicialTestFuncs(quality_testfuncs_cell{test_id}, marg_testfuncs_cell{test_id});
    end

    initial_constr = cell(marg_num, 1);

    for marg_id = 1:marg_num
        initial_constr{marg_id} = struct( ...
            'vertex_indices', coup_cell{marg_id}.vertex_indices, ...
            'points', coup_cell{marg_id}.points, ...
            'testfuncs_vals', coup_cell{marg_id}.testfuncs_vals, ...
            'costfunc_vals', coup_cell{marg_id}.costfunc_vals, ...
            'min_val', -inf);
    end

    output = MT.run(initial_constr, tolerance);

    LSIP_primal = MT.Runtime.PrimalSolution;
    LSIP_dual = MT.Runtime.DualSolution;
    LSIP_LB = MT.Runtime.LSIP_LB;
    LSIP_UB = MT.Runtime.LSIP_UB;

    MT_LB = MT.getMTLowerBound();

    log_text = sprintf('test %d: LB = %10.4f\n', ...
        test_id, ...
        MT_LB);

    fprintf(main_log_file, log_text);
    fprintf(log_text);

    output_cell{test_id} = output;
    LSIP_primal_cell{test_id} = LSIP_primal;
    LSIP_dual_cell{test_id} = LSIP_dual;
    LSIP_LB_list(test_id) = LSIP_LB;
    LSIP_UB_list(test_id) = LSIP_UB;
    MT_LB_list(test_id) = MT_LB;

    save(CONFIG.SAVEPATH_OUTPUTS_LB, ...
        'output_cell', ...
        'LSIP_primal_cell', ...
        'LSIP_dual_cell', ...
        'LSIP_LB_list', ...
        'LSIP_UB_list', ...
        'MT_LB_list', ...
        '-v7.3');
end

fprintf(main_log_file, '--- experiment ends ---\n\n');
fclose(main_log_file);