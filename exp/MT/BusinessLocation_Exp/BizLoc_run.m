% Compute the lower and upper bounds for the business outlet locations
% problem

CONFIG = BizLoc_config();

load(CONFIG.SAVEPATH_INPUTS);

options = struct;
options.log_file = CONFIG.LOGPATH_LSIP_MAIN;
options.global_formualtion = 'DLOG';
options.discmeas_cand_index = marg_num;
options.display = true;
options.sanitation_threshold = 1e-8;
options.reduce = struct;
options.reduce.thres = 5;
options.reduce.max_iter = 2000;
options.reduce.freq = 100;

global_options = struct;
global_options.OutputFlag = 1;
global_options.LogToConsole = 1;
global_options.LogFile = CONFIG.LOGPATH_LSIP_GLOBAL;

LP_options = struct;
LP_options.OutputFlag = 1;
LP_options.LogToConsole = 1;
LP_options.LogFile = CONFIG.LOGPATH_LSIP_LP;

tolerance = 2e-4;

output_cell = cell(test_num, 1);
LSIP_primal_cell = cell(test_num, 1);
LSIP_dual_cell = cell(test_num, 1);
LSIP_LB_list = zeros(test_num, 1);
LSIP_UB_list = zeros(test_num, 1);

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


for test_id = 1:test_num

    if test_id == 1
        for marg_id = 1:marg_num
            marg_cell{marg_id}.setSimplicialTestFuncs( ...
                marg_testfuncs_cell{test_id}{marg_id}{:});
        end
    
        MT = MT2DBizLoc_ParTrans(marg_cell, costfuncs, ...
            quality_vertices, ...
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
        initial_constr{marg_id} = coup_cell{marg_id};
        initial_constr{marg_id}.min_val = -inf;
    end

    output = MT.run(initial_constr, tolerance);

    
    LSIP_primal = MT.Runtime.PrimalSolution;
    LSIP_dual = MT.Runtime.DualSolution;
    LSIP_LB = MT.Runtime.LSIP_LB;
    LSIP_UB = MT.Runtime.LSIP_UB;

    log_text = sprintf('Test %d: LSIP_LB = %10.4f, LSIP_UB = %10.4f\n', ...
        test_id, LSIP_LB, LSIP_UB);

    fprintf(main_log_file, log_text);
    fprintf(log_text);

    output_cell{test_id} = output;
    LSIP_primal_cell{test_id} = LSIP_primal;
    LSIP_dual_cell{test_id} = LSIP_dual;
    LSIP_LB_list(test_id) = LSIP_LB;
    LSIP_UB_list(test_id) = LSIP_UB;

    save(CONFIG.SAVEPATH_OUTPUTS, ...
        'output_cell', ...
        'LSIP_primal_cell', ...
        'LSIP_dual_cell', ...
        'LSIP_LB_list', ...
        'LSIP_UB_list', ...
        '-v7.3');
end

fprintf(main_log_file, '--- experiment ends ---\n\n');
fclose(main_log_file);