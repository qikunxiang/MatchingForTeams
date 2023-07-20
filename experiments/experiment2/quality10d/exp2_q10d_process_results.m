input_file_path = ['experiments/experiment2/quality10d/' ...
    'exp2_q10d_inputs_T%02d.mat'];
result_file_path = 'experiments/experiment2/quality10d/%s_T%02d_M%03d.mat';
save_file_path = ['experiments/experiment2/quality10d/' ...
    'exp2_q10d_results.mat'];

MMOT_result_file_prefix = 'MMOT_results/RST_INC';

trial_num = 10;
marg_list = [4; 6; 8; 10; 12; 14];

marg_testfunc_num_mat = zeros(trial_num, 100);
for trial_id = 1:trial_num
    input_file = load(sprintf(input_file_path, trial_id), ...
        'testfuncs_cell');

    for marg_id = 1:max(marg_list)
        marg_testfunc_num_mat(trial_id, marg_id) = ...
            length(input_file.testfuncs_cell{marg_id}{1});
    end
end

marg_test_num = length(marg_list);

MMOT = struct;

MMOT.ERR = nan(marg_test_num, trial_num);
MMOT.RELERR = nan(marg_test_num, trial_num);
MMOT.THEB = nan(marg_test_num, trial_num);
MMOT.TIME = nan(marg_test_num, trial_num);
MMOT.LPTIME = nan(marg_test_num, trial_num);
MMOT.GLTIME = nan(marg_test_num, trial_num);
MMOT.ITER = nan(marg_test_num, trial_num);
MMOT.SPAR = nan(marg_test_num, trial_num);
MMOT.THSPAR = nan(marg_test_num, trial_num);

for marg_test_id = 1:marg_test_num
    marg_num = marg_list(marg_test_id);

    for trial_id = 1:trial_num
        
        result_file_path_filled = sprintf(result_file_path, ...
            MMOT_result_file_prefix, trial_id, marg_num);

        if exist(result_file_path_filled, 'file')
            result_file = load(result_file_path_filled, ...
                'MT_LB', 'LSIP_dual', 'MT_diff_mean', 'MT_THEB', 'output');
        
            MMOT.ERR(marg_test_id, trial_id) = result_file.MT_diff_mean;
            MMOT.RELERR(marg_test_id, trial_id) = ...
                result_file.MT_diff_mean / abs(result_file.MT_LB);
            MMOT.THEB(marg_test_id, trial_id) = ...
                result_file.MT_THEB;
            MMOT.TIME(marg_test_id, trial_id) = ...
                result_file.output.total_time;
            MMOT.LPTIME(marg_test_id, trial_id) = ...
                result_file.output.LP_time;
            MMOT.GLTIME(marg_test_id, trial_id) = ...
                result_file.output.global_time;
            MMOT.ITER(marg_test_id, trial_id) = result_file.output.iter;
            MMOT.SPAR(marg_test_id, trial_id) = ...
                length(result_file.LSIP_dual.Probabilities);
            MMOT.THSPAR(marg_test_id, trial_id) = ...
                sum(marg_testfunc_num_mat(trial_id, 1:marg_num) - 1) + 2;
        end
    end
end

save(save_file_path, 'MMOT', 'marg_list');