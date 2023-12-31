input_file_path = ['experiments/experiment2/' ...
    'exp2_inputs_T%02d.mat'];
result_file_path = 'experiments/experiment2/%s_T%02d_M%03d.mat';
save_file_path = 'experiments/experiment2/exp2_results.mat';

result_file_prefix = 'results/RST_LOG_DLOG';

trial_num = 10;
marg_list = [4; 6; 8; 10; 12; 14; 16; 18; 20; 50; 80; 100];

marg_testfunc_num_mat = zeros(trial_num, 100);
quality_testfunc_num_list = zeros(trial_num, 1);
for trial_id = 1:trial_num
    input_file = load(sprintf(input_file_path, trial_id), ...
        'testfuncs_cell', 'quality_testfuncs');

    for marg_id = 1:max(marg_list)
        marg_testfunc_num_mat(trial_id, marg_id) = ...
            length(input_file.testfuncs_cell{marg_id}{1});
    end

    quality_testfunc_num_list(trial_id) = ...
        size(input_file.quality_testfuncs{1}, 1);
end

marg_test_num = length(marg_list);

ParTrans = struct;

ParTrans.ERR = nan(marg_test_num, trial_num);
ParTrans.RELERR = nan(marg_test_num, trial_num);
ParTrans.THEB = nan(marg_test_num, trial_num);
ParTrans.TIME = nan(marg_test_num, trial_num);
ParTrans.LPTIME = nan(marg_test_num, trial_num);
ParTrans.GLTIME = nan(marg_test_num, trial_num);
ParTrans.ITER = nan(marg_test_num, trial_num);
ParTrans.SPAR = nan(marg_test_num, trial_num);
ParTrans.THSPAR = nan(marg_test_num, trial_num);

for marg_test_id = 1:marg_test_num
    marg_num = marg_list(marg_test_id);

    for trial_id = 1:trial_num
        
        result_file_path_filled = sprintf(result_file_path, ...
            result_file_prefix, trial_id, marg_num);

        if exist(result_file_path_filled, 'file')
            result_file = load(result_file_path_filled, ...
                'MT_LB', 'LSIP_dual', 'MT_diff_mean', ...
                'MT_THEB', 'output');
        
            ParTrans.ERR(marg_test_id, trial_id) = ...
                result_file.MT_diff_mean;
            ParTrans.RELERR(marg_test_id, trial_id) = ...
                result_file.MT_diff_mean / abs(result_file.MT_LB);
            ParTrans.THEB(marg_test_id, trial_id) = ...
                result_file.MT_THEB;
            ParTrans.TIME(marg_test_id, trial_id) = ...
                result_file.output.total_time;
            ParTrans.LPTIME(marg_test_id, trial_id) = ...
                result_file.output.LP_time;
            ParTrans.GLTIME(marg_test_id, trial_id) = ...
                result_file.output.global_time;
            ParTrans.ITER(marg_test_id, trial_id) = ...
                result_file.output.iter;
            ParTrans.SPAR(marg_test_id, trial_id) = ...
                length(result_file.LSIP_dual{1}.Probabilities);
            ParTrans.THSPAR(marg_test_id, trial_id) = ...
                min(marg_testfunc_num_mat(trial_id, 1:marg_num) - 1) ...
                + quality_testfunc_num_list(trial_id) - 1 + 2;
        end
    end
end

save(save_file_path, 'ParTrans', 'marg_list');