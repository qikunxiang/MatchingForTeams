CONFIG = BM1D_config();

trial_num = 10;
marg_list = [4; 6; 8; 10; 12; 14; 16; 18; 20; 50; 80; 100];

marg_testfunc_num_mat = zeros(trial_num, 100);
quality_testfunc_num_list = zeros(trial_num, 1);
for trial_id = 1:trial_num
    input_file_path = sprintf(CONFIG.SAVEPATH_INPUTS, ...
        trial_id);
    input_file = load(input_file_path, ...
        'testfuncs_cell', 'quality_testfuncs');

    for marg_id = 1:max(marg_list)
        marg_testfunc_num_mat(trial_id, marg_id) = ...
            length(input_file.testfuncs_cell{marg_id}{1});
    end

    quality_testfunc_num_list(trial_id) = ...
        size(input_file.quality_testfuncs{1}, 1);
end

marg_test_num = length(marg_list);

RST = struct;

RST.ERR = nan(marg_test_num, trial_num);
RST.RELERR = nan(marg_test_num, trial_num);
RST.THEB = nan(marg_test_num, trial_num);
RST.TIME = nan(marg_test_num, trial_num);
RST.LPTIME = nan(marg_test_num, trial_num);
RST.GLTIME = nan(marg_test_num, trial_num);
RST.ITER = nan(marg_test_num, trial_num);
RST.SPAR = nan(marg_test_num, trial_num);
RST.THSPAR = nan(marg_test_num, trial_num);

for marg_test_id = 1:marg_test_num
    marg_num = marg_list(marg_test_id);

    for trial_id = 1:trial_num
        result_file_path = sprintf(CONFIG.SAVEPATH_OUTPUTS, ...
            trial_id, marg_num);

        if exist(result_file_path, 'file')
            result_file = load(result_file_path, ...
                'MT_LB', 'LSIP_dual', 'MT_diff_mean', 'MT_THEB', 'output');
        
            RST.ERR(marg_test_id, trial_id) = ...
                result_file.MT_diff_mean;
            RST.RELERR(marg_test_id, trial_id) = ...
                result_file.MT_diff_mean / abs(result_file.MT_LB);
            RST.THEB(marg_test_id, trial_id) = ...
                result_file.MT_THEB;
            RST.TIME(marg_test_id, trial_id) = ...
                result_file.output.total_time;
            RST.LPTIME(marg_test_id, trial_id) = ...
                result_file.output.LP_time;
            RST.GLTIME(marg_test_id, trial_id) = ...
                result_file.output.global_time;
            RST.ITER(marg_test_id, trial_id) = ...
                result_file.output.iter;
            RST.SPAR(marg_test_id, trial_id) = ...
                length(result_file.LSIP_dual{1}.Probabilities);
            RST.THSPAR(marg_test_id, trial_id) = ...
                min(marg_testfunc_num_mat(trial_id, 1:marg_num) - 1) ...
                + quality_testfunc_num_list(trial_id) - 1 + 2;
        end
    end
end

save(CONFIG.SAVEPATH_SUMMARY, ...
    'RST', ...
    'marg_list', ...
    '-v7.3');
