result_file_path = ['experiments/experiment2/quality2d/' ...
    'exp2_q2d_results.mat'];

load(result_file_path);

fprintf('--------------  First table begins ----------------\n');

for row_id = 1:length(marg_list)
    marg_str = sprintf('%3d & ', marg_list(row_id));

    rst_mean_err = mean(MMOT.ERR(row_id, :), 'all');
    rst_max_err = max(MMOT.ERR(row_id, :));
    rst_mean_time = round(mean(MMOT.TIME(row_id, :), 'all'));
    rst_max_time = round(max(MMOT.TIME(row_id, :)));

    if any(isnan([rst_mean_err; rst_max_err; rst_mean_time; rst_max_time]))
        MMOT_str = sprintf('%8s & %8s & %6s & %6s & ', ...
            '--', '--', '--', '--');
    else
        MMOT_str = sprintf('%8.3f & %8.3f & %6d & %6d & ', ...
            rst_mean_err * 1e4, rst_max_err * 1e4, ...
            rst_mean_time, rst_max_time);
    end

    rst_mean_err = mean(ParTrans.ERR(row_id, :), 'all');
    rst_max_err = max(ParTrans.ERR(row_id, :));
    rst_mean_time = round(mean(ParTrans.TIME(row_id, :), 'all'));
    rst_max_time = round(max(ParTrans.TIME(row_id, :)));

    if any(isnan([rst_mean_err; rst_max_err; rst_mean_time; rst_max_time]))
        ParTrans_str = sprintf('%8s & %8s & %6s & %6s \\\\', ...
            '--', '--', '--', '--');
    else
        ParTrans_str = sprintf('%8.3f & %8.3f & %6d & %6d \\\\', ...
            rst_mean_err * 1e4, rst_max_err * 1e4, ...
            rst_mean_time, rst_max_time);
    end

    fprintf('%s%s%s\n', marg_str, MMOT_str, ParTrans_str);
end


fprintf('---------------  First table ends -----------------\n\n');


fprintf('--------------  Second table begins ---------------\n');


for row_id = 1:length(marg_list)
    marg_str = sprintf('%3d & ', marg_list(row_id));

    rst_mean_sp = mean(MMOT.SPAR(row_id, :), 'all');
    rst_max_sp = max(MMOT.SPAR(row_id, :));
    rst_th_sp = max(MMOT.THSPAR(row_id, :));

    if any(isnan([rst_mean_sp; rst_max_sp]))
        MMOT_str = sprintf('%8s & %8s & %8s & ', '--', '--', '--');
    else
        MMOT_str = sprintf('%8.1f & %8d & %8d & ', ...
            rst_mean_sp, rst_max_sp, rst_th_sp);
    end

    rst_mean_sp = mean(ParTrans.SPAR(row_id, :), 'all');
    rst_max_sp = max(ParTrans.SPAR(row_id, :));
    rst_th_sp = max(ParTrans.THSPAR(row_id, :));

    if any(isnan([rst_mean_sp; rst_max_sp]))
        ParTrans_str = sprintf('%8s & %8s & %8s \\\\', '--', '--', '--');
    else
        ParTrans_str = sprintf('%8.1f & %8d & %8d \\\\', ...
            rst_mean_sp, rst_max_sp, rst_th_sp);
    end

    fprintf('%s%s%s\n', marg_str, MMOT_str, ParTrans_str);
end


fprintf('---------------  Second table ends ----------------\n');