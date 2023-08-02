result_file_path = ['experiments/experiment2/' ...
    'exp2_results.mat'];

load(result_file_path);

for row_id = 1:length(marg_list)
    marg_str = sprintf('%3d & ', marg_list(row_id));

    rst_mean_err = mean(ParTrans.ERR(row_id, :), 'all');
    rst_max_err = max(ParTrans.ERR(row_id, :));
    rst_mean_time = round(mean(ParTrans.TIME(row_id, :), 'all'));
    rst_max_time = round(max(ParTrans.TIME(row_id, :)));
    rst_mean_sp = mean(ParTrans.SPAR(row_id, :), 'all');
    rst_max_sp = max(ParTrans.SPAR(row_id, :));
    rst_th_sp = max(ParTrans.THSPAR(row_id, :));

    if any(isnan([rst_mean_err; rst_max_err; rst_mean_time; ...
            rst_max_time; rst_mean_sp; rst_max_sp; rst_th_sp]))
        ParTrans_str = sprintf(['%8s & %8s & %6s & %6s & ' ...
            '%8s & %8s & %8s \\\\'], ...
            '--', '--', '--', '--', '--', '--', '--');
    else
        ParTrans_str = sprintf(['%8.3f & %8.3f & %6d & %6d & ' ...
            '%8.1f & %8d & %8d \\\\'], ...
            rst_mean_err * 1e4, rst_max_err * 1e4, ...
            rst_mean_time, rst_max_time, ...
            rst_mean_sp, rst_max_sp, rst_th_sp);
    end

    fprintf('%s%s\n', marg_str, ParTrans_str);
end