CONFIG = BM1D_config();

load(CONFIG.SAVEPATH_SUMMARY);

for row_id = 1:length(marg_list)
    marg_str = sprintf('%3d & ', marg_list(row_id));

    rst_mean_err = mean(RST.ERR(row_id, :), 'all');
    rst_max_err = max(RST.ERR(row_id, :));
    rst_mean_LP_time = mean(RST.LPTIME(row_id, :), 'all') / marg_list(row_id);
    rst_max_LP_time = max(RST.LPTIME(row_id, :)) / marg_list(row_id);
    rst_mean_global_time = mean(RST.GLTIME(row_id, :), 'all') / marg_list(row_id);
    rst_max_global_time = max(RST.GLTIME(row_id, :)) / marg_list(row_id);
    rst_mean_total_time = mean(RST.TIME(row_id, :), 'all') / marg_list(row_id);
    rst_max_total_time = max(RST.TIME(row_id, :)) / marg_list(row_id);
    rst_mean_sp = mean(RST.SPAR(row_id, :), 'all');
    rst_max_sp = max(RST.SPAR(row_id, :));
    rst_th_sp = max(RST.THSPAR(row_id, :));

    if any(isnan([rst_mean_err; rst_max_err; rst_mean_LP_time; rst_max_LP_time; rst_mean_global_time; rst_max_global_time; ...
            rst_mean_total_time; rst_max_total_time; rst_mean_sp; rst_max_sp; rst_th_sp]))
        str = sprintf('%8s & %8s & %6s & %6s & %6s & %6s & %6s & %6s \\\\', ...
            '--', '--', '--', '--', '--', '--', '--', '--');
    else
        str = sprintf('%8.3f & %8.3f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f \\\\', ...
            rst_mean_err * 1e4, ...
            rst_max_err * 1e4, ...
            rst_mean_LP_time, ...
            rst_max_LP_time, ...
            rst_mean_global_time, ...
            rst_max_global_time, ...
            rst_mean_total_time, ...
            rst_max_total_time);
    end

    fprintf('%s%s\n', marg_str, str);
end

fprintf('\n\n');

for row_id = 1:length(marg_list)
    marg_str = sprintf('%3d & ', marg_list(row_id));

    rst_mean_sp = mean(RST.SPAR(row_id, :), 'all');
    rst_max_sp = max(RST.SPAR(row_id, :));
    rst_th_sp = max(RST.THSPAR(row_id, :));

    if any(isnan([rst_mean_sp; rst_max_sp; rst_th_sp]))
        str = sprintf('%8s & %8s & %8s \\\\', '--', '--', '--');
    else
        str = sprintf('%8.1f & %8d & %8d \\\\', ...
            rst_mean_sp, ...
            rst_max_sp, ...
            rst_th_sp);
    end

    fprintf('%s%s\n', marg_str, str);
end