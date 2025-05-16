% Plot the computed lower and upper bounds

CONFIG = BizLoc_config();

load(CONFIG.SAVEPATH_INPUTS);
load(CONFIG.SAVEPATH_OUTPUTS_UB);

MT_UB1_err = zeros(test_num, 2);
MT_UB2_err = zeros(test_num, 2);

for test_id = 1:test_num
    qq1 = quantile(MT_UB_cell{test_id, 1}, [0.025, 0.975]);
    qq2 = quantile(MT_UB_cell{test_id, 2}, [0.025, 0.975]);

    MT_UB1_err(test_id, 1) = MT_UB_mean_list(test_id, 1) - qq1(1);
    MT_UB1_err(test_id, 2) = qq1(2) - MT_UB_mean_list(test_id, 1);
    MT_UB2_err(test_id, 1) = MT_UB_mean_list(test_id, 2) - qq2(1);
    MT_UB2_err(test_id, 2) = qq2(2) - MT_UB_mean_list(test_id, 2);
end

decivar_num_list = zeros(test_num, 1);

for test_id = 1:test_num
    for marg_id = 1:marg_num
        decivar_num_list(test_id) = decivar_num_list(test_id) + size(marg_testfuncs_cell{test_id}{marg_id}{1}, 1) - 1;
    end

    quality_testfuncs_num = size(quality_testfuncs_cell{test_id}{1}, 1) - 1;
    decivar_num_list(test_id) = decivar_num_list(test_id) + marg_num * (quality_testfuncs_num + 1);
end

x_tick = [2e2, 3e2, 4e2, 5e2, 6e2, 7e2, 8e2, 9e2, 1e3, 2e3, 3e3, 4e3, 5e3, 6e3, 7e3, 8e3, 9e3, 1e4, 2e4, 3e4, 4e4];
x_tickangle = 70;

line_width = 1.25;

% figure of all lower and upper bounds

figure('Position', [100, 100, 420, 300]);
ha = tight_subplot(1, 1, [0, 0], [0.170, 0.020], [0.085, 0.025]);
axes(ha(1));

hold on;

handle_LB = plot(decivar_num_list, MT_LB_list, ...
    'Marker', 'o', 'Color', 'blue', 'LineStyle', '--', 'LineWidth', line_width);
handle_UB1 = plot(decivar_num_list, MT_UB_mean_list(:, 1), ...
    'Marker', 'x', 'Color', 'red', 'LineStyle', ':', 'LineWidth', line_width);
handle_UB2 = plot(decivar_num_list, MT_UB_mean_list(:, 2), ...
    'Marker', '+', 'Color', [214/255, 148/255, 0/255], 'LineStyle', ':', 'LineWidth', line_width);

box on;
grid on;

legend([handle_UB1, handle_UB2, handle_LB], ...
    {'$\hat{\alpha}_{\mathsf{MT}}^{\mathsf{UB}}$', ...
    '$\tilde{\alpha}_{\mathsf{MT}}^{\mathsf{UB}}$', ...
    '$\alpha_{\mathsf{MT}}^{\mathsf{LB}}$'}, ...
    'Location', 'northeast', ...
    'Interpreter', 'latex', ...
    'FontSize', 13);

xlabel('number of decision variables');
ylabel('objective');

set(gca, 'XScale', 'log');
set(gca, 'XLim', [min(x_tick), max(x_tick)]);
set(gca, 'XTick', x_tick);
xtickangle(x_tickangle);

% figure of sub-optimalities and error bounds

figure('Position', [520, 100, 420, 300]);
ha = tight_subplot(1, 1, [0, 0], [0.170, 0.020], [0.095, 0.025]);
axes(ha(1));

hold on;

handle_com1 = errorbar(decivar_num_list, MT_UB_mean_list(:, 1) - MT_LB_list, ...
    MT_UB1_err(:, 1), MT_UB1_err(:, 2), ...
    'Marker', 'none', 'Color', 'red', 'LineStyle', ':', 'LineWidth', line_width);
handle_com2 = errorbar(decivar_num_list, MT_UB_mean_list(:, 2) - MT_LB_list, ...
    MT_UB2_err(:, 1), MT_UB2_err(:, 2), ...
    'Marker', 'none', 'Color', [214, 148, 0] / 255, 'LineStyle', ':', 'LineWidth', line_width);
handle_th = plot(decivar_num_list, MT_THEB_list, ...
    'Marker', 'diamond', 'Color', [17, 17, 80] / 255, 'LineStyle', '-', 'LineWidth', line_width);

box on;
grid on;

legend([handle_com1, handle_com2, handle_th], ...
    {'$\hat{\epsilon}_{\mathsf{sub}}$', ...
    '$\tilde{\epsilon}_{\mathsf{sub}}$', ...
    '$\epsilon_{\mathsf{theo}}$'}, ...
    'Location', 'northeast', ...
    'Interpreter', 'latex', ...
    'FontSize', 13);

set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gca, 'XLim', [min(x_tick), max(x_tick)]);
set(gca, 'XTick', x_tick);
xtickangle(x_tickangle);
set(gca, 'YLim', [1e-3, 1e1]);

xlabel('number of decision variables');
ylabel('sub-optimality');

bound_conservative = MT_THEB_list ./ MT_diff_mean_list;
fprintf('The a priori upper bound is %d to %d times more conservative than the computed bound 1.\n', ...
    round(min(bound_conservative(:, 1))), round(max(bound_conservative(:, 1))));
fprintf('The a priori upper bound is %d to %d times more conservative than the computed bound 2.\n', ...
    round(min(bound_conservative(:, 2))), round(max(bound_conservative(:, 2))));