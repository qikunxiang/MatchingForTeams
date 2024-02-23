% Plot the computed lower and upper bounds

CONFIG = BizLoc_config();

load(CONFIG.SAVEPATH_INPUTS);
load(CONFIG.SAVEPATH_OUTPUTS_UB);

xx = 1:test_num;

MT_UB1_err = zeros(test_num, 2);
MT_UB2_err = zeros(test_num, 2);

for test_id = 1:test_num
    qq1 = quantile(MT_UB_cell{test_id, 1}, [0.025, 0.975]);
    qq2 = quantile(MT_UB_cell{test_id, 2}, [0.025, 0.975]);

    % if the quantiles are less than the lower bound, the error bar will
    % not be shown; in this case, we make the quantiles equal to the lower
    % bound plus 1e-5
    if qq1(1) <= MT_LB_list(test_id)
        qq1(1) = MT_LB_list(test_id) + 1e-5;
    end

    if qq2(1) <= MT_LB_list(test_id)
        qq2(1) = MT_LB_list(test_id) + 1e-5;
    end

    MT_UB1_err(test_id, 1) = MT_UB_mean_list(test_id, 1) - qq1(1);
    MT_UB1_err(test_id, 2) = qq1(2) - MT_UB_mean_list(test_id, 1);
    MT_UB2_err(test_id, 1) = MT_UB_mean_list(test_id, 2) - qq2(1);
    MT_UB2_err(test_id, 2) = qq2(2) - MT_UB_mean_list(test_id, 2);
end

tick_label_cell = cell(test_num, 1);

for test_id = 1:test_num
    marg_testfuncs_num_sum = 0;

    for marg_id = 1:marg_num
        marg_testfuncs_num_sum = marg_testfuncs_num_sum ...
            + size(marg_testfuncs_cell{test_id}{marg_id}{1}, 1) - 1;
    end

    quality_testfuncs_num = size(quality_testfuncs_cell{test_id}{1}, 1) ...
        - 1;

    tick_label_cell{test_id} = sprintf('(%d, %d)', ...
        marg_testfuncs_num_sum, quality_testfuncs_num);
end

line_width = 1.25;

% figure of all lower and upper bounds

figure('Position', [100, 100, 400, 400]);
ha = tight_subplot(1, 1, [0, 0], [0.155, 0.015], [0.085, 0.025]);
axes(ha(1));

hold on;

handle_LB = plot(xx, MT_LB_list, 'Marker', 'o', 'Color', 'blue', ...
    'LineStyle', '--', 'LineWidth', line_width);
handle_UB1 = plot(xx, MT_UB_mean_list(:, 1), ...
    'Marker', 'x', 'Color', 'red', 'LineStyle', ':', ...
    'LineWidth', line_width);
handle_UB2 = plot(xx, MT_UB_mean_list(:, 2), ...
    'Marker', '+', 'Color', [214/255, 148/255, 0/255], ...
    'LineStyle', ':', 'LineWidth', line_width);

box on;
grid on;

legend([handle_UB1, handle_UB2, handle_LB], ...
    {'$\hat{\alpha}_{\mathsf{MT}}^{\mathsf{UB}}$', ...
    '$\tilde{\alpha}_{\mathsf{MT}}^{\mathsf{UB}}$', ...
    '$\alpha_{\mathsf{MT}}^{\mathsf{LB}}$'}, ...
    'Location', 'northeast', ...
    'Interpreter', 'latex', ...
    'FontSize', 13);

set(gca, 'XTick', xx);
set(gca, 'XTickLabel', tick_label_cell);
xtickangle(25);

xlabel('number of test functions');
ylabel('objective');

% figure of sub-optimalities and error bounds

xx = 1:test_num;

figure('Position', [500, 100, 400, 400]);
ha = tight_subplot(1, 1, [0, 0], [0.155, 0.015], [0.095, 0.025]);
axes(ha(1));

hold on;

handle_com1 = errorbar(xx, MT_UB_mean_list(xx, 1) - MT_LB_list(xx), ...
    MT_UB1_err(xx, 1), MT_UB1_err(xx, 2), ...
    'Marker', 'none', 'Color', 'red', 'LineStyle', ':', ...
    'LineWidth', line_width);
handle_com2 = errorbar(xx, MT_UB_mean_list(xx, 2) - MT_LB_list(xx), ...
    MT_UB2_err(xx, 1), MT_UB2_err(xx, 2), ...
    'Marker', 'none', 'Color', [214, 148, 0] / 255, ...
    'LineStyle', ':', 'LineWidth', line_width);
handle_th = plot(xx, MT_THEB_list(xx), ...
    'Marker', 'diamond', 'Color', [17, 17, 80] / 255, ...
    'LineStyle', '-', 'LineWidth', line_width);

box on;
grid on;

legend([handle_com1, handle_com2, handle_th], ...
    {'$\hat{\epsilon}_{\mathsf{sub}}$', ...
    '$\tilde{\epsilon}_{\mathsf{sub}}$', ...
    '$\epsilon_{\mathsf{theo}}$'}, ...
    'Location', 'northeast', ...
    'Interpreter', 'latex', ...
    'FontSize', 13);

set(gca, 'XTick', xx);
set(gca, 'XTickLabel', tick_label_cell);
xtickangle(25);
set(gca, 'YScale', 'log');
set(gca, 'YLim', [1e-3, 1e1]);

xlabel('number of test functions');
ylabel('sub-optimality');
