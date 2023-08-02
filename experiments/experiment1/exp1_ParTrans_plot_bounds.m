rst_file_path = 'experiments/experiment1/exp1_ParTrans_rst.mat';
input_file_path = 'experiments/experiment1/exp1_inputs.mat';

load(rst_file_path, 'MT_LB_list', ...
    'MT_UB_cell', 'MT_UB_mean_list', ...
    'MT_diff_cell', 'MT_diff_mean_list', ...
    'MT_OTEB_list', 'MT_THEB_list');

load(input_file_path, 'testfunc_num_mat', 'quality_testfuncs_num_list');

test_num = length(MT_LB_list);

xx = 1:test_num;

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

tick_label_cell = cell(test_num, 1);

for test_id = 1:test_num
    tick_label_cell{test_id} = sprintf('(%d, %d)', ...
        sum(testfunc_num_mat(test_id, :) - 1), ...
        quality_testfuncs_num_list(test_id) - 1);
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
xtickangle(45);
set(gca, 'YLim', [-13, -5]);

xlabel('number of test functions');
ylabel('objective');


% figure of lower and upper bounds in settings 4 to 7

xx = 4:test_num;

figure('Position', [500, 100, 400, 400]);
ha = tight_subplot(1, 1, [0, 0], [0.155, 0.015], [0.11, 0.025]);
axes(ha(1));

hold on;

handle_LB = plot(xx, MT_LB_list(xx), 'Marker', 'o', 'Color', 'blue', ...
    'LineStyle', '--', 'LineWidth', line_width);
handle_UB1 = errorbar(xx, MT_UB_mean_list(xx, 1), ...
    MT_UB1_err(xx, 1), MT_UB1_err(xx, 2), ...
    'Marker', 'none', 'Color', 'red', 'LineStyle', ':', ...
    'LineWidth', line_width);
handle_UB2 = errorbar(xx, MT_UB_mean_list(xx, 2), ...
    MT_UB2_err(xx, 1), MT_UB2_err(xx, 2), ...
    'Marker', 'none', 'Color', [214, 148, 0] / 255, ...
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
set(gca, 'XTickLabel', tick_label_cell(xx));
xtickangle(45);
set(gca, 'YLim', [-8.8, -8.5]);

xlabel('number of test functions');
ylabel('objective');


% figure of sub-optimalities and error bounds

xx = 1:test_num;

figure('Position', [900, 100, 400, 400]);
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
    {'$\hat{\epsilon}_{\mathsf{MT}}$', ...
    '$\tilde{\epsilon}_{\mathsf{MT}}$', ...
    '$\epsilon^\ddagger_{\overline{\overline{W}}_1}$'}, ...
    'Location', 'northeast', ...
    'Interpreter', 'latex', ...
    'FontSize', 13);

set(gca, 'XTick', xx);
set(gca, 'XTickLabel', tick_label_cell);
xtickangle(45);
set(gca, 'YScale', 'log');
set(gca, 'YLim', [1e-2, 1e3]);

xlabel('number of test functions');
ylabel('sub-optimality');
