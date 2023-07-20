result_file_path = ['experiments/experiment2/quality2d/' ...
    'exp2_q2d_results.mat'];

load(result_file_path);

MMOT_LP_time = mean(MMOT.LPTIME, 2, 'includemissing');
MMOT_GL_time = mean(MMOT.GLTIME, 2, 'includemissing');
ParTrans_LP_time = mean(ParTrans.LPTIME, 2, 'includemissing');
ParTrans_GL_time = mean(ParTrans.GLTIME, 2, 'includemissing');

figure('Position', [100, 100, 500, 300]);
ha = tight_subplot(1, 1, [0, 0], [0.105, 0.055], [0.065, 0.02]);
axes(ha);
hold on;

mlp = plot(marg_list, MMOT_LP_time, ':^', 'Color', 'red', ...
    'MarkerSize', 4);
mgl = plot(marg_list, MMOT_GL_time, '--o', 'Color', 'red', ...
    'MarkerSize', 4);
plp = plot(marg_list, ParTrans_LP_time, ':^', 'Color', 'blue', ...
    'MarkerSize', 4);
pgl = plot(marg_list, ParTrans_GL_time, '--o', 'Color', 'blue', ...
    'MarkerSize', 4);

box on;
set(gca, 'XLim', [0, max(marg_list)]);

legend([mlp, mgl, plp, pgl], ...
    {'MMOT: LP', 'MMOT: $\mathtt{Oracle}_1$', ...
    'Parametric: LP', 'Parametric: $\mathtt{Oracle}_2$'}, ...
    'Interpreter', 'latex', 'Location', 'northwest');

xlabel('number $N$ of agent categories', 'Interpreter', 'latex');
ylabel('time (sec)');


% magnification
plot_list = 1:8;

figure('Position', [600, 100, 500, 300]);
ha = tight_subplot(1, 1, [0, 0], [0.105, 0.055], [0.085, 0.02]);
axes(ha);
hold on;

mlp = plot(marg_list(plot_list), MMOT_LP_time(plot_list), ...
    ':^', 'Color', 'red', 'MarkerSize', 4);
mgl = plot(marg_list(plot_list), MMOT_GL_time(plot_list), ...
    '--o', 'Color', 'red', 'MarkerSize', 4);
plp = plot(marg_list(plot_list), ParTrans_LP_time(plot_list), ...
    ':^', 'Color', 'blue', 'MarkerSize', 4);
pgl = plot(marg_list(plot_list), ParTrans_GL_time(plot_list), ...
    '--o', 'Color', 'blue', 'MarkerSize', 4);

box on;
set(gca, 'XLim', [min(marg_list(plot_list)), max(marg_list(plot_list))]);

legend([mlp, mgl, plp, pgl], ...
    {'MMOT: LP', 'MMOT: $\mathtt{Oracle}_1$', ...
    'Parametric: LP', 'Parametric: $\mathtt{Oracle}_2$'}, ...
    'Interpreter', 'latex', 'Location', 'northwest');

xlabel('number $N$ of agent categories', 'Interpreter', 'latex');
ylabel('time (sec)');