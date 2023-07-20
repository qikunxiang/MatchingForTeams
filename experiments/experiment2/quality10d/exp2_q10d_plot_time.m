result_file_path = ['experiments/experiment2/quality10d/' ...
    'exp2_q10d_results.mat'];

load(result_file_path);

MMOT_LP_time = mean(MMOT.LPTIME, 2, 'includemissing');
MMOT_GL_time = mean(MMOT.GLTIME, 2, 'includemissing');

figure('Position', [100, 100, 500, 300]);
ha = tight_subplot(1, 1, [0, 0], [0.105, 0.055], [0.065, 0.02]);
axes(ha);
hold on;

mlp = plot(marg_list, MMOT_LP_time, ':^', 'Color', 'red', ...
    'MarkerSize', 4);
mgl = plot(marg_list, MMOT_GL_time, '--o', 'Color', 'red', ...
    'MarkerSize', 4);

box on;
set(gca, 'XLim', [0, max(marg_list)]);

legend([mlp, mgl], ...
    {'MMOT: LP', 'MMOT: $\mathtt{Oracle}_1$'}, ...
    'Interpreter', 'latex', 'Location', 'northwest');

xlabel('number $N$ of agent categories', 'Interpreter', 'latex');
ylabel('time (sec)');