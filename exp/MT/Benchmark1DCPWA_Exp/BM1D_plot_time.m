CONFIG = BM1D_config();

load(CONFIG.SAVEPATH_SUMMARY);

LP_time = mean(RST.LPTIME, 2, 'includemissing');
GL_time = mean(RST.GLTIME, 2, 'includemissing');

figure('Position', [100, 100, 500, 300]);
ha = tight_subplot(1, 1, [0, 0], [0.105, 0.055], [0.065, 0.02]);
axes(ha);
hold on;
plp = plot(marg_list, LP_time, ':^', 'Color', 'blue', ...
    'MarkerSize', 4);
pgl = plot(marg_list, GL_time, '--o', 'Color', 'red', ...
    'MarkerSize', 4);

box on;
set(gca, 'XLim', [0, max(marg_list)]);

legend([plp, pgl], ...
    {'LP', '$\mathtt{Oracle}$'}, ...
    'Interpreter', 'latex', 'Location', 'northwest');

xlabel('number $N$ of agent categories', 'Interpreter', 'latex');
ylabel('time (sec)');