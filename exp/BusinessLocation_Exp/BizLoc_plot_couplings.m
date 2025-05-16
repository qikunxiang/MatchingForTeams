% Plot the computed couplings

CONFIG = BizLoc_config();

load(CONFIG.SAVEPATH_INPUTS);
load(CONFIG.SAVEPATH_COUPLINGS);

x_axis_lim = [-2.1, 2.1];
y_axis_lim = [-3.1, 2.1];

figure('Position', [100, 100, 1100, 280]);
ha = tight_subplot(1, 5, [0, 0.02], [0.135, 0.010], [0.015, 0.005]);

MCsamp_num = size(MCsamp.ContinuousQualities, 1);

for marg_id = 1:marg_num
    axes(ha(marg_id));
    hold on;
    
    for samp_id = 1:MCsamp_num
        x_i = MCsamp.ContinuousInputs{marg_id}(samp_id, :)';
        z = MCsamp.DiscreteQualities(samp_id, :)';

        scatter(x_i(1), x_i(2), 8, 'black', 'filled', 'MarkerEdgeAlpha', 0.2, 'MarkerFaceAlpha', 0.4);
        scatter(z(1), z(2), 8, 'red', 'filled', 'MarkerEdgeAlpha', 0.2, 'MarkerFaceAlpha', 0.4);

        line([x_i(1), z(1)], [x_i(2), z(2)], 'Color', [0, 0, 1, 0.1]);
    end

    box on;

    set(gca, 'XLim', x_axis_lim);
    set(gca, 'YLim', y_axis_lim);
    set(gca, 'YTick', -3:1:2);
    xlabel(sprintf('$\\tilde{\\gamma}_{%d}$', marg_id), 'Interpreter', 'latex', 'FontSize', 17);
end