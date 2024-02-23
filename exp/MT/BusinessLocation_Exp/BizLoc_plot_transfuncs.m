% Plot the computed business location distributions

CONFIG = BizLoc_config();

load(CONFIG.SAVEPATH_INPUTS);
load(CONFIG.SAVEPATH_OUTPUTS_TRANSFUNCS);

plot_index = 5;

figure('Position', [100, 100, 1300, 280]);
ha = tight_subplot(1, marg_num, [0, 0.015], ...
    [0.15, 0.02], [0.01, 0.05]);

tf_min = min(min(MT_transfuncs_cell{plot_index}));
tf_max = max(max(MT_transfuncs_cell{plot_index}));

for marg_id = 1:marg_num
    axes(ha(marg_id));

    pc = pcolor(quality_transfuncs_grid_x, ...
        quality_transfuncs_grid_y, ...
        reshape(MT_transfuncs_cell{plot_index}(:, marg_id), ...
        size(quality_transfuncs_grid_x, 1), ...
        size(quality_transfuncs_grid_x, 2)));

    pc.FaceColor = 'interp';
    pc.EdgeColor = 'interp';
    clim([tf_min, tf_max]);

    colormap('jet');

    set(gca, 'XLim', [quality_x_min, quality_x_max]);
    set(gca, 'YLim', [quality_y_min, quality_y_max]);
    set(gca, 'XTick', -2:1:2);
    set(gca, 'YTick', -3:1:2);

    xlabel(sprintf('$\\tilde{\\varphi}_{%d}$', marg_id), ...
        'Interpreter', 'latex', 'FontSize', 15);
end

cb = colorbar(ha(end), 'manual');
cb.Position = [0.96, 0.15, 0.015, 0.83];