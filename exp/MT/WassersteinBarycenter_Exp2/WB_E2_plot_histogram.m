% Plot the histogram of computed approximate Wasserstein barycenter

CONFIG = WB_E2_config();

load(CONFIG.SAVEPATH_INPUTS);
load(CONFIG.SAVEPATH_OUTPUTS);

x_axis_lim = [-1.6, 1.6];
y_axis_lim = [-1.6, 1.6];

x_label_cell = cell(test_num, 1);

for test_id = 1:test_num
    marg_testfuncs_num_sum = 0;

    for marg_id = 1:marg_num
        marg_testfuncs_num_sum = marg_testfuncs_num_sum ...
            + size(marg_testfuncs_cell{test_id}{marg_id}{1}, 1) - 1;
    end

    quality_testfuncs_num = size(quality_testfuncs_cell{test_id}{1}, 1) ...
        - 1;

    x_label_cell{test_id} = sprintf('(%d, %d)', ...
        marg_testfuncs_num_sum, quality_testfuncs_num);
end

dens_max = max(max(vertcat(WB_histpdf_cell{:})));

figure('Position', [0, 100, 1600, 265]);
ha = tight_subplot(1, test_num, [0, 0.02], ...
    [0.125, 0.025], [0.015, 0.005]);


for test_id = 1:test_num
    axes(ha(test_id));
    hold on;

    plot_color = pcolor(quality_plot_hist_grid_x, ...
        quality_plot_hist_grid_y, ...
        WB_histpdf_cell{test_id}');
    plot_color.EdgeColor = 'interp';
    plot_color.FaceColor = 'interp';

    box on;
    colormap('hot');

    clim([0, dens_max]);
    set(gca, 'Color', 'black');
    set(gca, 'XLim', x_axis_lim);
    set(gca, 'YLim', y_axis_lim);

    xlabel(x_label_cell{test_id}, 'FontWeight', 'bold');
end