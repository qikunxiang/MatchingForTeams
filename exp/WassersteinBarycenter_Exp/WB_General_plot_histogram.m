% Plot the histogram of computed approximate Wasserstein barycenter

CONFIG = WB_General_config();

load(CONFIG.SAVEPATH_INPUTS);
load(CONFIG.SAVEPATH_OUTPUTS_UB);

x_axis_lim = [-1.6, 1.6];
y_axis_lim = [-1.6, 1.6];

decivar_num_list = zeros(test_num, 1);

for test_id = 1:test_num
    for marg_id = 1:marg_num
        decivar_num_list(test_id) = decivar_num_list(test_id) + size(marg_testfuncs_cell{test_id}{marg_id}{1}, 1) - 1;
    end

    quality_testfuncs_num = size(quality_testfuncs_cell{test_id}{1}, 1) - 1;
    decivar_num_list(test_id) = decivar_num_list(test_id) + marg_num * (quality_testfuncs_num + 1);
end


figure('Position', [0, 100, 1500, 255]);
ha = tight_subplot(1, test_num, [0, 0.018], [0.125, 0.010], [0.013, 0.0035]);


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

    dens_max = max(max(vertcat(WB_histpdf_cell{test_id})));
    clim([0, dens_max]);
    set(gca, 'Color', 'black');
    set(gca, 'XLim', x_axis_lim);
    set(gca, 'YLim', y_axis_lim);

    xlabel(sprintf('$n=%d$', decivar_num_list(test_id)), 'Interpreter', 'latex', 'FontSize', 15);
end