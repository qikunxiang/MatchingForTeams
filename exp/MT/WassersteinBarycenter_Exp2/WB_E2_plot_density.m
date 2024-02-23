% Plot the density functions of the reference measure and the marginals
% from the same location-scatter family

CONFIG = WB_E2_config();

load(CONFIG.SAVEPATH_INPUTS);

plot_testfuncs = false;
plot_test_id = 3;

x_axis_lim = [-3.2, 3.2];
y_axis_lim = [-3.2, 3.2];

% instantiate the marginals
marg_cell = cell(marg_num, 1);
marg_dens_cell = cell(marg_num, 1);

for marg_id = 1:marg_num
    Meas = ProbMeas2D_MixNorm(marg_vertices_cell{marg_id}, ...
        marg_triangles_cell{marg_id}, ...
        marg_mixnorm_cell{marg_id});

    marg_dens_cell{marg_id} = reshape( ...
        Meas.densityFunction(quality_plot_hist_grid), ...
        quality_hist_x_num, quality_hist_y_num);

    marg_cell{marg_id} = Meas;
end

dens_max = max(max(vertcat(marg_dens_cell{:})));

figure('Position', [100, 100, 1000, 830]);
ha = tight_subplot(4, 5, [0.05, 0.016], [0.05, 0.008], [0.015, 0.06]);

for marg_id = 1:marg_num
    axes(ha(marg_id));
    hold on;
    plot_color = pcolor(quality_plot_hist_grid_x, ...
        quality_plot_hist_grid_y, ...
        marg_dens_cell{marg_id});
    plot_color.EdgeColor = 'interp';
    plot_color.FaceColor = 'interp';

    if plot_testfuncs
        Meas = marg_cell{marg_id};
        Meas.setSimplicialTestFuncs( ...
            marg_testfuncs_cell{plot_test_id}{marg_id}{:});
        triplot(Meas.SimplicialTestFuncs.Triangles, ...
            Meas.SimplicialTestFuncs.Vertices(:, 1), ...
            Meas.SimplicialTestFuncs.Vertices(:, 2));
    end

    box on;
    colormap('hot');

    clim([0, dens_max]);
    set(gca, 'Color', 'black');
    set(gca, 'XLim', x_axis_lim);
    set(gca, 'YLim', y_axis_lim);

    xlabel(sprintf('$\\mu_{%d}$', marg_id), ...
        'Interpreter', 'latex', 'FontSize', 14);
end


cb = colorbar('manual');
cb.Position = [0.955, 0.05, 0.015, 0.942];
