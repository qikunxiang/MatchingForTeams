% Plot the density functions of the marginals

CONFIG = BizLoc_config();

load(CONFIG.SAVEPATH_INPUTS);

plot_testfuncs = false;
plot_test_id = 3;

% instantiate the marginals
marg_cell = cell(marg_num, 1);
marg_dens_cell = cell(marg_num, 1);

dens_min = inf;
dens_max = 0;

for marg_id = 1:marg_num
    Meas = ProbMeas2D_CPWADens(marg_vertices_cell{marg_id}, ...
        marg_triangles_cell{marg_id}, ...
        marg_density_cell{marg_id});
    marg_cell{marg_id} = Meas;

    marg_dens_cell{marg_id} = reshape( ...
        Meas.densityFunction(marg_plot_grid_cell{marg_id}), ...
        plot_y_num, plot_x_num);

    dens_min = min(dens_min, min(min(marg_dens_cell{marg_id})));
    dens_max = max(dens_max, max(max(marg_dens_cell{marg_id})));
end

x_axis_lim = [-2.1, 2.1];
y_axis_lim = [-3.1, 2.1];


figure('Position', [100, 100, 1250, 280]);
ha = tight_subplot(1, 6, [0, 0.015], [0.12, 0.02], [0.015, 0.04]);

axes(ha(1));
hold on;
box on;

city = polyshape([quality_x_min, quality_y_min; ...
    quality_x_max, quality_y_min; ...
    quality_x_max, quality_y_max; ...
    quality_x_min, quality_y_max]);

plot(city, 'FaceColor', '#4DBEEE', 'EdgeColor', '#0072BD', ...
    'LineWidth', 5);

rad = 0.15;
ang = linspace(0, 2 * pi, 200)';

for sta_id = 1:size(costfuncs.stations)
    line(costfuncs.stations(sta_id, 1) + cos(ang) * rad, ...
        costfuncs.stations(sta_id, 2) + sin(ang) * rad, ...
        'Color', '#77AC30', 'LineWidth', 3);
end

xx = linspace(0, 1, 100)';

line(xx * 1.2 - 1.35, xx * 0 - 1.5, 'Color', '#77AC30', 'LineWidth', 2);
line(xx * 0, xx * 1.2 - 1.35, 'Color', '#77AC30', 'LineWidth', 2);
line(xx * 0, xx * 1.2 + 0.15, 'Color', '#77AC30', 'LineWidth', 2);
line(xx * 1.2 + 0.15, xx * 0 + 1.5, 'Color', '#77AC30', 'LineWidth', 2);

arr_ps_l = polyshape([1.5, -1.1; 1.3, -1.6; 1.5, -1.5]);
arr_ps_r = polyshape([1.5, -1.1; 1.7, -1.6; 1.5, -1.5]);
plot(arr_ps_l, 'EdgeColor', '#77AC30', 'FaceColor', 'none', ...
    'LineWidth', 2);
plot(arr_ps_r, 'EdgeColor', '#77AC30', 'FaceColor', '#77AC30', ...
    'LineWidth', 2);
text(1.36, -0.9, 'N', 'FontWeight', 'bold', 'Color', '#77AC30', ...
    'FontSize', 16);

set(gca, 'XLim', x_axis_lim);
set(gca, 'YLim', y_axis_lim);
set(gca, 'YTick', -3:1:2);
xlabel('train stations', 'FontWeight', 'bold');

for marg_id = 1:marg_num
    axes(ha(marg_id + 1));
    hold on;
    plot_color = pcolor(marg_plot_grid_x_cell{marg_id}, ...
        marg_plot_grid_y_cell{marg_id}, ...
        marg_dens_cell{marg_id});
    plot_color.EdgeColor = 'interp';
    plot_color.FaceColor = 'interp';
    clim([dens_min, dens_max]);

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

    set(gca, 'XLim', x_axis_lim);
    set(gca, 'YLim', y_axis_lim);
    set(gca, 'YTick', -3:1:2);
    xlabel(sprintf('$\\mu_{%d}$', marg_id), ...
        'Interpreter', 'latex', 'FontSize', 14);
end

cb = colorbar(ha(end), 'manual');
cb.Position = [0.964, 0.12, 0.015, 0.86];