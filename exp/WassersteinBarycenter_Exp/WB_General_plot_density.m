% Plot the density functions of the reference measure and the marginals from the same location-scatter family

CONFIG = WB_General_config();

load(CONFIG.SAVEPATH_INPUTS);

plot_testfuncs = false;
plot_test_id = 1;

x_axis_lim = [-2.45, 2.45];
y_axis_lim = [-2.45, 2.45];

marg_plot_x_num = 300;
marg_plot_y_num = 300;

marg_plot_x_min = -2.45;
marg_plot_x_max = 2.45;
marg_plot_y_min = -2.45;
marg_plot_y_max = 2.45;

marg_plot_edge_x = linspace(marg_plot_x_min, marg_plot_x_max, marg_plot_x_num + 1)';
marg_plot_edge_y = linspace(marg_plot_y_min, marg_plot_y_max, marg_plot_y_num + 1)';
marg_plot_x = (marg_plot_edge_x(1:end - 1) + marg_plot_edge_x(2:end)) / 2;
marg_plot_y = (marg_plot_edge_y(1:end - 1) + marg_plot_edge_y(2:end)) / 2;
[marg_plot_grid_x, marg_plot_grid_y] = meshgrid(marg_plot_x, marg_plot_y);
marg_plot_grid = [marg_plot_grid_x(:), marg_plot_grid_y(:)];

% instantiate the marginals
marg_cell = cell(marg_num, 1);
marg_dens_cell = cell(marg_num, 1);

for marg_id = 1:marg_num
    Meas = ProbMeas2DCPWADens( ...
        marg_vertices_cell{marg_id}, ...
        marg_triangles_cell{marg_id}, ...
        marg_densities_cell{marg_id});

    marg_dens_cell{marg_id} = reshape(Meas.densityFunction(marg_plot_grid), marg_plot_x_num, marg_plot_y_num);

    marg_cell{marg_id} = Meas;
end

dens_max = max(max(vertcat(marg_dens_cell{:})));

figure('Position', [100, 100, 1600, 310]);
ha = tight_subplot(2, 10, [0.010, 0.0025], [0.0135, 0.008], [0.0025, 0.03]);

for marg_id = 1:marg_num
    axes(ha(marg_id));
    hold on;
    plot_color = pcolor(marg_plot_grid_x, marg_plot_grid_y, marg_dens_cell{marg_id});
    plot_color.EdgeColor = 'interp';
    plot_color.FaceColor = 'interp';

    if plot_testfuncs
        Meas = marg_cell{marg_id};
        Meas.setSimplicialTestFuncs(marg_testfuncs_cell{plot_test_id}{marg_id}{:});
        triplot(Meas.SimplicialTestFuncs.Triangles, ...
            Meas.SimplicialTestFuncs.Vertices(:, 1), Meas.SimplicialTestFuncs.Vertices(:, 2), ...
            'Color', 'yellow');
    end

    box on;
    colormap('hot');

    clim([0, dens_max * 0.8]);
    set(gca, 'Color', 'black');
    set(gca, 'XLim', x_axis_lim);
    set(gca, 'YLim', y_axis_lim);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);

    % xlabel(sprintf('$\\mu_{%d}$', marg_id), 'Interpreter', 'latex', 'FontSize', 14);
end

cb = colorbar('manual');
cb.Position = [0.973, 0.0135, 0.010, 0.978];
