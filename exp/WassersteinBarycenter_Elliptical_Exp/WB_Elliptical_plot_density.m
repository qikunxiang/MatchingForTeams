% Plot the density functions of the reference measure and the marginals from the same elliptical family

CONFIG = WB_Elliptical_config();

load(CONFIG.SAVEPATH_INPUTS);

plot_testfuncs = false;
plot_test_id = 3;

x_axis_lim = [-4, 4];
y_axis_lim = [-4, 4];

marg_plot_x_num = 800;
marg_plot_y_num = 800;

marg_plot_x_min = -4;
marg_plot_x_max = 4;
marg_plot_y_min = -4;
marg_plot_y_max = 4;

marg_plot_edge_x = linspace(marg_plot_x_min, marg_plot_x_max, marg_plot_x_num + 1)';
marg_plot_edge_y = linspace(marg_plot_y_min, marg_plot_y_max, marg_plot_y_num + 1)';
marg_plot_x = (marg_plot_edge_x(1:end - 1) + marg_plot_edge_x(2:end)) / 2;
marg_plot_y = (marg_plot_edge_y(1:end - 1) + marg_plot_edge_y(2:end)) / 2;
[marg_plot_grid_x, marg_plot_grid_y] = meshgrid(marg_plot_x, marg_plot_y);
marg_plot_grid = [marg_plot_grid_x(:), marg_plot_grid_y(:)];

Meas_ref = ProbMeas2DMixNorm(meas_ref_vertices, meas_ref_triangles, meas_ref_mixnorm);

meas_ref_dens = reshape(Meas_ref.densityFunction(marg_plot_grid), marg_plot_x_num, marg_plot_y_num);

% instantiate the marginals
marg_cell = cell(marg_num, 1);
marg_dens_cell = cell(marg_num, 1);

for marg_id = 1:marg_num
    Meas = ProbMeas2DMixNorm(marg_vertices_cell{marg_id}, ...
        marg_triangles_cell{marg_id}, ...
        marg_mixnorm_cell{marg_id});

    marg_dens_cell{marg_id} = reshape(Meas.densityFunction(marg_plot_grid), marg_plot_x_num, marg_plot_y_num);

    marg_cell{marg_id} = Meas;
end

dens_max = max(max(max(meas_ref_dens)), max(max(vertcat(marg_dens_cell{:}))));

figure('Position', [100, 100, 1450, 265]);
ha = tight_subplot(1, 6, [0, 0.003], [0.101, 0.003], [0.0005, 0.027]);

axes(ha(1));
hold on;
plot_color = pcolor(marg_plot_grid_x, marg_plot_grid_y, meas_ref_dens);
plot_color.EdgeColor = 'interp';
plot_color.FaceColor = 'interp';

if plot_testfuncs
    Meas_ref.setSimplicialTestFuncs(meas_ref_testfuncs_cell{plot_test_id}{:});
    triplot(Meas_ref.SimplicialTestFuncs.Triangles, ...
        Meas_ref.SimplicialTestFuncs.Vertices(:, 1), Meas_ref.SimplicialTestFuncs.Vertices(:, 2));
end

clim([0, dens_max]);
set(gca, 'Color', 'black');
set(gca, 'XLim', x_axis_lim);
set(gca, 'YLim', y_axis_lim);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);

box on;
colormap('hot');

xlabel('$\mu_{0}$', 'Interpreter', 'latex', 'FontSize', 20);

for marg_id = 1:marg_num
    axes(ha(marg_id + 1));
    hold on;
    plot_color = pcolor(marg_plot_grid_x, marg_plot_grid_y, marg_dens_cell{marg_id});
    plot_color.EdgeColor = 'interp';
    plot_color.FaceColor = 'interp';

    if plot_testfuncs
        Meas = marg_cell{marg_id};
        Meas.setSimplicialTestFuncs(marg_testfuncs_cell{plot_test_id}{marg_id}{:});
        triplot(Meas.SimplicialTestFuncs.Triangles, ...
            Meas.SimplicialTestFuncs.Vertices(:, 1), Meas.SimplicialTestFuncs.Vertices(:, 2));
    end

    box on;
    colormap('hot');

    clim([0, dens_max]);
    set(gca, 'Color', 'black');
    set(gca, 'XLim', x_axis_lim);
    set(gca, 'YLim', y_axis_lim);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);

    xlabel(sprintf('$\\mu_{%d}$', marg_id), 'Interpreter', 'latex', 'FontSize', 20);
end

cb = colorbar('manual');
cb.Position = [0.9765, 0.101, 0.008, 0.897];
