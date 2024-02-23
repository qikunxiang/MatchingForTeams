% Plot the density functions of the reference measure and the marginals
% from the same location-scatter family

CONFIG = WB_E1_config();

load(CONFIG.SAVEPATH_INPUTS);

plot_testfuncs = false;
plot_test_id = 3;

x_axis_lim = [-6.5, 6.5];
y_axis_lim = [-6.5, 6.5];

Meas_ref = ProbMeas2D_MixNorm(meas_ref_vertices, ...
    meas_ref_triangles, meas_ref_mixnorm);

meas_ref_dens = reshape( ...
    Meas_ref.densityFunction(quality_plot_hist_grid), ...
    quality_hist_x_num, quality_hist_y_num);

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

dens_max = max(max(max(meas_ref_dens)), max(max(vertcat( ...
    marg_dens_cell{:}))));

figure('Position', [100, 100, 1620, 265]);
ha = tight_subplot(1, 6, [0, 0.012], [0.125, 0.025], [0.009, 0.04]);

axes(ha(1));
hold on;
plot_color = pcolor(quality_plot_hist_grid_x, ...
    quality_plot_hist_grid_y, ...
    meas_ref_dens);
plot_color.EdgeColor = 'interp';
plot_color.FaceColor = 'interp';

if plot_testfuncs
    Meas_ref.setSimplicialTestFuncs( ...
        meas_ref_testfuncs_vertices_cell{plot_test_id}, ...
        meas_ref_testfuncs_triangles_cell{plot_test_id}, false);
    triplot(Meas_ref.SimplicialTestFuncs.Triangles, ...
        Meas_ref.SimplicialTestFuncs.Vertices(:, 1), ...
        Meas_ref.SimplicialTestFuncs.Vertices(:, 2));
end

clim([0, dens_max]);
set(gca, 'Color', 'black');
set(gca, 'XLim', x_axis_lim);
set(gca, 'YLim', y_axis_lim);

box on;
colormap('hot');

xlabel('$\mu_{0}$', 'Interpreter', 'latex', 'FontSize', 14);

for marg_id = 1:marg_num
    axes(ha(marg_id + 1));
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
cb.Position = [0.965, 0.125, 0.015, 0.85];
