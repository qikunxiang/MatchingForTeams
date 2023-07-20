rst_file_path = 'experiments/experiment1/exp1_MMOT_rst.mat';
input_file_path = 'experiments/experiment1/exp1_inputs.mat';

load(input_file_path);

if ~exist('MCsamps', 'var')
    load(rst_file_path, 'MCsamp_cell');
    test_num = length(testfuncs_cell);
    plotted_test = test_num;
    MCsamps = MCsamp_cell{plotted_test};
    clear MCsamp_cell;
end

disc_samp_indices = 21:30;
cont_samp_indices = 31:40;

color_cell = ...
    {[0.8627, 0.0784, 0.2353]; ...
    [1.0000, 0.0784, 0.5765]; ...
    [1.0000, 0.2706, 0]; ...
    [0.9686, 0.0980, 0.0039]; ...
    [0.4863, 0.9882, 0]; ...
    [0, 0.9804, 0.6039]; ...
    [0, 0.5451, 0.5451]; ...
    [0, 0, 0.8039]; ...
    [0.2941, 0, 0.5098]; ...
    [0.8549, 0.4392, 0.8392]};

color_map_length = 256;

color_map_cell = cell(marg_num, 1);

for marg_id = 1:marg_num
    interp_list = linspace(0, 1, color_map_length - 1)';
    color_map_cell{marg_id} = [1, 1, 1; ...
        (1 - interp_list) .* [1, 1, 1] + interp_list ...
        .* color_cell{marg_id}];
end

color_map = vertcat(color_map_cell{:});

% figure showing the continuous-to-discrete coupling

figure('Position', [100, 100, 500, 500]);
ha = tight_subplot(1, 1, [0, 0], [0.04, 0.05], [0.04, 0.01]);
axes(ha);
hold on;
box on;

colormap(color_map);

for marg_id = 1:marg_num
    marg = ProbMeas2D_CPWADens(marg_vert_cell{marg_id}, ...
        marg_tri_cell{marg_id}, marg_dens_cell{marg_id});
    marg_plot_x = plot_pt_x_cell{marg_id};
    marg_plot_y = plot_pt_y_cell{marg_id};
    marg_density = reshape(marg.Dens.Func(plot_pts_cell{marg_id}), ...
        length(marg_plot_x), length(marg_plot_y));
    marg_density_color = floor(marg_density ...
        / max(max(marg_density)) * color_map_length) + 1 ...
        + (marg_id - 1) * color_map_length;
    marg_inside = reshape(plot_inside_cell{ marg_id}, ...
        length(marg_plot_x), length(marg_plot_y));

    [grid_x, grid_y] = meshgrid(marg_plot_x, marg_plot_y);

    pc_plot = pcolor(grid_x, grid_y, marg_density_color);
    pc_plot.EdgeColor = 'none';
    pc_plot.AlphaData = marg_inside;
    pc_plot.FaceColor = 'texturemap';
    pc_plot.FaceAlpha = 'texturemap';
end

for poly_id = 1:length(quality_poly_cell)
    q = polyshape(quality_poly_cell{poly_id}(:, 1), ...
        quality_poly_cell{poly_id}(:, 2));
    plot(q, 'FaceColor', 'black');
end

for samp_id = disc_samp_indices
    qpt = MCsamps.DiscreteQualities(samp_id, :);

    for marg_id = 1:marg_num
        pt = MCsamps.ContinuousPoints{marg_id}(samp_id, :);
        scatter(pt(1), pt(2), 10, 'black', 'filled');

        line([qpt(1), pt(1)], [qpt(2), pt(2)], ...
            'Color', color_cell{marg_id});
    end
end

for samp_id = disc_samp_indices
    qpt = MCsamps.DiscreteQualities(samp_id, :);
    scatter(qpt(1), qpt(2), 15, 'black', 'filled');
end

set(gca, 'XLim', [-18, 18]);
set(gca, 'YLim', [-18, 22]);

title('$\hat{\gamma}_1,\;\ldots\;,\hat{\gamma}_{10}$', ...
    'Interpreter', 'latex', 'FontSize', 20);


% figure showing the continuous-to-continuous coupling

figure('Position', [600, 100, 500, 500]);
ha = tight_subplot(1, 1, [0, 0], [0.04, 0.05], [0.04, 0.01]);
axes(ha);
hold on;
box on;

colormap(color_map);

for marg_id = 1:marg_num
    marg = ProbMeas2D_CPWADens(marg_vert_cell{marg_id}, ...
        marg_tri_cell{marg_id}, marg_dens_cell{marg_id});
    marg_plot_x = plot_pt_x_cell{marg_id};
    marg_plot_y = plot_pt_y_cell{marg_id};
    marg_density = reshape(marg.Dens.Func(plot_pts_cell{marg_id}), ...
        length(marg_plot_x), length(marg_plot_y));
    marg_density_color = floor(marg_density ...
        / max(max(marg_density)) * color_map_length) + 1 ...
        + (marg_id - 1) * color_map_length;
    marg_inside = reshape(plot_inside_cell{ marg_id}, ...
        length(marg_plot_x), length(marg_plot_y));

    [grid_x, grid_y] = meshgrid(marg_plot_x, marg_plot_y);

    pc_plot = pcolor(grid_x, grid_y, marg_density_color);
    pc_plot.EdgeColor = 'none';
    pc_plot.AlphaData = marg_inside;
    pc_plot.FaceColor = 'texturemap';
    pc_plot.FaceAlpha = 'texturemap';
end

for poly_id = 1:length(quality_poly_cell)
    q = polyshape(quality_poly_cell{poly_id}(:, 1), ...
        quality_poly_cell{poly_id}(:, 2));
    plot(q, 'FaceColor', 'black');
end

for samp_id = cont_samp_indices
    qpt = MCsamps.DiscreteQualities(samp_id, :);

    for marg_id = 1:marg_num
        pt = MCsamps.ContinuousPoints{marg_id}(samp_id, :);
        scatter(pt(1), pt(2), 10, 'black', 'filled');

        line([qpt(1), pt(1)], [qpt(2), pt(2)], ...
            'Color', color_cell{marg_id});
    end
end

for samp_id = cont_samp_indices
    qpt = MCsamps.DiscreteQualities(samp_id, :);
    scatter(qpt(1), qpt(2), 15, 'black', 'filled');
end

set(gca, 'XLim', [-18, 18]);
set(gca, 'YLim', [-18, 22]);

title('$\tilde{\gamma}_1,\;\ldots\;,\tilde{\gamma}_{10}$', ...
    'Interpreter', 'latex', 'FontSize', 20);