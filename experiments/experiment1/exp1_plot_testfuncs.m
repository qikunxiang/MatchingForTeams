input_file_path = 'experiments/experiment1/exp1_inputs.mat';

load(input_file_path);

test_num = length(testfuncs_cell);

plotted_test = 4;

text_pos = ...
    [0, 15;
    0, 8; ...
    -12, 12; ...
    -10, -1; ...
    -10, -13; ...
    -2.7, -13; ...
    0, -9; ...
    10, -9; ...
    12, 2; ...
    12, 15];
quality_text_pos = [0, -0.5];

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

figure('Position', [100, 100, 500, 500]);
ha = tight_subplot(1, 1, [0, 0], [0.04, 0.01], [0.04, 0.01]);
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

    triplot(vertcat(testfuncs_cell{plotted_test}{marg_id}{2}{:}), ...
        testfuncs_cell{plotted_test}{marg_id}{1}(:, 1), ...
        testfuncs_cell{plotted_test}{marg_id}{1}(:, 2), ...
        'Color', [0.5, 0.5, 0.5]);
end

for poly_id = 1:length(quality_poly_cell)
    q = polyshape(quality_poly_cell{poly_id}(:, 1), ...
        quality_poly_cell{poly_id}(:, 2));
    plot(q, 'FaceColor', 'black');
end

triplot(quality_testfuncs_cell{plotted_test}{2}, ...
    quality_testfuncs_cell{plotted_test}{1}(:, 1), ...
    quality_testfuncs_cell{plotted_test}{1}(:, 2), ...
    'Color', [0.5, 0.5, 0.5]);

set(gca, 'XLim', [-18, 18]);
set(gca, 'YLim', [-18, 22]);

for marg_id = 1:marg_num
    text(text_pos(marg_id, 1) - 0.7, text_pos(marg_id, 2), ...
        sprintf('$\\mu_{%d}$', marg_id), ...
        'Interpreter', 'latex', 'FontSize', 25);
end

text(quality_text_pos(1) - 0.7, quality_text_pos(2), ...
    '$\mathcal{Z}$', ...
    'Interpreter', 'latex', 'FontSize', 25);