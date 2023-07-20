rst_file_path = 'experiments/experiment1/exp1_MMOT_rst.mat';
input_file_path = 'experiments/experiment1/exp1_inputs.mat';

load(input_file_path, ...
    'plot_pt_x_cell', 'plot_pt_y_cell', 'plot_pts_cell', ...
    'plot_inside_cell', 'plot_quality_pt_x', 'plot_quality_pt_y', ...
    'plot_quality_pts', 'plot_quality_inside', 'quality_poly_cell');

if ~exist('parfunc', 'var')
    load(rst_file_path, 'parfunc_cell', 'transfunc_cell');

    [test_num, marg_num] = size(parfunc_cell);
    plotted_step = test_num;
    
    parfunc = parfunc_cell(plotted_step, :);
    transfunc = transfunc_cell{plotted_step};

    clear parfunc_cell;
    clear transfunc_cell;
end

quality_poly_num = length(quality_poly_cell);

color_map = jet;

% figure containing the transfer functions
figure('Position', [100, 100, 800, 300]);
ha = tight_subplot(2, 5, [0.09, 0.01], [0.06, 0.02], [0.01, 0.10]);

min_transfunc = min(min(transfunc));
max_transfunc = max(max(transfunc));

for marg_id = 1:marg_num
    axes(ha(marg_id)); %#ok<LAXES> 
    hold on;
    box on;
    colormap(color_map);

    grid_x_num = length(plot_quality_pt_x);
    grid_y_num = length(plot_quality_pt_y);
    transfunc_val = reshape(transfunc(:, marg_id), ...
        grid_x_num, grid_y_num);
    transfunc_inside = reshape(plot_quality_inside, grid_x_num, ...
        grid_y_num);

    [grid_x, grid_y] = meshgrid(plot_quality_pt_x, plot_quality_pt_y);

    pc_plot = pcolor(grid_x, grid_y, transfunc_val);
    pc_plot.EdgeColor = 'none';
    pc_plot.AlphaData = transfunc_inside;
    pc_plot.FaceColor = 'texturemap';
    pc_plot.FaceAlpha = 'texturemap';

    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'XLim', [min(plot_quality_pt_x) - 0.2, ...
        max(plot_quality_pt_x) + 0.2]);
    set(gca, 'YLim', [min(plot_quality_pt_y) - 0.2, ...
        max(plot_quality_pt_y) + 0.2]);
    xlabel(sprintf('$\\tilde{\\varphi}_{%d}$', marg_id), ...
        'Interpreter', 'latex');
end

cb = colorbar;
cb.Position = [0.92, 0.03, 0.035, 0.94];
