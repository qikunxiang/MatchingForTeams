input_file_path = 'experiments/experiment1/exp1_inputs.mat';
OT_file_path = 'experiments/experiment1/exp1_OT.mat';

load(input_file_path);
load(OT_file_path);

test_num = length(testfuncs_cell);

plotted_test = 4;

color_map = jet;

figure('Position', [100, 100, 1500, 600]);
ha = tight_subplot(2, 5, [0.04, 0.01], [0.04, 0.01], [0.01, 0.01]);

Laguerre_subsample = 1:20:10000;
edge_color = [0.4, 0.4, 0.4];
line_style = ':';
line_width = 1.25;

for marg_id = 1:marg_num

    axes(ha(marg_id)); %#ok<LAXES> 
    hold on;
    box on;
    
    colormap(color_map);

    marg = ProbMeas2D_CPWADens(marg_vert_cell{marg_id}, ...
        marg_tri_cell{marg_id}, marg_dens_cell{marg_id});
    marg.setSimplicialTestFuncs(testfuncs_cell{plotted_test}{marg_id}{:});
    marg_plot_x = plot_pt_x_cell{marg_id};
    marg_plot_y = plot_pt_y_cell{marg_id};

    atoms = marg.SimplicialTestFuncs.Vertices;
    probs = marg.SimplicialTestFuncs.Integrals;

    marg.setOTWeights(atoms, probs, ...
        OT_cell{plotted_test}{marg_id}{:});

    color_indices = ceil(probs / max(probs) * size(color_map, 1));

    fill(marg.OT.Laguerre.x(Laguerre_subsample, :), ...
        marg.OT.Laguerre.y(Laguerre_subsample, :), color_indices, ...
        'FaceAlpha', 0.6, 'EdgeColor', edge_color, ...
        'LineStyle', line_style, 'LineWidth', line_width);

    if marg_id == 2
        fill([6; 2; 6], [12; 8; 4], 'white', ...
            'EdgeColor', edge_color, ...
            'LineStyle', line_style, 'LineWidth', line_width);
        fill([-6; -2; -6], [12; 8; 4], 'white', ...
            'EdgeColor', edge_color, ...
            'LineStyle', line_style, 'LineWidth', line_width);
        fill([5.95; 6.05; 6.05; 5.95], [11.95; 11.95; 4.05; 4.05], ...
            'white', 'EdgeColor', 'none');
        fill([-5.95; -6.05; -6.05; -5.95], [11.95; 11.95; 4.05; 4.05], ...
            'white', 'EdgeColor', 'none');
    elseif marg_id == 6
        fill([4; 0; -4], [-11; -13; -11], 'white', ...
            'EdgeColor', edge_color, ...
            'LineStyle', line_style, 'LineWidth', line_width);
        fill([4; 0; -4], [-15; -13; -15], 'white', ...
            'EdgeColor', edge_color, ...
            'LineStyle', line_style, 'LineWidth', line_width);
        fill([3.95; 3.95; -3.95; -3.95], ...
            [-10.95; -11.05; -11.05; -10.95], ...
            'white', 'EdgeColor', 'none');
        fill([3.95; 3.95; -3.95; -3.95], ...
            [-14.95; -15.05; -15.05; -14.95], ...
            'white', 'EdgeColor', 'none');
    elseif marg_id == 8
        fill([14; 10; 14], [-8; -11; -14], 'white', ...
            'EdgeColor', edge_color, ...
            'LineStyle', line_style, 'LineWidth', line_width);
        fill([6; 10; 6], [-8; -11; -14], 'white', ...
            'EdgeColor', edge_color, ...
            'LineStyle', line_style, 'LineWidth', line_width);
        fill([13.95; 14.05; 14.05; 13.95], ...
            [-8.05; -8.05; -13.95; -13.95], ...
            'white', 'EdgeColor', 'none');
        fill([5.95; 6.05; 6.05; 5.95], ...
            [-8.05; -8.05; -13.95; -13.95], ...
            'white', 'EdgeColor', 'none');
    end

    for atom_id = 1:size(atoms, 1)
        scatter(atoms(atom_id, 1), atoms(atom_id, 2), 10, 'black', ...
            'filled', 'o', 'MarkerFaceColor', ...
            color_map(color_indices(atom_id, :), :), ...
            'MarkerEdgeColor', 'black');
    end


    leftbound_x = min(plot_pt_x_cell{marg_id});
    rightbound_x = max(plot_pt_x_cell{marg_id});
    bottombound_y = min(plot_pt_y_cell{marg_id});
    topbound_y = max(plot_pt_y_cell{marg_id});
    square_length = max(rightbound_x - leftbound_x, ...
        topbound_y - bottombound_y) + 0.6;
    pad_x = (square_length - (rightbound_x - leftbound_x)) / 2;
    pad_y = (square_length - (topbound_y - bottombound_y)) / 2;
    
    set(gca, 'XLim', [leftbound_x - pad_x, rightbound_x + pad_x]);
    set(gca, 'YLim', [bottombound_y - pad_y, topbound_y + pad_y]);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);

    xlabel(sprintf('$\\mu_{%d}$', marg_id), 'Interpreter', 'latex');
end
