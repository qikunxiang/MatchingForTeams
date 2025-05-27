% Plot the computed business location distributions

CONFIG = BizLoc_config();

load(CONFIG.SAVEPATH_INPUTS);
load(CONFIG.SAVEPATH_OT);
load(CONFIG.SAVEPATH_OUTPUTS_UB);

disc_atoms = cell(test_num, 1);
disc_probs = cell(test_num, 1);

cont_atoms = cell(test_num, 1);
cont_probs = cell(test_num, 1);
cont_densi = cell(test_num, 1);

% if a single cell in the histogram has probability at least 0.1, it is regarded as an atom at the center of the cell
cont_atom_thres = 5e-3;

cell_area = (quality_hist_edge_x(2) - quality_hist_edge_x(1)) * (quality_hist_edge_y(2) - quality_hist_edge_y(1));

for test_id = 1:test_num
    cand = OT_info_cell{test_id}.discrete.candidate;
    disc_atoms{test_id} = cand.Atoms;
    disc_probs{test_id} = cand.Probabilities;

    cont_densi{test_id} = MT_histpdf_cell{test_id}';
    cell_probs = cont_densi{test_id} * cell_area;

    [a_r, a_c, a_v] = find(cell_probs >= cont_atom_thres);
    a_lin = sub2ind([quality_hist_x_num, quality_hist_y_num], a_r, a_c);
    cont_atoms{test_id} = [quality_plot_hist_grid_x(a_lin), quality_plot_hist_grid_y(a_lin)];
    cont_probs{test_id} = cell_probs(a_lin);
    cont_densi{test_id}(a_lin) = 0;
end

probs_agg = [vertcat(disc_probs{:}); vertcat(cont_probs{:})];
prob_max = max(probs_agg);
prob_min = min(probs_agg);
dens_max = max(max(vertcat(cont_densi{:})));

figure('Position', [100, 100, 1350, 550]);
ha = tight_subplot(2, test_num, [0.035, 0.012], [0.065, 0.004], [0.012, 0.004]);
theo_sparsity = inf(test_num, 1);

for test_id = 1:test_num
    axes(ha(test_id));
    hold on;
    box on;
    bbc = bubblechart(disc_atoms{test_id}(:, 1), disc_atoms{test_id}(:, 2), disc_probs{test_id}, ...
        'red', 'MarkerEdgeAlpha', 0.6, 'MarkerFaceAlpha', 0.3);
    bubblelim([prob_min, prob_max * 2.7]);

    set(gca, 'XLim', [quality_x_min - 0.2, quality_x_max + 0.2]);
    set(gca, 'YLim', [quality_y_min - 0.2, quality_y_max + 0.2]);
    set(gca, 'XTick', -2:1:2);
    set(gca, 'YTick', -3:1:2);

    axes(ha(test_id + test_num));
    hold on;
    box on;

    pc = pcolor(quality_plot_hist_grid_x, quality_plot_hist_grid_y, cont_densi{test_id});
    pc.FaceColor = 'interp';
    pc.EdgeColor = 'interp';

    clim([0, dens_max * 0.3]);
    colormap(1 - (gray()) .^ 0.7);
    

    bubblechart(cont_atoms{test_id}(:, 1), cont_atoms{test_id}(:, 2), cont_probs{test_id}, ...
        'red', 'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.6);
    bubblelim([prob_min, prob_max * 2.7]);

    for atom_id = 1:length(cont_probs{test_id})
        text(cont_atoms{test_id}(atom_id, 1) - 0.15, cont_atoms{test_id}(atom_id, 2), ...
            sprintf('%.3f', cont_probs{test_id}(atom_id)), 'Color', 'red', 'FontSize', 7);
    end


    set(gca, 'XLim', [quality_x_min - 0.2, quality_x_max + 0.2]);
    set(gca, 'YLim', [quality_y_min - 0.2, quality_y_max + 0.2]);
    set(gca, 'XTick', -2:1:2);
    set(gca, 'YTick', -3:1:2);

    marg_testfuncs_num_sum = 0;

    for marg_id = 1:marg_num
        tf_num = size(marg_testfuncs_cell{test_id}{marg_id}{1}, 1) - 1;
        marg_testfuncs_num_sum = marg_testfuncs_num_sum + tf_num;
        theo_sparsity(test_id) = min(theo_sparsity(test_id), tf_num);
    end

    quality_testfuncs_num = size(quality_testfuncs_cell{test_id}{1}, 1) - 1;
    theo_sparsity(test_id) = theo_sparsity(test_id) + quality_testfuncs_num + 2;
    xlabel(sprintf('$n=%d$', marg_testfuncs_num_sum + marg_num * (quality_testfuncs_num + 1)), ...
        'FontSize', 18, 'Interpreter', 'latex');

    fprintf(['test %d: theoretical sparsity = %d, ' 'actual sparsity = %d\n'], ...
        test_id, theo_sparsity(test_id), size(disc_atoms{test_id}, 1));
end