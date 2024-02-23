% Plot the histogram of computed approximate Wasserstein barycenter

CONFIG = WB_E1_config();

load(CONFIG.SAVEPATH_INPUTS);
load(CONFIG.SAVEPATH_OUTPUTS);
load(CONFIG.SAVEPATH_FIXEDPOINT);

% retrieve the 2-Wasserstein barycenter via transforming the reference
% measure

trans_mat = chol(WB_fp_cov_mat)';
shift_vec = WB_fp_mean_vec;
meas_WB_fp_vertices = meas_ref_vertices * trans_mat' + shift_vec';
meas_WB_fp_triangles = meas_ref_triangles;

comp_num = length(meas_ref_mixnorm.weights);
meas_WB_fp_mixnorm = struct;
meas_WB_fp_mixnorm.weights = meas_ref_mixnorm.weights;
meas_WB_fp_mixnorm.components = cell(comp_num, 1);

for comp_id = 1:comp_num
    meas_WB_fp_mixnorm.components{comp_id} = struct;
    meas_WB_fp_mixnorm.components{comp_id}.mean_vec = trans_mat ...
        * meas_ref_mixnorm.components{comp_id}.mean_vec + shift_vec;
    meas_WB_fp_mixnorm.components{comp_id}.cov_mat = trans_mat ...
        * meas_ref_mixnorm.components{comp_id}.cov_mat * trans_mat';
end

x_axis_lim = [-6.5, 6.5];
y_axis_lim = [-6.5, 6.5];

x_label_cell = cell(test_num, 1);

for test_id = 1:test_num
    marg_testfuncs_num_sum = 0;

    for marg_id = 1:marg_num
        marg_testfuncs_num_sum = marg_testfuncs_num_sum ...
            + size(marg_testfuncs_cell{test_id}{marg_id}{1}, 1) - 1;
    end

    quality_testfuncs_num = size(quality_testfuncs_cell{test_id}{1}, 1) ...
        - 1;

    x_label_cell{test_id} = sprintf('(%d, %d)', ...
        marg_testfuncs_num_sum, quality_testfuncs_num);
end

Meas_WB = ProbMeas2D_MixNorm(meas_WB_fp_vertices, ...
    meas_WB_fp_triangles, meas_WB_fp_mixnorm);

meas_WB_dens = reshape( ...
    Meas_WB.densityFunction(quality_plot_hist_grid), ...
    quality_hist_x_num, quality_hist_y_num);

dens_max = max(max(max(meas_WB_dens)), max(max(vertcat( ...
    WB_histpdf_cell{:}))));

figure('Position', [0, 100, 1830, 265]);
ha = tight_subplot(1, test_num + 1, [0, 0.01], ...
    [0.125, 0.025], [0.008, 0.005]);


for test_id = 1:test_num
    axes(ha(test_id));
    hold on;

    plot_color = pcolor(quality_plot_hist_grid_x, ...
        quality_plot_hist_grid_y, ...
        WB_histpdf_cell{test_id}');
    plot_color.EdgeColor = 'interp';
    plot_color.FaceColor = 'interp';

    box on;
    colormap('hot');

    clim([0, dens_max * 0.6]);
    set(gca, 'Color', 'black');
    set(gca, 'XLim', x_axis_lim);
    set(gca, 'YLim', y_axis_lim);

    xlabel(x_label_cell{test_id}, 'FontWeight', 'bold');
end


axes(ha(end));
hold on;
plot_color = pcolor(quality_plot_hist_grid_x, ...
    quality_plot_hist_grid_y, ...
    meas_WB_dens);
plot_color.EdgeColor = 'interp';
plot_color.FaceColor = 'interp';

clim([0, dens_max * 0.6]);
set(gca, 'Color', 'black');
set(gca, 'XLim', x_axis_lim);
set(gca, 'YLim', y_axis_lim);
xlabel('fixed-point', 'FontWeight', 'bold');

box on;
colormap('hot');