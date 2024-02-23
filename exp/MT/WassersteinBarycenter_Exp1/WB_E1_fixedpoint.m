% Compute the 2-Wasserstein barycenter using the fixed-point scheme

CONFIG = WB_E1_config();

load(CONFIG.SAVEPATH_INPUTS);

mean_vec_cell = cell(marg_num, 1);
cov_mat_cell = cell(marg_num, 1);

for marg_id = 1:marg_num
    Meas = ProbMeas2D_MixNorm(marg_vertices_cell{marg_id}, ...
        marg_triangles_cell{marg_id}, ...
        marg_mixnorm_cell{marg_id});

    mean_vec_cell{marg_id} = Meas.meanVector();
    cov_mat_cell{marg_id} = Meas.covarianceMatrix();
end

[WB_fp_mean_vec, WB_fp_cov_mat, WB_fp_min_cost, WB_fp_cost_iter] ...
    = WB_locationscatter(mean_vec_cell, cov_mat_cell, marg_weights, 1e-10);

fprintf('Minimum cost computed by the fixed-point scheme: %.4f\n', ...
    WB_fp_min_cost);

save(CONFIG.SAVEPATH_FIXEDPOINT, ...
    'WB_fp_mean_vec', ...
    'WB_fp_cov_mat', ...
    'WB_fp_min_cost', ...
    'WB_fp_cost_iter', ...
    '-v7.3');