% Evaluate the performance of the multi-marginal optimal transport (MMOT) method of Neufeld, Xiang via semi-discrete Wasserstein-2 OT

CONFIG = WB_General_config();

load(CONFIG.SAVEPATH_INPUTS);

optimization_options = struct;
optimization_options.OptimalityTolerance = 1e-5;
optimization_options.Display = 'iter-detailed';

marg_cell = cell(marg_num, 1);

for marg_id = 1:marg_num
    marg_cell{marg_id} = ProbMeas2DCPWADens( ...
        marg_vertices_cell{marg_id}, ...
        marg_triangles_cell{marg_id}, ...
        marg_densities_cell{marg_id});
end


MT = MatchTeam2DWassersteinBarycenter( ...
    marg_cell, ...
    marg_weights, ...
    {quality_vertices}, ...
    [], []);

NX_file = load(CONFIG.LOADPATH_RESULTS_MMOTNEUFELDXIANG);
WB_atoms = NX_file.WB_atoms;
WB_probs = NX_file.WB_probs;

atoms_positive = WB_probs >= 1e-6;
WB_atoms = WB_atoms(atoms_positive, :);
WB_probs = WB_probs(atoms_positive);
WB_probs = WB_probs / sum(WB_probs);

[objective, W2_dist_list] = MT.evaluateObjectiveViaW2( ...
    WB_atoms, WB_probs, optimization_options);

fprintf('objective = %10.4f\n', objective);

save(CONFIG.SAVEPATH_EVALUATION_MMOTNEUFELDXIANG, ...
    'objective', ...
    'W2_dist_list', ...
    '-v7.3');