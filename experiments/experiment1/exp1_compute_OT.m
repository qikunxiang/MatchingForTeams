load_file_path = 'experiments/experiment1/exp1_inputs.mat';
save_path_path = 'experiments/experiment1/exp1_OT.mat';

rng(1000, 'combRecursive');

load(load_file_path);

options = struct;
options.log_file = 'exp1_OT.log';
options.OT = struct;
options.OT.optimization_options = struct;
options.OT.optimization_options.Display = 'iter-detailed';

exp_log_file_path = 'exp1_compute_OT.log';

test_num = length(testfuncs_cell);

OT_cell = cell(test_num, 1);

exp_log_file = fopen(exp_log_file_path, 'a');

if exp_log_file < 0
    error('cannot open log file');
end

fprintf(exp_log_file, '--- experiment starts ---\n');

for test_id = 1:test_num
    marg_cell = cell(marg_num, 1);

    for marg_id = 1:marg_num
        marg_cell{marg_id} = ProbMeas2D_CPWADens( ...
            marg_vert_cell{marg_id}, ...
            marg_tri_cell{marg_id}, ...
            marg_dens_cell{marg_id});
    end

    MT = MT2DQuad_MMOT(marg_cell, marg_weights, quality_poly_cell, ...
            options);
    MT.setSimplicialTestFuncs(testfuncs_cell{test_id});

    OT_cell{test_id} = MT.performReassembly();

    fprintf(exp_log_file, ...
        'test %d: semi-discrete OT completed\n', test_id);

    save(save_path_path, 'OT_cell');
end

fprintf(exp_log_file, '--- experiment ends ---\n\n');
fclose(exp_log_file);