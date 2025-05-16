% Generate csv files containing samples from the input measures that will subsequently be used to test other methods for computing
% Wasserstein barycenter

CONFIG = WB_General_config();

load(CONFIG.SAVEPATH_INPUTS);

% generate 60000 samples per input measure
samp_num = 60000;

% generate enough samples for 500 epochs
epoch_num = 500;

RS = RandStream('mrg32k3a', 'Seed', 3500);

for marg_id = 1:marg_num
    Meas = ProbMeas2DCPWADens( ...
        marg_vertices_cell{marg_id}, ...
        marg_triangles_cell{marg_id}, ...
        marg_densities_cell{marg_id});

    RS.Substream = marg_id;

    for epoch_id = 1:epoch_num
        samps = Meas.randSample(samp_num, RS);
        writematrix(samps, sprintf(CONFIG.SAVEFORMAT_INPUTS_SAMPLES, marg_id, epoch_id));
    end
end