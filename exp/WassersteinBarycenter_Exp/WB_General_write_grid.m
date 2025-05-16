% Generate csv files containing grid points for the barycenter will subsequently be used to test other methods for computing
% Wasserstein barycenter

CONFIG = WB_General_config();

load(CONFIG.SAVEPATH_INPUTS);

grid = quality_testfuncs_cell{end}{1};

writematrix([size(grid); grid], CONFIG.SAVEPATH_INPUTS_GRID);