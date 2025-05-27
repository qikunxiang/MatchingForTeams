function CONFIG = WB_General_config()
% Place to store global configurations of the Wasserstein barycenter experiment with measures with continuous piece-wise affine
% densities
% Output:
%   CONFIG: a struct containing configurations as fields

CONFIG = global_config();

% path to save files
CONFIG.SAVEPATH = [CONFIG.SAVEPATH_ROOT, 'MT/WassersteinBarycenter_Exp/'];

% path to save files (local)
CONFIG.SAVEPATH_LOCAL = '../LocalSaveFiles/MT/WassersteinBarycenter_Exp/';

% path to log files
CONFIG.LOGPATH = [CONFIG.LOGPATH_ROOT, 'MT/WassersteinBarycenter_Exp/'];

% root folder for gurobi log files
CONFIG.LOGPATH_GUROBI = [CONFIG.LOGPATH, 'Gurobi/'];

% if the directory does not exist, create it first
if ~exist(CONFIG.SAVEPATH, 'dir')
    mkdir(CONFIG.SAVEPATH);
end

% if the directory does not exist, create it first
if ~exist(CONFIG.SAVEPATH_LOCAL, 'dir')
    mkdir(CONFIG.SAVEPATH_LOCAL);
end

% if the directory does not exist, create it first
if ~exist(CONFIG.LOGPATH, 'dir')
    mkdir(CONFIG.LOGPATH);
end

% if the directory does not exist, create it first
if ~exist(CONFIG.LOGPATH_GUROBI, 'dir')
    mkdir(CONFIG.LOGPATH_GUROBI);
end

% file name of the inputs
CONFIG.FILENAME_INPUTS = 'inputs';
CONFIG.SAVEPATH_INPUTS = [CONFIG.SAVEPATH, CONFIG.FILENAME_INPUTS, '.mat'];

% file name of the optimal transport outputs
CONFIG.FILENAME_OUTPUTS_OT = 'outputs_OT';
CONFIG.SAVEPATH_OUTPUTS_OT = [CONFIG.SAVEPATH, CONFIG.FILENAME_OUTPUTS_OT, '.mat'];

% file name of the Wasserstein-2 optimal transport outputs
CONFIG.FILENAME_OUTPUTS_W2OT = 'outputs_W2OT';
CONFIG.SAVEPATH_OUTPUTS_W2OT = [CONFIG.SAVEPATH, CONFIG.FILENAME_OUTPUTS_W2OT, '.mat'];

% file name of the outputs for lower bounds
CONFIG.FILENAME_OUTPUTS_LB = 'outputs_LB';
CONFIG.SAVEPATH_OUTPUTS_LB = [CONFIG.SAVEPATH, CONFIG.FILENAME_OUTPUTS_LB, '.mat'];

% file name of the outputs for upper bounds
CONFIG.FILENAME_OUTPUTS_UB = 'outputs_UB';
CONFIG.SAVEPATH_OUTPUTS_UB = [CONFIG.SAVEPATH, CONFIG.FILENAME_OUTPUTS_UB, '.mat'];

% file name of the outputs for upper bounds based on Wasserstein-2 optimal transport
CONFIG.FILENAME_OUTPUTS_W2OTUB = 'outputs_W2OTUB';
CONFIG.SAVEPATH_OUTPUTS_W2OTUB = [CONFIG.SAVEPATH, CONFIG.FILENAME_OUTPUTS_W2OTUB, '.mat'];

% file name of the grid points for the barycenter
CONFIG.FILENAME_INPUTS_GRID = 'grid_points';
CONFIG.SAVEPATH_INPUTS_GRID = [CONFIG.SAVEPATH_LOCAL, CONFIG.FILENAME_INPUTS_GRID, '.csv'];

% file format string for the path to store samples from the input measures (used for comparison with other methods)
CONFIG.FILEFORMAT_INPUTS_SAMPLES = 'input_measure_%d_samples_epoch_%d';
CONFIG.SAVEFORMAT_INPUTS_SAMPLES = [CONFIG.SAVEPATH_LOCAL, CONFIG.FILEFORMAT_INPUTS_SAMPLES, '.csv'];

% file format string for the path to store samples from the output barycenter resulted from the Neural Wasserstein Barycenter (NWB)
% method of Fan, Taghvaei, Chen
CONFIG.FILEFORMAT_RESULTS_NWBFANTAGHVAEICHEN = 'results_NWBFanTaghvaeiChen_samples_epoch_%d';
CONFIG.PATHFORMAT_RESULTS_NWBFANTAGHVAEICHEN = [CONFIG.SAVEPATH_LOCAL, CONFIG.FILEFORMAT_RESULTS_NWBFANTAGHVAEICHEN, '.csv'];

% file name of the outputs for evaluating the Neural Wasserstein Barycenter (NWB) method of Fan, Taghvaei, Chen via semi-discrete
% Wasserstein-2 optimal transport
CONFIG.FILENAME_EVALUATION_NWBFANTAGHVAEICHEN = 'evaluation_NWBFanTaghvaeiChen';
CONFIG.SAVEPATH_EVALUATION_NWBFANTAGHVAEICHEN = ...
    [CONFIG.SAVEPATH, CONFIG.FILENAME_EVALUATION_NWBFANTAGHVAEICHEN, '.mat'];

% file format string for the path to store samples from the output barycenter resulted from the Wasserstein Iterative Network (WIN)
% method of Korotin, Egiazarian, Li, Burnaev
CONFIG.FILEFORMAT_RESULTS_WINKOROTINEGIAZARIANLIBURNAEV = 'results_WINKorotinEgiazarianLiBurnaev_samples_epoch_%d';
CONFIG.PATHFORMAT_RESULTS_WINKOROTINEGIAZARIANLIBURNAEV = ...
    [CONFIG.SAVEPATH_LOCAL, CONFIG.FILEFORMAT_RESULTS_WINKOROTINEGIAZARIANLIBURNAEV, '.csv'];

% file name of the outputs for evaluating the Wasserstein Iterative Network (WIN) method of Korotin, Egiazarian, Li, Burnaev via 
% semi-discrete Wasserstein-2 optimal transport
CONFIG.FILENAME_EVALUATION_WINKOROTINEGIAZARIANLIBURNAEV = 'evaluation_WINKorotinEgiazarianLiBurnaev';
CONFIG.SAVEPATH_EVALUATION_WINKOROTINEGIAZARIANLIBURNAEV = ...
    [CONFIG.SAVEPATH, CONFIG.FILENAME_EVALUATION_WINKOROTINEGIAZARIANLIBURNAEV, '.mat'];

% file format string for the path to store samples from the output barycenter resulted from the Continuous Wasserstein-2 Barycenter
% method without minimax optimization (CW2B) of Korotin, Li, Solomon, Burnaev
CONFIG.FILEFORMAT_RESULTS_CW2BKOROTINLISOLOMONBURNAEV = 'results_CW2BKorotinLiSolomonBurnaev_samples_epoch_%d';
CONFIG.PATHFORMAT_RESULTS_CW2BKOROTINLISOLOMONBURNAEV = ...
    [CONFIG.SAVEPATH_LOCAL, CONFIG.FILEFORMAT_RESULTS_CW2BKOROTINLISOLOMONBURNAEV, '.csv'];

% file name of the outputs for evaluating the Continuous Wasserstein-2 Barycenter method without minimax optimization (CW2B) of 
% Korotin, Li, Solomon, Burnaev via semi-discrete Wasserstein-2 optimal transport
CONFIG.FILENAME_EVALUATION_CW2BKOROTINLISOLOMONBURNAEV = 'evaluation_CW2BKorotinLiSolomonBurnaev';
CONFIG.SAVEPATH_EVALUATION_CW2BKOROTINLISOLOMONBURNAEV = ...
    [CONFIG.SAVEPATH, CONFIG.FILENAME_EVALUATION_CW2BKOROTINLISOLOMONBURNAEV, '.mat'];

% file format string for the path to store atom weights of the output barycenter resulted from the Parallel Streaming Wasserstein
% Barycenter (PSWB) method of Staib, Claici, Solomon, Jegelka
CONFIG.FILEFORMAT_RESULTS_PSWBSTAIBCLAICISOLOMONJEGELKA = 'results_PSWBStaibClaiciSolomonJegelka_barycenter_weights_epoch_%d';
CONFIG.PATHFORMAT_RESULTS_PSWBSTAIBCLAICISOLOMONJEGELKA = ...
    [CONFIG.SAVEPATH_LOCAL, CONFIG.FILEFORMAT_RESULTS_PSWBSTAIBCLAICISOLOMONJEGELKA, '.csv'];

% file name of the outputs for evaluating the Parallel Streaming Wasserstein Barycenter (PSWB) method of Staib, Claici, Solomon, 
% Jegelka via semi-discrete Wasserstein-2 optimal transport
CONFIG.FILENAME_EVALUATION_PSWBSTAIBCLAICISOLOMONJEGELKA = 'evaluation_PSWBStaibClaiciSolomonJegelka';
CONFIG.SAVEPATH_EVALUATION_PSWBSTAIBCLAICISOLOMONJEGELKA = ...
    [CONFIG.SAVEPATH, CONFIG.FILENAME_EVALUATION_PSWBSTAIBCLAICISOLOMONJEGELKA, '.mat'];

% file name of the output barycenter resulted from the multi-marginal optimal transport (MMOT) method of Neufeld, Xiang
CONFIG.FILENAME_RESULTS_MMOTNEUFELDXIANG = 'results_MMOTNeufeldXiang';
CONFIG.LOADPATH_RESULTS_MMOTNEUFELDXIANG = [CONFIG.SAVEPATH_LOCAL, CONFIG.FILENAME_RESULTS_MMOTNEUFELDXIANG, '.mat'];


% file name of the outputs for evaluating the multi-marginal optimal transport (MMOT) method of Neufeld, Xiang via semi-discrete 
% Wasserstein-2 optimal transport
CONFIG.FILENAME_EVALUATION_MMOTNEUFELDXIANG = 'evaluation_MMOTNeufeldXiang';
CONFIG.SAVEPATH_EVALUATION_MMOTNEUFELDXIANG = [CONFIG.SAVEPATH, CONFIG.FILENAME_EVALUATION_MMOTNEUFELDXIANG, '.mat'];

% file name of the outputs for evaluating our algorithm via semi-discrete Wasserstein-2 optimal transport
CONFIG.FILENAME_EVALUATION_OURALGO = 'evaluation_OurAlgo';
CONFIG.SAVEPATH_EVALUATION_OURALGO = [CONFIG.SAVEPATH, CONFIG.FILENAME_EVALUATION_OURALGO, '.mat'];

% file name of the optimal transport logs
CONFIG.LOGNAME_OT = 'OT';
CONFIG.LOGPATH_OT = [CONFIG.LOGPATH, CONFIG.LOGNAME_OT, '.log'];

% file name of the Wasserstein-2 optimal transport logs
CONFIG.LOGNAME_W2OT = 'W2OT';
CONFIG.LOGPATH_W2OT = [CONFIG.LOGPATH, CONFIG.LOGNAME_W2OT, '.log'];

% file name of the main logs
CONFIG.LOGNAME_MAIN = 'main';
CONFIG.LOGPATH_MAIN = [CONFIG.LOGPATH, CONFIG.LOGNAME_MAIN, '.log'];

% file name of the linear semi-infinite programming logs
CONFIG.LOGNAME_LSIP_MAIN = 'LSIP_main';
CONFIG.LOGPATH_LSIP_MAIN = [CONFIG.LOGPATH, CONFIG.LOGNAME_LSIP_MAIN, '.log'];

% file name of the linear programming logs
CONFIG.LOGNAME_LSIP_LP = 'LSIP_LP';
CONFIG.LOGPATH_LSIP_LP = [CONFIG.LOGPATH_GUROBI, CONFIG.LOGNAME_LSIP_LP, '.log'];

% file name of the global optimization logs
CONFIG.LOGNAME_LSIP_GLOBAL = 'LSIP_global';
CONFIG.LOGPATH_LSIP_GLOBAL = [CONFIG.LOGPATH, CONFIG.LOGNAME_LSIP_GLOBAL, '.log'];

% file name of the reassembly based upper bounds
CONFIG.LOGNAME_UB = 'UB';
CONFIG.LOGPATH_UB = [CONFIG.LOGPATH, CONFIG.LOGNAME_UB, '.log'];

% file name of the Wasserstein-2 optimal transport based upper bounds
CONFIG.LOGNAME_W2OTUB = 'W2OTUB';
CONFIG.LOGPATH_W2OTUB = [CONFIG.LOGPATH, CONFIG.LOGNAME_W2OTUB, '.log'];

% file name of the logs for evaluating other methods
CONFIG.LOGNAME_EVALUATION = 'evaluation';
CONFIG.LOGPATH_EVALUATION = [CONFIG.LOGPATH, CONFIG.LOGNAME_EVALUATION, '.log'];

end

