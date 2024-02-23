function CONFIG = BM1D_config()
% Place to store global configurations of the one-dimensional benchmark
% problem
% Output:
%   CONFIG: a struct containing configurations as fields

CONFIG = global_config();

% path to save files
CONFIG.SAVEPATH = [CONFIG.SAVEPATH_ROOT, ...
    'MT/Benchmark1DCPWA_Exp/'];

% path to log files
CONFIG.LOGPATH = [CONFIG.LOGPATH_ROOT, ...
    'MT/Benchmark1DCPWA_Exp/'];

% root folder for gurobi log files
CONFIG.LOGPATH_GUROBI = [CONFIG.LOGPATH, 'Gurobi/'];

% if the directory does not exist, create it first
if ~exist(CONFIG.SAVEPATH, 'dir')
    mkdir(CONFIG.SAVEPATH);
end

% if the directory does not exist, create it first
if ~exist(CONFIG.LOGPATH, 'dir')
    mkdir(CONFIG.LOGPATH);
end

% if the directory does not exist, create it first
if ~exist(CONFIG.LOGPATH_GUROBI, 'dir')
    mkdir(CONFIG.LOGPATH_GUROBI);
end

% file name (template) of the inputs
CONFIG.FILENAME_INPUTS = 'inputs_T%02d';
CONFIG.SAVEPATH_INPUTS = ...
    [CONFIG.SAVEPATH, CONFIG.FILENAME_INPUTS, '.mat'];

% file name (template) of the outputs
CONFIG.FILENAME_OUTPUTS = 'outputs_T%02d_M%03d';
CONFIG.SAVEPATH_OUTPUTS = ...
    [CONFIG.SAVEPATH, CONFIG.FILENAME_OUTPUTS, '.mat'];

% file name of the summary
CONFIG.FILENAME_SUMMARY = 'summary';
CONFIG.SAVEPATH_SUMMARY = ...
    [CONFIG.SAVEPATH, CONFIG.FILENAME_SUMMARY, '.mat'];

% file name of the main logs
CONFIG.LOGNAME_MAIN = 'main';
CONFIG.LOGPATH_MAIN = ...
    [CONFIG.LOGPATH, CONFIG.LOGNAME_MAIN, '.log'];

% file name of the linear semi-infinite programming logs
CONFIG.LOGNAME_LSIP_MAIN = 'LSIP_main';
CONFIG.LOGPATH_LSIP_MAIN = ...
    [CONFIG.LOGPATH, CONFIG.LOGNAME_LSIP_MAIN, '.log'];

% file name of the linear programming logs
CONFIG.LOGNAME_LSIP_LP = 'LSIP_LP';
CONFIG.LOGPATH_LSIP_LP = ...
    [CONFIG.LOGPATH_GUROBI, CONFIG.LOGNAME_LSIP_LP, '.log'];

% file name of the global optimization logs
CONFIG.LOGNAME_LSIP_GLOBAL = 'LSIP_global';
CONFIG.LOGPATH_LSIP_GLOBAL = ...
    [CONFIG.LOGPATH_GUROBI, CONFIG.LOGNAME_LSIP_GLOBAL, '.log'];

end

