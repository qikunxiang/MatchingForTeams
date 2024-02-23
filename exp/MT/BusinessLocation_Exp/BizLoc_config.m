function CONFIG = BizLoc_config()
% Place to store global configurations of the business outlet locations
% experiment
% Output:
%   CONFIG: a struct containing configurations as fields

CONFIG = global_config();

% path to save files
CONFIG.SAVEPATH = [CONFIG.SAVEPATH_ROOT, ...
    'MT/BusinessLocation_Exp/'];

% path to log files
CONFIG.LOGPATH = [CONFIG.LOGPATH_ROOT, ...
    'MT/BusinessLocation_Exp/'];

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

% file name of the inputs
CONFIG.FILENAME_INPUTS = 'inputs';
CONFIG.SAVEPATH_INPUTS = ...
    [CONFIG.SAVEPATH, CONFIG.FILENAME_INPUTS, '.mat'];

% file name of the outputs
CONFIG.FILENAME_OUTPUTS = 'outputs_LSIP';
CONFIG.SAVEPATH_OUTPUTS = ...
    [CONFIG.SAVEPATH, CONFIG.FILENAME_OUTPUTS, '.mat'];

% file name of the optimal transport outputs
CONFIG.FILENAME_OT = 'OT';
CONFIG.SAVEPATH_OT = ...
    [CONFIG.SAVEPATH, CONFIG.FILENAME_OT, '.mat'];

% file name of the upper bounds
CONFIG.FILENAME_OUTPUTS_UB = 'outputs_UB';
CONFIG.SAVEPATH_OUTPUTS_UB = ...
    [CONFIG.SAVEPATH, CONFIG.FILENAME_OUTPUTS_UB, '.mat'];

% file name of the transfer functions
CONFIG.FILENAME_OUTPUTS_TRANSFUNCS = 'outputs_transfuncs';
CONFIG.SAVEPATH_OUTPUTS_TRANSFUNCS = ...
    [CONFIG.SAVEPATH, CONFIG.FILENAME_OUTPUTS_TRANSFUNCS, '.mat'];

% file name of the couplings for plotting
CONFIG.FILENAME_COUPLINGS = 'couplings';
CONFIG.SAVEPATH_COUPLINGS = ...
    [CONFIG.SAVEPATH, CONFIG.FILENAME_COUPLINGS, '.mat'];

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

% file name of the optimal transport logs
CONFIG.LOGNAME_OT = 'OT';
CONFIG.LOGPATH_OT = ...
    [CONFIG.LOGPATH, CONFIG.LOGNAME_OT, '.log'];

% file name of the upper bounds logs
CONFIG.LOGNAME_UB = 'UB';
CONFIG.LOGPATH_UB = ...
    [CONFIG.LOGPATH, CONFIG.LOGNAME_UB, '.log'];

% file name of the transfer functions logs
CONFIG.LOGNAME_TRANSFUNCS = 'transfuncs';
CONFIG.LOGPATH_TRANSFUNCS = ...
    [CONFIG.LOGPATH, CONFIG.LOGNAME_TRANSFUNCS, '.log'];

end

