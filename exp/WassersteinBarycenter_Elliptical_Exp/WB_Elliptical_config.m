function CONFIG = WB_Elliptical_config()
% Place to store global configurations of the Wasserstein barycenter example within an elliptical family
% Output:
%   CONFIG: a struct containing configurations as fields

CONFIG = global_config();

% path to save files
CONFIG.SAVEPATH = [CONFIG.SAVEPATH_ROOT, 'MT/WassersteinBarycenter_Elliptical_Exp/'];

% path to log files
CONFIG.LOGPATH = [CONFIG.LOGPATH_ROOT, 'MT/WassersteinBarycenter_Elliptical_Exp/'];

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
CONFIG.SAVEPATH_INPUTS = [CONFIG.SAVEPATH, CONFIG.FILENAME_INPUTS, '.mat'];

% file name of the outputs of the fixed-point scheme
CONFIG.FILENAME_OUTPUTS_FIXEDPOINT = 'outputs_fixedpoint';
CONFIG.SAVEPATH_OUTPUTS_FIXEDPOINT = [CONFIG.SAVEPATH, CONFIG.FILENAME_OUTPUTS_FIXEDPOINT, '.mat'];

% file name of the optimal transport outputs
CONFIG.FILENAME_OUTPUTS_OT = 'outputs_OT';
CONFIG.SAVEPATH_OUTPUTS_OT = [CONFIG.SAVEPATH, CONFIG.FILENAME_OUTPUTS_OT, '.mat'];

% file name of the outputs for lower bounds
CONFIG.FILENAME_OUTPUTS_LB = 'outputs_LB';
CONFIG.SAVEPATH_OUTPUTS_LB = [CONFIG.SAVEPATH, CONFIG.FILENAME_OUTPUTS_LB, '.mat'];

% file name of the outputs for upper bounds
CONFIG.FILENAME_OUTPUTS_UB = 'outputs_UB';
CONFIG.SAVEPATH_OUTPUTS_UB = [CONFIG.SAVEPATH, CONFIG.FILENAME_OUTPUTS_UB, '.mat'];

% file name of the optimal transport logs
CONFIG.LOGNAME_OT = 'OT';
CONFIG.LOGPATH_OT = [CONFIG.LOGPATH, CONFIG.LOGNAME_OT, '.log'];

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

end
