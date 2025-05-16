% Print the objectives and the sub-optimalities of the approximate Wasserstein barycenters computed by our algorithm and other
% algorithms

CONFIG = WB_General_config();

load(CONFIG.SAVEPATH_INPUTS);

PSWB = load(CONFIG.SAVEPATH_EVALUATION_PSWBSTAIBCLAICISOLOMONJEGELKA);
NWB = load(CONFIG.SAVEPATH_EVALUATION_NWBFANTAGHVAEICHEN);
CW2B = load(CONFIG.SAVEPATH_EVALUATION_CW2BKOROTINLISOLOMONBURNAEV);
WIN = load(CONFIG.SAVEPATH_EVALUATION_WINKOROTINEGIAZARIANLIBURNAEV);
MMOT = load(CONFIG.SAVEPATH_EVALUATION_MMOTNEUFELDXIANG);
OURS = load(CONFIG.SAVEPATH_EVALUATION_OURALGO);


WB_LB = MT_LB_list(end);

fprintf('%25s%20s%20s\n', 'Algorithm', 'Objective', 'Sub-optimality')
fprintf('%25s%20.6f%20.4e\n', 'Staib et al. (2017)', PSWB.objective, PSWB.objective - WB_LB);
fprintf('%25s%20.6f%20.4e\n', 'Fan et al. (2021)', NWB.objective, NWB.objective - WB_LB);
fprintf('%25s%20.6f%20.4e\n', 'Korotin et al. (2021)', CW2B.objective, CW2B.objective - WB_LB);
fprintf('%25s%20.6f%20.4e\n', 'Korotin et al. (2022)', WIN.objective, WIN.objective - WB_LB);
fprintf('%25s%20.6f%20.4e\n', 'Neufeld and Xiang (2022)', MMOT.objective, MMOT.objective - WB_LB);
fprintf('%25s%20.6f%20.4e\n', 'our algorithm', OURS.objective, OURS.objective - WB_LB);