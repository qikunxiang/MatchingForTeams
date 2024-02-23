# Feasible approximation of matching equilibria for large-scale matching for teams problems

+ By Ariel Neufeld and Qikun Xiang
+ Article link (arXiv): https://arxiv.org/abs/2308.03550

## Description of files

+ cutplane/      contains classes related to the cutting-plane algorithms used for solving linear semi-infinite programming problems
	- **LSIPMinCuttingPlaneAlgo.m**: abstract class for cutting-plane algorithms used for solving linear semi-infinite programming (minimization) problems
	- **MT2DBizLoc_ParTrans.m**: class for the business location problem involving two-dimensional type spaces and two-dimensional quality space solved via the parametric formulation based approach
   - **MT2DQuad.m**: abstract class for the matching for teams problem involving two-dimensional type spaces and two-dimensional quality space
   - **MT2DQuad\_MMOT.m**: class for the matching for teams problem involving two-dimensional type spaces and two-dimensional quality space solved via the multi-marginal optimal transport (MMOT) based approach
   - **MT2DQuad\_ParTrans.m**: class for the matching for teams problem involving two-dimensional type spaces and two-dimensional quality space solved via the parametric formulation based approach
   - **MT1DCPWA.m**: abstract class for the matching for teams problem involving one-dimensional type spaces 
   - **MT1DCPWA\_MMOT.m**: class for the matching for teams problem involving one-dimensional type spaces solved via the multi-marginal optimal transport (MMOT) based approach	
   - **MT1DCPWA\_ParTrans.m**: class for the matching for teams problem involving one-dimensional type spaces solved via the parametric formulation based approach
   - **OTDiscrete.m**: class for the classical two-marginal discrete optimal transport problem solved via a cutting-plane/constraint-generation algorithm

+ probmeas/		contains classes and functions related to probability measures
 	- **ProbMeas2D\_ConvexPolytope.m**: abstract class for probability measures supported in a two-dimensional convex polytope
 	- **HasTractableQuadraticIntegrals**: abstract class for providing two function interfaces for computing the mean vector and the covariance matrix of a probability measure
 	- **ProbMeas2D\_AffDens.m**: class for probability measures supported in a two-dimensional convex polytope with affine density functions
 	- **ProbMeas2D\_CPWADens.m**: class for probability measures supported in a union of triangles with continuous piece-wise affine (CPWA) density functions
 	- **ProbMeas2D\_MixNorm.m**: class for probability measures supported in a union of triangles with truncated mixture of Gaussians density functions
 	- **ProbMeas1D\_Interval.m**: abstract class for probability measures supported in a one-dimensional compact interval
 	- **ProbMeas1D\_CPWADens.m**: class for probability measures supported in a one-dimensional with continuous piece-wise affine density functions
 	- **check\_triangulation\_simplicial\_cover.m**: function used for checking whether a collection of triangles satisfy some conditions
 	- **cleanup\_trangles.m**: function that removes triangles in a triangulation that are close to degeneracy
 	- **comonotone\_coupling.m**: function that computes the comonotone coupling (i.e., coupling formed with the comonotone copula) of one-dimensional probability measures
 	- **discrete\_OT.m**: function that computes an optimal coupling between two discrete probability measures via linear programming
 	- **discrete\_reassembly.m**: function that computes a reassembly of a discrete probabilty measure with discrete marginals
 	- **discrete\_reassembly\_1D.m**: function that computes a reassembly of a multi-dimensional discrete probability measure with one-dimensional discrete marginals
 	- **hyperbolic\_Dirichlet\_tessellation\_polar.m**: function used to express a cell in a hyperbolic Dirichlet tessellation via its polar form

+ experiments/            contains the scripts to run the numerical experiments (see below for detailed instructions)

+ utils/          contains external libraries and auxiliary functions
   - /tight\_subplot/: used for creating figures with narrow margins
   - **WB_locationscatter**: auxiliary function used for (approximately) computing the 2-Wasserstein barycenter of probability measures that belong to the same location-scatter family (see Álvarez-Esteban, P. C., Del Barrio, E., Cuesta-Albertos, J. A., & Matrán, C. (2016). A fixed-point approach to barycenters in Wasserstein space. *Journal of Mathematical Analysis and Applications*, 441(2), 744-762.)

## Instructions to run the numerical experiments

### Configurations

+ All folders and subfolders must be added to the MATLAB search path. 
+ Gurobi optimization (version 9.5.0 or above) must be installed on the machine and relevant files must be added to the search path. 
+ In the file **exp/global\_config.m**, the values of *CONFIG.SAVEPATH_ROOT* and *CONFIG.LOGPATH_ROOT* can be modified to change the directories where save files and log files are saved. 

### Experiment: Business location distribution problem

#### Step 1: generate the input file
+ Run **exp/BusinessLocation\_Exp/BizLoc\_prepare.m** to generate a file containing the setting of the business location distribution problem as well as the test functions used in the approximation schemes.
+ Run **exp/BusinessLocation\_Exp/BizLoc\_plot\_density.m** to plot the density functions of the input measures. 

#### Step 2: approximate the matching for teams problem
+ Run **exp/BusinessLocation\_Exp/BizLoc\_OT.m** to generate a file containing the information characterizing optimal couplings between the continuous probability measures and their discrete approximations.
+ Run **exp/BusinessLocation\_Exp/BizLoc\_run.m** to solve the linear semi-infinite programming problem via the cutting-plane algorithm. 
+ Run **exp/BusinessLocation\_Exp/BizLoc\_run\_UB.m** to compute the upper bounds via Monte Carlo integration. 
+ Run **exp/BusinessLocation\_Exp/BizLoc\_run\_transfuncs.m** to evaluate the computed approximately optimal transfer functions at a grid in the quality space. 

#### Step 3: visualize the results
+ Run **exp/BusinessLocation\_Exp/BizLoc\_plot\_bounds.m** to plot the computed upper and lower bounds, sub-optimality estimates, and their a priori theoretical upper bounds.
+ Run **exp/BusinessLocation\_Exp/BizLoc\_plot\_business.m** to plot the computed approximately optimal business location distributions. 
+ Run **exp/BusinessLocation\_Exp/BizLoc\_run\_couplings.m** to generate random samples from the computed approximately optimal couplings for plotting. 
+ Run **exp/BusinessLocation\_Exp/BizLoc\_plot\_couplings.m** to plot the computed approximately optimal couplings.
+ Run **exp/BusinessLocation\_Exp/BizLoc\_plot\_transfuncs.m** to plot the computed approximately optimal transfer functions.


### Experiment: 2-Wasserstein barycenter (location-scatter family)

#### Step 1: generate the input file
+ Run **exp/WassersteinBarycenter_Exp1/WB\_prepare.m** to generate a file containing the input measures as well as the test functions used in the approximation schemes.
+ Run **exp/WassersteinBarycenter_Exp1/WB\_plot\_density.m** to plot the density functions of the input measures. 

#### Step 2: approximate the matching for teams problem
+ Run **exp/WassersteinBarycenter_Exp1/WB\_OT.m** to generate a file containing the information characterizing optimal couplings between the continuous probability measures and their discrete approximations.
+ Run **exp/WassersteinBarycenter_Exp1/WB\_run.m** to solve the linear semi-infinite programming problem via the cutting-plane algorithm and compute the upper bounds via Monte Carlo integration. 


#### Step 3: visualize the results
+ Run **exp/WassersteinBarycenter_Exp1/WB\_fixedpoint.m** to compute the "true" 2-Wasserstein barycenter via the fixed-point scheme of Álvarez-Esteban, P. C., Del Barrio, E., Cuesta-Albertos, J. A., & Matrán, C. (2016).
+ Run **exp/WassersteinBarycenter_Exp1/WB\_plot\_bounds.m** to plot the computed upper and lower bounds, sub-optimality estimates, and their a priori theoretical upper bounds.
+ Run **exp/WassersteinBarycenter_Exp1/WB\_plot\_histogram.m** to plot the computed approximate 2-Wasserstein barycenters and the "true" 2-Wasserstein barycenter computed via the fixed-point scheme. 

### Experiment: 2-Wasserstein barycenter (general case)

#### Step 1: generate the input file
+ Run **exp/WassersteinBarycenter_Exp2/WB\_prepare.m** to generate a file containing the input measures as well as the test functions used in the approximation schemes.
+ Run **exp/WassersteinBarycenter_Exp2/WB\_plot\_density.m** to plot the density functions of the input measures. 

#### Step 2: approximate the matching for teams problem
+ Run **exp/WassersteinBarycenter_Exp2/WB\_OT.m** to generate a file containing the information characterizing optimal couplings between the continuous probability measures and their discrete approximations.
+ Run **exp/WassersteinBarycenter_Exp2/WB\_run.m** to solve the linear semi-infinite programming problem via the cutting-plane algorithm and compute the upper bounds via Monte Carlo integration. 


#### Step 3: visualize the results
+ Run **exp/WassersteinBarycenter_Exp2/WB\_plot\_bounds.m** to plot the computed upper and lower bounds, sub-optimality estimates, and their a priori theoretical upper bounds.
+ Run **exp/WassersteinBarycenter_Exp2/WB\_plot\_histogram.m** to plot the computed approximate 2-Wasserstein barycenters.


### Experiment 3

#### Step 1: generate the input files
+ Run **exp/Benchmark1DCPWA_Exp/BM1D\_prepare.m** to generate files containing the problem instances as well as the test functions used in the approximation schemes.

#### Step 2: approximate the matching for teams problems
+ Run **exp/Benchmark1DCPWA_Exp/BM1D\_run.m** to approximate the problem instances via the parametric formulation based approach. 

#### Step 3: visualize the results
+ Run **exp/Benchmark1DCPWA_Exp/BM1D\_process\_results.m** to compile all results into a single data file.
+ Run **exp/Benchmark1DCPWA_Exp/BM1D\_print\_results.m** to print a table containing the computed sub-optimality estimates, the running time, and the sparsity of support. 
+ Run **exp/Benchmark1DCPWA_Exp/BM1D\_plot\_time.m** to plot the running time of the linear programming solver and the global minimization oracle.
