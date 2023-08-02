# Feasible approximation of matching equilibria in the matching for teams problem beyond discrete measures

+ By Ariel Neufeld and Qikun Xiang
+ Article link (arXiv): [TBD]

## Description of files

+ cutplane/      contains classes related to the cutting-plane algorithms used for solving linear semi-infinite programming problems
	- **LSIPMinCuttingPlaneAlgo.m**: abstract class for cutting-plane algorithms used for solving linear semi-infinite programming (minimization) problems
   - **MT2DQuad.m**: abstract class for the matching for teams problem involving two-dimensional type spaces and two-dimensional quality space (studied in Experiment 1)
   - **MT2DQuad\_MMOT.m**: class for the matching for teams problem involving two-dimensional type spaces and two-dimensional quality space solved via the multi-marginal optimal transport (MMOT) based approach
   - **MT2DQuad\_ParTrans.m**: class for the matching for teams problem involving two-dimensional type spaces and two-dimensional quality space solved via the parametric formulation based approach
   - **MT1DCPWA.m**: abstract class for the matching for teams problem involving one-dimensional type spaces (studied in Experiment 2)
   - **MT1DCPWA\_MMOT.m**: class for the matching for teams problem involving one-dimensional type spaces solved via the multi-marginal optimal transport (MMOT) based approach	
   - **MT1DCPWA\_ParTrans.m**: class for the matching for teams problem involving one-dimensional type spaces solved via the parametric formulation based approach

+ probmeas/		contains classes and functions related to probability measures
 	- **ProbMeas2D\_ConvexPolytope.m**: abstract class for probability measures supported in a two-dimensional convex polytope
 	- **ProbMeas2D\_AffDens.m**: class for probability measures supported in a two-dimensional convex polytope with affine density functions
 	- **ProbMeas2D\_CPWADens.m**: class for probability measures supported in a two-dimensional convex polytope with continuous piece-wise affine (CPWA) density functions (used in Experiment 1)
 	- **ProbMeas1D\_Interval.m**: abstract class for probability measures supported in a one-dimensional compact interval
 	- **ProbMeas1D\_CPWADens.m**: class for probability measures supported in a one-dimensional with continuous piece-wise affine density functions (used in Experiment 2)
 	- **check\_triangulation\_simplicial\_cover.m**: function used for checking whether a collection of triangles satisfy some conditions
 	- **comonotone\_coupling.m**: function that computes the comonotone coupling (i.e., coupling formed with the comonotone copula) of one-dimensional probability measures
 	- **discrete\_OT.m**: function that computes an optimal coupling between two discrete probability measures via linear programming
 	- **discrete\_reassembly.m**: function that computes a reassembly of a discrete probabilty measure with discrete marginals
 	- **discrete\_reassembly\_1D.m**: function that computes a reassembly of a multi-dimensional discrete probability measure with one-dimensional discrete marginals
 	- **hyperbolic\_Dirichlet\_tessellation\_polar.m**: function used to express a cell in a hyperbolic Dirichlet tessellation via its polar form

+ experiments/            contains the scripts to run the numerical experiments (see below for detailed instructions)

+ utils/          contains external libraries
    - utils/tight\_subplot/:             used for creating figures with narrow margins

## Instructions to run the numerical experiments

### Configurations

+ All folders and subfolders must be added to the MATLAB search path. 
+ Gurobi optimization (version 9.5.0 or above) must be installed on the machine and relevant files must be added to the search path. 
+ Since running both Experiment 1 and Experiment 2 will generate around 36GB of data files, one needs to make sure that sufficient disk space is available. 

### Experiment 1

#### Step 1: generate the input file
+ Run **experiments/experiment1/exp1\_prepare.m** to generate a file containing the setting of the matching for teams problem as well as the test functions used in the approximation schemes.
+ Run **experiments/experiment1/exp1\_plot\_testfuncs.m** to plot the agent type spaces, the type distributions, the quality space, and the test functions. 

#### Step 2: compute optimal couplings
+ Run **experiments/experiment1/exp1\_compute\_OT.m** to generate a file containing the information characterizing optimal couplings between the continuous probability measures and their discrete approximations.
+ Run **experiments/experiment1/exp1\_plot\_Laguerre.m** to plot the Laguerre diagrams characterizing the computed optimal couplings. 

#### Step 3: approximate the matching for teams problem via the parametric formulation based approach
+ Run **experiments/experiment1/exp1\_ParTrans\_run.m** to solve the linear semi-infinite programming problem via the cutting-plane algorithm as well as compute the upper bounds via Monte Carlo integration. *Warning: this step will generate around 3.7GB of data files.*
+ Run **experiments/experiment1/exp1\_ParTrans\_plot\_bounds.m** to plot the computed upper and lower bounds, sub-optimality estimates, and their a priori theoretical upper bounds.
+ Run **experiments/experiment1/exp1\_ParTrans\_plot\_couplings.m** and **experiments/experiment1/exp1\_ParTrans\_plot\_funcs.m** to plot the computed approximately optimal transfer functions as well as the computed approximately optimal couplings.

### Experiment 2

#### Step 1: generate the input files
+ Run **experiments/experiment2/exp2\_prepare.m** to generate files containing the settings of the matching for teams problems as well as the test functions used in the approximation schemes.

#### Step 2: approximate the matching for teams problems
+ Run **experiments/experiment2/exp2\_ParTrans\_run.m** to approximate the matching for teams problems via the parametric formulation based approach. *Warning: this step will generate around 32.2GB of data files.*

#### Step 3: visualize the results
+ Run **experiments/experiment2/exp2\_process\_results.m** to compile all results into a single data file.
+ Run **experiments/experiment2/exp2\_print\_results.m** to print a table containing the computed sub-optimality estimates, the running time, and the sparsity of support. 
+ Run **experiments/experiment2/exp2\_plot\_time.m** to plot the running time of the linear programming solver and the global minimization oracle.
