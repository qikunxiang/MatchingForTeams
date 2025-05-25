# Feasible approximation of matching equilibria for large-scale matching for teams problems
+ By Ariel Neufeld and Qikun Xiang
+ Article link (arXiv): https://arxiv.org/abs/2308.03550

## Description of files

+ **cutplane/** contains classes related to the cutting-plane algorithms used for solving linear semi-infinite programming problems

	- **LSIPMinCuttingPlaneAlgo.m**: abstract class for cutting-plane algorithms used for solving linear semi-infinite programming (minimization) problems

	- **MatchTeam1DCPWA.m**: class for the matching for teams problem with one-dimensional type spaces

	- **MatchTeam2DBusinessLocation.m**: class for the business location problem with two-dimensional type spaces and two-dimensional quality space

	- **MatchTeam2DWassersteinBarycenter.m**: class for the 2-Wasserstein barycenter problem with two-dimensional input measures solved via the parametric formulation based approach

	- **OTDiscrete.m**: class for the classical two-marginal discrete optimal transport problem solved via a cutting-plane/constraint-generation algorithm

+ **probmeas/** contains classes and functions related to probability measures
	- **HasTractableQuadraticIntegrals.m**: abstract class for providing two function interfaces for computing the mean vector and the covariance matrix of a probability measure

	- **ProbMeas1DCPWADens.m**: class for probability measures supported in a one-dimensional compact interval with continuous piece-wise affine density functions

	- **ProbMeas1DInterval.m**: abstract class for probability measures supported in a one-dimensional compact interval

	- **ProbMeas1DMixNorm.m**: class for one-dimensional mixture of Gaussian probability measures truncated to a one-dimensional compact interval 
 	
	- **ProbMeas2DAffDens.m**: class for probability measures supported in a two-dimensional convex polytope with affine density functions
 	
	- **ProbMeas2DConvexPolytope.m**: abstract class for probability measures supported in a two-dimensional convex polytope
 	
	- **ProbMeas2DCPWADens.m**: class for probability measures supported in a union of triangles with continuous piece-wise affine (CPWA) density functions
 	
	- **ProbMeas2DMixNorm.m**: class for probability measures supported in a union of triangles with truncated mixture of Gaussians density functions

+ **mex/** contains the C++ code used when solving the 2-Wasserstein barycenter problem in two dimensions

    - **Makefile**: used for compiling the C++ code into mex files for MATLAB

    - **mesh\_intersect\_power\_diagram.m**: empty MATLAB file for creating the MATLAB interface for the mex function mesh\_intersect\_power\_diagram; this function computes a power diagram corresponding to a collection of circles as well as the polygonal complex formed by the intersection of this power diagram and a given triangular mesh; it is used when computing semi-discrete 2-Wasserstein optimal transport for the class ProbMeas2DCPWADens

	- **minkowski\_sum\_2d.m**: empty MATLAB file for creating the MATLAB interface for the mex function minkowski\_sum\_2d; this function computes the Minkowski sum of two-dimensional polygons

	- **power\_diagram\_intersection.m**: empty MATLAB file for creating the MATLAB interface for the mex function power\_diagram\_intersection; this function computes power diagrams corresponding to multiple collections of circles as well as the polygonal complex formed by their intersections

    - **src/** contains the C++ source code

		- **mesh\_intersect\_power\_diagram.hpp**: C++ header file for the function mesh\_intersect\_power\_diagram

		- **mesh\_intersect\_power\_diagram.cpp**: C++ source file containing the function mesh\_intersect\_power\_diagram

        - **power\_diagram\_intersection.hpp**: C++ header file for the function power\_diagram\_intersection

        - **power\_diagram\_intersection.cpp**: C++ source file containing the function power\_diagram\_intersection

        - **mex\_wrapper\_mesh\_intersect\_power\_diagram.cpp**: mex wrapper class for the C++ function mesh\_intersect\_power\_diagram implemented in **power\_diagram\_intersection.cpp**

		- **mex\_wrapper\_minkowski\_sum\_2d.cpp**: mex wrapper class for the C++ function minkowski\_sum\_2d implemented in **minkowski\_sum\_2d.cpp**

        - **mex\_wrapper\_power\_diagram\_intersection.cpp**: mex wrapper class for the C++ function power\_diagram\_intersection implemented in **power\_diagram\_intersection.cpp**

        - **test\_power\_diagram\_intersection.cpp**: C++ code for testing the function power\_diagram\_intersection

		- **test\_minkowski\_sum\_2d.cpp**: C++ code for testing the function minkowski\_sum\_2d

        - **test\_mesh\_intersect\_power\_diagram.cpp**: C++ code for testing the function mesh\_intersect\_power\_diagram

+ **exp/** contains the scripts to run the numerical experiments (see below for detailed instructions)
	- **BusinessLocation\_EXP/**: contains all the relevant scripts used for the business location distribution experiment

	- **WassersteinBarycenter\_Exp/**: contains all the relevant scripts for the 2-Wasserstein barycenter experiment with general input measures

	- **WassersteinBarycenter\_Elliptical\_Exp/**: contains all the relevant scripts used for the 2-Wasserstein barycenter experiment where the input measures are from the same elliptical family

	- **Benchmark1DCPWA\_Exp/**: contains all the relevant scripts used for the experiment involving one-dimensional type measures

+ **utils/** contains external libraries and auxiliary functions
	- **tight\_subplot/**: used for creating figures with narrow margins

	- **check\_triangulation\_simplicial\_cover.m**: function used for checking whether a collection of triangles satisfy some conditions

	- **cleanup\_trangles.m**: function that removes triangles in a triangulation that are close to degeneracy

	- **comonotone\_coupling.m**: function that computes the comonotone coupling (i.e., coupling formed with the comonotone copula) of one-dimensional probability measures

	- **discrete\_OT.m**: function that computes an optimal coupling between two discrete probability measures via linear programming

	- **discrete\_reassembly\_1D.m**: function that computes a reassembly of a multi-dimensional discrete probability measure with one-dimensional discrete marginals

	- **discrete\_reassembly.m**: function that computes a reassembly of a discrete probabilty measure with discrete marginals

	- **WB\_elliptical\_fixedpoint.m**: auxiliary function used for (approximately) computing the 2-Wasserstein barycenter of probability measures that belong to the same elliptical family (see Álvarez-Esteban, P. C., Del Barrio, E., Cuesta-Albertos, J. A., & Matrán, C. (2016). A fixed-point approach to barycenters in Wasserstein space. *Journal of Mathematical Analysis and Applications*, 441(2), 744-762.)

# Instructions to run the numerical experiments

## Configurations

+ All folders and subfolders must be added to the MATLAB search path. 
+ Gurobi optimization (version 12.0.0 or above) must be installed on the machine and relevant files must be added to the search path. 
+ The following lines in the file **exp/global_config.m** can be edited to change the directories for the save files and log files if necessary.

		% root folder for all save files
		CONFIG.SAVEPATH_ROOT = '~/';

		% root folder for all log files
		CONFIG.LOGPATH_ROOT = '~/Logs/';

## Experiment: business location distribution problem

+ All the relevant scripts used in this experiment are located in **exp/BusinessLocation\_Exp/**.

### Step 1: generate the input file
+ Run **BizLoc\_prepare.m** to generate a file containing the setting of the business location distribution problem as well as the test functions used in the approximation schemes.

+ Run **BizLoc\_plot\_density.m** to plot the density functions of the input measures. 

### Step 2: compute the approximate matching equilibria

+ Run **BizLoc\_OT.m** to generate a file containing the information characterizing optimal couplings between the continuous probability measures and their discrete approximations.

+ Run **BizLoc\_run.m** to solve the linear semi-infinite programming problem via the cutting-plane algorithm. 

+ Run **BizLoc\_run\_UB.m** to compute the upper bounds via Monte Carlo integration. 

+ Run **BizLoc\_run\_transfuncs.m** to evaluate the computed approximately optimal transfer functions at a grid in the quality space. 

### Step 3: visualize the results

+ Run **BizLoc\_plot\_bounds.m** to plot the computed upper and lower bounds, sub-optimality estimates, and their a priori theoretical upper bounds.

+ Run **BizLoc\_plot\_business.m** to plot the computed approximately optimal business location distributions. 

+ Run **BizLoc\_run\_couplings.m** to generate random samples from the computed approximately optimal couplings for plotting. 

+ Run **BizLoc\_plot\_couplings.m** to plot the computed approximately optimal couplings.

+ Run **BizLoc\_plot\_transfuncs.m** to plot the computed approximately optimal transfer functions.

## Experiment: 2-Wasserstein barycenter

+ All the relevant files used in this experiment are located in **exp/WassersteinBarycenter_Exp/**.

+ Part of the code uses mex functions which need to be compiled from C++. The Makefile for the compilation process only supports macOS (both Intel or Apple Silicon are supported). The Makefile needs to be modified for other platforms.

+ The Computational Geometry Algorithms Library (CGAL) must be installed. The installation can be done either directly from the [GitHub repository](https://github.com/CGAL/cgal/), or via a package manager such as [Homebrew](https://brew.sh) or [Conda](https://anaconda.org/anaconda/conda).

+ Uncomment the line with `MEX_SUFFIX` in the following part of **mex/Makefile** depending on whether the machine has an Intel or Apple Silicon chip.

		# Suffix of mex files on this operating system
		# macOS with Apple Silicon
		#MEX_SUFFIX = mexmaca64
		
		# macOS with Intel
		#MEX_SUFFIX = mexmaci64


+ Modify the value of the variable `MATLAB_PATH` in **mex/Makefile** to the location of the MATLAB app on the machine.

		# Path to MATLAB
		MATLAB_PATH = /Applications/MATLAB_R2024b.app
		
+ Modify the following lines in **mex/Makefile** to the locations of the header files and the library files of CGAL (e.g., under `/usr/local`, `/opt/homebrew`, or `/opt/miniconda3`)

		# Include path (for header files)
		INCLUDE_PATH = /usr/local/include

		# Library path (for linking to libraries)
		LIBRARY_PATH = /usr/local/lib

### Step 0: compile the mex functions
+ In Terminal, use `cd` to change the current working directory to the **mex/** directory.
+ Execute the following commands in Terminal to build the mex files.

		mkdir build
		make

### Step 1: generate the input file
+ Run **WB\_General\_prepare.m** to generate a file containing the input measures as well as the test functions used in the approximation schemes.

+ Run **WB\_General\_plot\_density.m** to plot the density functions of the input measures. 

### Step 2: compute the approximate 2-Wasserstein barycenter
+ Run **WB\_General\_OT.m** to generate a file containing the information characterizing optimal couplings between the continuous probability measures and their discrete approximations.

+ Run **WB\_General\_run\_LB.m** to solve the linear semi-infinite programming problem via the cutting-plane algorithm and compute the lower bounds.

+ Run **WB\_General\_run\_UB.m** to compute the upper bounds via Monte Carlo integration.


### Step 3: visualize the results
+ Run **WB\_General\_plot\_bounds.m** to plot the computed upper and lower bounds, sub-optimality estimates, and their a priori theoretical upper bounds.

+ Run **WB\_General\_plot\_histogram.m** to plot the histograms of the computed approximate 2-Wasserstein barycenters.

### Step 4: compute the approximate 2-Wasserstein barycenter via 2-Wasserstein optimal couplings (optional)
+ Run **WB\_General\_W2OT.m** to generate a file containing the information characterizing 2-Wasserstein optimal couplings between the continuous probability measures and the discrete approximate barycenter.

+ Run **WB\_General\_run\_W2OTUB.m** to compute the upper bounds via Monte Carlo integration.


### Step 5: visualize the results via 2-Wasserstein optimal couplings (optional)
+ Run **WB\_General\_plot\_bounds\_with\_W2OT.m** to plot the upper and lower bounds, sub-optimality estimates computed via 2-Wasserstein optimal couplings, and their a priori theoretical upper bounds.

+ Run **WB\_General\_plot\_histogram\_with\_W2OT.m** to plot the histograms of the approximate 2-Wasserstein barycenters computed via 2-Wasserstein optimal couplings.


## Experiment: 2-Wasserstein barycenter of measures from the same elliptical family

+ All the relevant files used in this experiment are located in **exp/WassersteinBarycenter\_Elliptical\_Exp/**.

+ Part of the code uses mex functions which need to be compiled from C++. The Makefile for the compilation process only supports macOS (both Intel or Apple Silicon are supported). The Makefile needs to be modified for other platforms.
The configuration is identical to the **Experiment: 2-Wasserstein barycenter** section above.

### Step 0: compile the mex functions
+ In Terminal, use `cd` to change the current working directory to the **mex/** directory.
+ Execute the following commands in Terminal to build the mex files.

		mkdir build
		make

### Step 1: generate the input file
+ Run **WB\_prepare.m** to generate a file containing the input measures as well as the test functions used in the approximation schemes.

+ Run **WB\_Elliptical\_fixedpoint.m** to compute the "true" 2-Wasserstein barycenter via the fixed-point scheme of Álvarez-Esteban, P. C., Del Barrio, E., Cuesta-Albertos, J. A., & Matrán, C. (2016).

+ Run **WB\_Elliptical\_plot\_density.m** to plot the density functions of the "true" 2-Wasserstein barycenter and the input measures. 

### Step 2: approximate the matching for teams problem

+ Run **WB\_Elliptical\_OT.m** to generate a file containing the information characterizing optimal couplings between the continuous probability measures and their discrete approximations.

+ Run **WB\_Elliptical\_run\_LB.m** to solve the linear semi-infinite programming problem via the cutting-plane algorithm and compute the lower bounds

+ Run **WB\_Elliptical\_run\_UB.m** to compute the upper bounds via Monte Carlo integration. 


### Step 3: visualize the results

+ Run **WB\_Elliptical\_plot\_bounds.m** to plot the computed upper and lower bounds, sub-optimality estimates, and their a priori theoretical upper bounds.

+ Run **WB\_Elliptical\_plot\_histogram.m** to plot the computed approximate 2-Wasserstein barycenters and the "true" 2-Wasserstein barycenter computed via the fixed-point scheme. 


## Experiment: matching with one-dimensional type spaces

+ All the relevant files used in this experiment are located in **exp/Benchmark1DCPWA\_Exp/**.

### Step 1: generate the input files
+ Run **BM1D\_prepare.m** to generate files containing the problem instances as well as the test functions used in the approximation schemes.

### Step 2: approximate the matching for teams problems
+ Run **BM1D\_run.m** to approximate the problem instances via the parametric formulation based approach. 

### Step 3: visualize the results
+ Run **BM1D\_process\_results.m** to compile all results into a single data file.

+ Run **BM1D\_print\_results.m** to print a table containing the computed sub-optimality estimates, the running time, and the sparsity of support. 
