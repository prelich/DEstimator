
Estimation of the Diffusion Constant from Intermittent Trajectories with Variable Position Uncertainties
Peter Relich, Mark Olah, Patrick Cutler, Keith Lidke

% REQUIRED SOFTWARE PACKAGES PACKAGES:
MATLAB (tested on 2015a)
Dip Image (www.diplib.org) - required for the DSim class.

---- CONTENTS ----
The software package contains 1) A Simulation class for generating photon emissions from a diffusing particle. 2) Estimation functions for returning likelihood distributions and associated quantities.  3) Figure scripts for generating the figures in the manuscript.
-------------------- 


--- MATLAB Implementation Read Me ----

Quick Start.
Unpack all files into a folder and add to your MATLAB path.

For Figure 1 --- run: Sample_Trajectory.m
For Figure 2 --- run: QuadraticApproximationPlot.m
For Figure 3 --- run: InfoValsPlot.m
For Figure 4 --- run: VarInfoPlot.m
For Figure 5 & 6 run: Trajectory_Sim.m   --- Warning: this simulation takes ~ 8 hours on a modern PC (~2012 architecture)

*** NOTE: Some hand manipulation was performed to align the axis labels on figures 3 & 4 so the scripts are not purely automated

-----------------------------------------------------------------------------------------------------------


File Organization:


Simulation Code:
****************
DSim.m - Class for the Diffusion Simulation with a lattice background.


Estimation Code:
****************
computeMLE.m -- Finds the MLE for the likelihood distribution from the recursive method.
LLH_recursive1D.m -- Implementation of the likelihood distribution (recursive method).
LLH_recursiveND.m -- Multi-dimensional wrapper for LLH_recursive1D for processing higher dimensional trajectories (2D and higher)
LLH_ObsI_recursive1D.m -- Modified recursive1D implementation that also returns the "observed information" at a given sampled D.
			*** NOTE: If the sampled D is not the MLE, the returned "observed information" is meaningless ***
recursiveFisher1D.m -- Returns the expected fisher information for a sampled D
waldInterval.m -- Returns the sigma deviation of an estimate given the observed information (function is invariant to the space of the observed information).


Figure Generation Code:
***********************
Sample_Trajectory.m
QuadraticApproximationPlot.m
InfoValsPlot.m
VarInfoPlot.m
Trajectory_Sim.m