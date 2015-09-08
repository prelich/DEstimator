This is the help file for the DEstimator package.  Which accompanies the paper:

"Estimation of the Diffusion Constant from Intermittent Trajectories with Variable Localization Precision"

May 2015.

Peter K. Relich (physx.grad@gmail.com)_
Mark J. Olah (mjo@cs.unm.edu)
Partick J. Cutler
Keith A. Lidke (klidke@unm.edu)

Getting started using the DEstimator Package:

1. Unzip the DEstimator package into a folder that can be accessed by matlab.


2. In Matlab, go into the folder that the DEstimator package was unzipped to and run the file:

>> startupDEstimator.m

This will add the correct paths to matlab.


3. Create an instance of a DEstimator object in the matlab workspace:
>> est = DEstimator();

See also:

Class overview:
>> help DEstimator

Class constructor useage:
>> help DEstimator.DEstimator


4. Import your trajectory into the matlab workspace -- the trajectory should be in the format:
 i.)	Matrix of observed_positions (N x dim) where N is the number of localizations, dim is the dimension of the data (units are arbitrary distance units)
 ii.)	Vector of observation_times (N x 1), must be increasing. (Given in arbitrary time units)
 iii.)	Matrix of observation_variances (N x dim), the variance in the measurement (estimation) of the observation. (These are in dist_units^2)
 iv.)	Scalar exposureTime (dT), giving the exposure time of an individual frame (in same time units used for times vector)

See Also:
>> help DEstimator.DEstimator

Alternatively, the DEstimator class can create a simulated 1D trajectory:
>> [observed_positions, observation_times, observation_variances] = A.simulate1D(D,N,V,dT);

Please see the test script: test_DEstimator.m for example useage of this function.
See Also:
>> help DEstimator.simulate1D


5. The trajectory needs to be loaded into the the DEstimator object:
>> success = est.initializeTrack(observed_positions, observation_times, obeservation_variances, exposureTime)

If success = 1, then the track is loaded correctly.  An estimator can be reset 
to work with a different trajectory by calling est.initializeTrack() again.

See also:
>> help DEstimator.initializeTrack

6. Once a trajectory is loaded, the DEstimator object can compute the maximum likelihood estimate of D (D_mle) with the following 
command (using matlab's fminsearch):
>> D_mle = est.MLE();

See also:
>> help DEstimator.MLE

More generally, the log-likelihood can be computed for 1 or more user determined diffusion constants of interest (D_vals):
>> LLHvec = est.LLH(D_vals);

LLHvec is the log likelihood vector and D_vals is a user input vector of sampled D values.

See also:
>> help DEstimator.LLH

Please see the test_DEstimator.m script for a quick demonstration of this class.

For more details on the DEstimator class, see the file DEstimator.m which has extensive documentation, 
and many useful functions for simulating and plotting output and testing
the recursive, Markov, and Laplace methods using various implementations.

