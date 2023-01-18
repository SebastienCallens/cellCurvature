# cellCurvature
A collection of Matlab tools for analysing cell-curvature interactions, used in the publication "Emergent collective organization of bone cells in complex curvature fields", Callens et al., (bioRxiv, 2020)

# Required software tools
Matlab (these codes were tested in 2018b and 2020a). Codes should work as soon as the software folder is copied in the working directory of Matlab.
Running times for all codes on a standard desktop computer should not take more than a couple of seconds (at least below one minute). A set of example data is included. 

# Overview of the software folder

- example_data: data folder containing image data for a set of experiments (D8 convex structures & D5 convex unduloid for nuclei)
- additional MAT-files: curvature maps and principal direction maps
- Matlab scripts:
	- frequencyMap: generates frequency maps of actin signal of periodic units 
	- intensityVSCurvatureMap: creates heat map of intensity vs k1 and k2, uses associated excel data sheets as input 
	- intensityVSDistanceBatch: creates curve of intensity vs distance to k2<0 
	- midlineIntensity: creates profile of actin intensity along the center line of a substrate 
	- orientationROI: example script to compute SF orientation in a user-defined ROI 
	- orientationROI_vs_PD: example script to compute SF orientation in a user-defined ROI and compare it to the principal directions 
	- runx2IntBatch: computes normalized mean DAPI and RUNX2 intensities in masked images of the convex unduloid substrate (D8), at predefined ROIs 
	- sheetSagging: script to compute vertical displacement of cell sheet and anchor density on spherical or cylindrical substrates. Requires user-input to draw mask.
	- nucleusHeat: script to compute density maps and plots of density vs distance to k2<0 for nuclei centroids. 
- Matlab functions are called in the associated scripts
