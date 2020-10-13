These are scripts and MATLAB functions that we used to compute reliability estimates in Compere et al (2020).

These documents are distributed exclusively by the PICAN group Greg Siegle, Ph.D. (gsiegle@pitt.edu) Director. Laurie Compere (ljc44@pitt.edu), repository author. Program Contact: pican@upmc.edu

==============================================
Disclaimer: - They are our group's code made to work with our idiosyncratic datasets. They were not made with public consumption in mind. Rather, they were functional for us. They will likely require a bunch of modification to work with your data. That said, we are assuming that having our code is easier for you than starting from scratch. 

The scripts assume the data is in AFNI or .nii format. 
The scripts that start with rl_* are scripts that you will need to modify to run on your own data. 
The others are functions called by those scripts. 
Note - you'll also need to download the afni MATLAB scripts to load afni data into matlab, from here: https://sscc.nimh.nih.gov/afni/matlab/ and add that folder to your path.

================================================
The following scripts should be run in order for activation tasks:

- rl_cbtsretssregs_auc.csh will compute regressions at the level subject with area under the curve (via multiple regression of a delta function across 8 TRs using 3dDeconvolve, with sum of betas as the parameter retained) and peak amplitude from the same 
regressions

- rl_cbtsretssregs_canonical_amplitude.csh compute regressions at the level subject with amplitude of a canonically shaped BOLD signal (using AFNI’s 3dDeconvolve with a tent function)

- rl_compute_gamma_variates.m compute regressions at the level subject with a gamma variate model with parameters for onset-delay, rise-decay rate,  and height.

- rl_compute_ICCs.m will compute ICCs and semi partial correlations with and without covariates map at the second level and store them in a corr structure saved in the current directory

- rl_write_ICC_maps_computed_with_matlab.m will write those maps under a format of .BRIK & .HEAD used by afni

- rl_choose_your_reliability_threshold.m will then read those .BRIK and .HEAD maps and count how many voxels are reliable in each ROI at different thresholds (whithout cluster correction) and compute mean, median and standard deviation values of ICC per ROI and display this information in the command window

- rl_which_model_best.m will also read those .BRIK and .HEAD maps and compare distributions of semi partial correlations reliability estimates provided by each model to determine which models provide significantly better reliability

- rl_reliabilityclusters_adjusted.csh will compute the number of voxel appropriate for cluster correction on each ROI and each first level measure's ICC map and create masks to apply cluster correction in each case

- rl_choose_your_reliability_threshold_cluster_corrected_a.m will read the ICC maps and count how many voxels are reliable in each ROI applying cluster correction adjusted at each first level measure and reliability threshold and display this information in the command window.

================================================
The following scripts should be run in order for rtfMRI-nf task:

-rl_nf_preprocessing.csh will run two kind of preprocessing on the raw data: one preprocessing emulating the real-time processing run by Turbo BrainVoyager called "BV style", and a more standard preprocessing usually run when analyzing our data post hoc, called "standard".

- rl_nf_preprocessing_nf_data.m will emulate the analysis done in real time by Turbo BrainVoyager on the data preprocessed with the "BV style" pipeline. It will take out the first 3 dummy scans, keep only the time course in voxels during happy blocks minus the mean value of the previous rest block adn create stim files for data preprocessed with the BS style and standard pipelines independantly.

- rl_nf_ssregs_auc.csh will compute regressions at the level subject with area under the curve (via multiple regression of a delta function across 8 TRs using 3dDeconvolve, with sum of betas as the parameter retained) and peak amplitude from the same 
regressions.

- rl_nf_ssregs_canonical_amplitude.csh compute regressions at the level subject with amplitude of a canonically shaped BOLD signal (using AFNI’s 3dDeconvolve with a tent function).

- rl_nf_compute_gamma_variates.m compute regressions at the level subject with a gamma variate model with parameters for onset-delay, rise-decay rate,  and height.

- rl_nf_compute_ICCs.m will compute ICCs and semi partial correlations map with and without covariates at the second level and store them in a corr structure saved in the current directory

- rl_nf_write_ICC_maps_computed_with_matlab.m will write those maps under a format of .BRIK & .HEAD used by afni

- rl_nf_choose_your_reliability_threshold.m will then read those .BRIK and .HEAD maps and count how many voxels are reliable in each ROI at different thresholds and compute  mean, standard deviation and median values of ICC per ROI and display this information in the command window

- rl_nf_which_model_best.m will also read those .BRIK and .HEAD maps and compare distributions of semi partial correlations reliability estimates provided by each model to determine which models provide significantly better reliability

- rl_nf_reliabilityclusters_adjusted.csh will compute the number of voxel appropriate for cluster correction on each ROI and each first level measure's ICC map and create masks to apply cluster correction in each case

- rl_nf_choose_your_reliability_threshold_cluster_corrected_a.m will read the ICC maps and count how many voxels are reliable in each ROI applying cluster correction adjusted at each first level measure and reliability threshold and display this information in the command window.

================================================
Copyleft and Legalese

The software is distributed under copyleft protections such that we request all code in this repository including this readme document  be distributed with it. We further request that citations to this work be as:
  Compere, L., Siegle, G.J., Young, K. (2019). Software for computing fMRI reliability.

Permission to use this software is granted subject to the following restrictions and understandings:
1) There is no warrantee or statement that the operation of this software will be error free. The authors and the University of Pittsburgh are under no obligation to provide any services by way of maintenance, update or otherwise.
2) Any user of such software agrees to indemnify and hold harmless the authors and the University of Pittsburgh from all claims arising out of the use of this software or arising out of any accident, injury or damage and from all costs, counsel fees, and liabilities incurred in or about any such claim, action, or proceeding brought thereon
3) As the software is still experimental, Greg Siegle requests that he be consulted on projects using the software.
4) Users are requested to inform Greg Siegle of noteworthy uses of this software.
