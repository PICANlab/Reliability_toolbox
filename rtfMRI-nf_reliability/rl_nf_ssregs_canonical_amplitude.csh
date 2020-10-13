#!/bin/csh
# does single subject regressions with modelling a canonical hemodynamic response for the neurofeedback task
# as implemented on the dataset from Young et al (2017) AJP
#   Note that this paper combined data from two groups, and as such there is code below
#   for integrating across the 2 groups
# This script assumes we have the AFNI library in our path
# This code has been described in Compere et al. (2020)
# the script assumes we have regressors made which are 1 at trial onsets and 0 elsewhere for a given condition
# e.g., 
#  1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0....
# and also censor files with entries for every TR in which 0=censored. 1=uncensored
# (our scripts were via create_censor_file.m)

set fdir = ssregs_canonical_amplitude
cd data_for_voxel_wise_reliability_right_data/Active
# BV style analysis
foreach sub(`ls`)
    cd $sub/Preprocessing_BV_style
    set functional = Preprocessed_all_runs.$sub+tlrc
    foreach visit(Visit_1 Visit_2)
        cd $visit/$sub.results

        3dDeconvolve \
        -input $functional \
        -mask ../../../../../../Lamygdala_resampled+tlrc \
        -num_stimts 1 \
        -nfirst 1 \
        -stim_times 1  '1D: 0 40 80 120' 'BLOCK(40,1)' -stim_label 1 happy \
        -polort 0 \
        -fout -rout -tout \
        -full_first \
        -bucket $fdir/`echo $functional| cut -c 23-27`_BV_style-happy_canonical.ssreg \
        -GOFORIT 2
        cd ../..
    end
    cd ../Preprocessing_standard
    set functional = all_runs.$sub+tlrc
    foreach visit(Visit_1 Visit_2)
        cd $visit/$sub.results
        ## Normal post hoc analysis
        3dDeconvolve \
        -input $functional \
        -mask ../../../../../../Lamygdala_resampled+tlrc \
        -num_stimts 3 \
        -nfirst 1 \
        -stim_times 1 '1D: 40 160 280 400' 'BLOCK(40,1)' -stim_label 1 happy \
        -stim_times 2 '1D: 0 120 240 360 480' 'BLOCK(40,1)' -stim_label 2 rest \
        -stim_times 3 '1D: 80 200 320 440' 'BLOCK(40,1)' -stim_label 3 count \
        -censor "censor_file.1D" \
        -ortvec "dfile_rall.1D" motion \
        -polort 0 \
        -fout -rout -tout \
        -full_first \
        -bucket $fdir/`echo $functional| cut -c 10-14`_standard-happy_canonical.ssreg \
        -GOFORIT 2
        cd ../..
    end
cd ../..
end
cd ../Control
foreach sub(`ls`)
# BV style analysis
    cd $sub/Preprocessing_BV_style
    set functional = Preprocessed_all_runs.$sub+tlrc
    foreach visit(Visit_1 Visit_2)
        cd $visit/$sub.results
        3dDeconvolve \
        -input $functional \
        -mask ../../../../../../Lamygdala_resampled+tlrc \
        -num_stimts 1 \
        -nfirst 1 \
        -stim_times 1  '1D: 0 40 80 120' 'BLOCK(40,1)' -stim_label 1 happy \
        -polort 0 \
        -fout -rout -tout \
        -full_first \
        -bucket $fdir/`echo $functional| cut -c 23-27`_BV_style-happy_canonical.ssreg \
        -GOFORIT 2
        cd ../..
    end
    cd ../Preprocessing_standard
    set functional = all_runs.$sub+tlrc
    foreach visit(Visit_1 Visit_2)
        cd $visit/$sub.results
        ## Normal post hoc analysis
        3dDeconvolve \
        -input $functional \
        -mask ../../../../../../Lamygdala_resampled+tlrc \
        -num_stimts 3 \
        -nfirst 1 \
        -stim_times 1 '1D: 40 160 280 400' 'BLOCK(40,1)' -stim_label 1 happy \
        -stim_times 2 '1D: 0 120 240 360 480' 'BLOCK(40,1)' -stim_label 2 rest \
        -stim_times 3 '1D: 80 200 320 440' 'BLOCK(40,1)' -stim_label 3 count \
        -censor "censor_file.1D" \
        -ortvec "dfile_rall.1D" motion \
        -polort 0 \
        -fout -rout -tout \
        -full_first \
        -bucket $fdir/`echo $functional| cut -c 10-14`_standard-happy_canonical.ssreg \
        -GOFORIT 2
        cd ../..
    end
    cd ../..
end
