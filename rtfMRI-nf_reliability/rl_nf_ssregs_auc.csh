#!/bin/csh
# does single subject regressions with area under the curve modelling for the neurofeedback task
# as implemented on the dataset from Young et al (2017) AJP
#   Note that this paper combined data from two groups, and as such there is code below
#   for integrating across the 2 groups
# This script assumes we have the AFNI library in our path
# This code has been described in Compere et al. (2020)

# the script assumes we have regressors made which are 1 at trial onsets and 0 elsewhere for a given condition
# e.g., 
#  1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0....
# and also censor files with entries for every TR in which 0=censored. 1=uncensored
# Creates ROI and resamples it so it has the same size of the functionals
whereami -mask_atlas_region TT_Daemon:left:amygdala -prefix Lamygdala
3dresample -prefix Lamygdala_resampled -master data_for_voxel_wise_reliability_right_data/Active/AC820/Preprocessing_standard/Visit_1/AC820.results/all_runs.AC820+tlrc -inset Lamygdala+tlrc
set fdir = ssregs_auc
cd data_for_voxel_wise_reliability_right_data/Active
# BV style analysis
foreach sub(`ls`)
    cd $sub/Preprocessing_BV_style
    set functional = Preprocessed_all_runs.$sub+tlrc
    foreach visit(Visit_1 Visit_2)
        cd $visit/$sub.results
        3drefit -TR 2 $functional

        3dDeconvolve \
        -input $functional \
        -mask ../../../../../../Lamygdala_resampled+tlrc \
        -num_stimts 1 \
        -nfirst 1 \
        -stim_file 1  "../../../../../../stim_file_for_preprocessed_data.1D"   -stim_label 1 happy \
        -stim_minlag 1 0 \
        -stim_maxlag 1 19 \
        -iresp 1 $fdir/`echo $functional| cut -c 23-27`_BV_style-happy.irf \
        -polort 0 \
        -fout -rout -tout \
        -full_first \
        -bucket $fdir/`echo $functional| cut -c 23-27`_BV_style-happy.ssreg \
        -GOFORIT 2
        cd $fdir
        foreach irf ("`echo $functional| cut -c 23-27`_BV_style-happy")
            3dcalc -prefix ${irf}-AUC -expr "a+b+c+d+e+f+g+h+i+j+k+l+m+n+o+p+q+r+s+t" -a0 "${irf}.irf+tlrc" -b1 "${irf}.irf+tlrc" -c2 "${irf}.irf+tlrc" -d3 "${irf}.irf+tlrc" -e4 "${irf}.irf+tlrc" -f5 "${irf}.irf+tlrc" -g6 "${irf}.irf+tlrc" -h7 "${irf}.irf+tlrc" -i8 "${irf}.irf+tlrc" -j9 "${irf}.irf+tlrc" -k10 "${irf}.irf+tlrc" -l11 "${irf}.irf+tlrc" -m12 "${irf}.irf+tlrc" -n13 "${irf}.irf+tlrc" -o14 "${irf}.irf+tlrc" -p15 "${irf}.irf+tlrc" -q16 "${irf}.irf+tlrc" -r17 "${irf}.irf+tlrc" -s18 "${irf}.irf+tlrc" -t19 "${irf}.irf+tlrc"
        end

        # get the peak of the HRF
        foreach irf ("`echo $functional| cut -c 23-27`_BV_style-happy")
            3dcalc -prefix ${irf}-AMP -expr "max(a,max(b,max(c,max(d,max(e,max(f,max(g,max(h,max(i,max(j,max(k,max(l,max(m,max(n,max(o,max(p,max(q,max(r,max(s,t)))))))))))))))))))" -a0 "${irf}.irf+tlrc" -b1 "${irf}.irf+tlrc" -c2 "${irf}.irf+tlrc" -d3 "${irf}.irf+tlrc" -e4 "${irf}.irf+tlrc" -f5 "${irf}.irf+tlrc" -g6 "${irf}.irf+tlrc" -h7 "${irf}.irf+tlrc" -i8 "${irf}.irf+tlrc" -j9 "${irf}.irf+tlrc" -k10 "${irf}.irf+tlrc" -l11 "${irf}.irf+tlrc" -m12 "${irf}.irf+tlrc" -n13 "${irf}.irf+tlrc" -o14 "${irf}.irf+tlrc" -p15 "${irf}.irf+tlrc" -q16 "${irf}.irf+tlrc" -r17 "${irf}.irf+tlrc" -s18 "${irf}.irf+tlrc" -t19 "${irf}.irf+tlrc"
        end
        cd ../../..
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
        -stim_file 1  "../../../../../../stim_file_for_original_data_happy.1D"   -stim_label 1 happy \
        -stim_minlag 1 0 \
        -stim_maxlag 1 19 \
        -stim_file 2  "../../../../../../stim_file_for_original_data_rest.1D"   -stim_label 2 rest \
        -stim_minlag 2 0 \
        -stim_maxlag 2 19 \
        -stim_file 3  "../../../../../../stim_file_for_original_data_count.1D"   -stim_label 3 count \
        -stim_minlag 3 0 \
        -stim_maxlag 3 19 \
        -censor "censor_file.1D" \
        -ortvec "dfile_rall.1D" motion \
        -iresp 1 $fdir/`echo $functional| cut -c 10-14`_standard-happy.irf \
        -polort 0 \
        -fout -rout -tout \
        -full_first \
        -bucket $fdir/`echo $functional| cut -c 10-14`_standard-happy.ssreg \
        -GOFORIT 2

        cd $fdir

        foreach irf ("`echo $functional| cut -c 10-14`_standard-happy")
            3dcalc -prefix ${irf}-AUC -expr "a+b+c+d+e+f+g+h+i+j+k+l+m+n+o+p+q+r+s+t" -a0 "${irf}.irf+tlrc" -b1 "${irf}.irf+tlrc" -c2 "${irf}.irf+tlrc" -d3 "${irf}.irf+tlrc" -e4 "${irf}.irf+tlrc" -f5 "${irf}.irf+tlrc" -g6 "${irf}.irf+tlrc" -h7 "${irf}.irf+tlrc" -i8 "${irf}.irf+tlrc" -j9 "${irf}.irf+tlrc" -k10 "${irf}.irf+tlrc" -l11 "${irf}.irf+tlrc" -m12 "${irf}.irf+tlrc" -n13 "${irf}.irf+tlrc" -o14 "${irf}.irf+tlrc" -p15 "${irf}.irf+tlrc" -q16 "${irf}.irf+tlrc" -r17 "${irf}.irf+tlrc" -s18 "${irf}.irf+tlrc" -t19 "${irf}.irf+tlrc"
        end

       # get the peak of the HRF
        foreach irf ("`echo $functional| cut -c 10-14`_standard-happy")
           3dcalc -prefix ${irf}-AMP -expr "max(a,max(b,max(c,max(d,max(e,max(f,max(g,max(h,max(i,max(j,max(k,max(l,max(m,max(n,max(o,max(p,max(q,max(r,max(s,t)))))))))))))))))))" -a0 "${irf}.irf+tlrc" -b1 "${irf}.irf+tlrc" -c2 "${irf}.irf+tlrc" -d3 "${irf}.irf+tlrc" -e4 "${irf}.irf+tlrc" -f5 "${irf}.irf+tlrc" -g6 "${irf}.irf+tlrc" -h7 "${irf}.irf+tlrc" -i8 "${irf}.irf+tlrc" -j9 "${irf}.irf+tlrc" -k10 "${irf}.irf+tlrc" -l11 "${irf}.irf+tlrc" -m12 "${irf}.irf+tlrc" -n13 "${irf}.irf+tlrc" -o14 "${irf}.irf+tlrc" -p15 "${irf}.irf+tlrc" -q16 "${irf}.irf+tlrc" -r17 "${irf}.irf+tlrc" -s18 "${irf}.irf+tlrc" -t19 "${irf}.irf+tlrc"
        end

        cd ../../..
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
        3drefit -TR 2 $functional
        3dDeconvolve \
        -input $functional \
        -mask ../../../../../../Lamygdala_resampled+tlrc \
        -num_stimts 1 \
        -nfirst 1 \
        -stim_file 1  "../../../../../../stim_file_for_preprocessed_data.1D"   -stim_label 1 happy \
        -stim_minlag 1 0 \
        -stim_maxlag 1 19 \
        -iresp 1 $fdir/`echo $functional| cut -c 23-27`_BV_style-happy.irf \
        -polort 0 \
        -fout -rout -tout \
        -full_first \
        -bucket $fdir/`echo $functional| cut -c 23-27`_BV_style-happy.ssreg \
        -GOFORIT 2

        cd $fdir

        foreach irf ("`echo $functional| cut -c 23-27`_BV_style-happy")
            3dcalc -prefix ${irf}-AUC -expr "a+b+c+d+e+f+g+h+i+j+k+l+m+n+o+p+q+r+s+t" -a0 "${irf}.irf+tlrc" -b1 "${irf}.irf+tlrc" -c2 "${irf}.irf+tlrc" -d3 "${irf}.irf+tlrc" -e4 "${irf}.irf+tlrc" -f5 "${irf}.irf+tlrc" -g6 "${irf}.irf+tlrc" -h7 "${irf}.irf+tlrc" -i8 "${irf}.irf+tlrc" -j9 "${irf}.irf+tlrc" -k10 "${irf}.irf+tlrc" -l11 "${irf}.irf+tlrc" -m12 "${irf}.irf+tlrc" -n13 "${irf}.irf+tlrc" -o14 "${irf}.irf+tlrc" -p15 "${irf}.irf+tlrc" -q16 "${irf}.irf+tlrc" -r17 "${irf}.irf+tlrc" -s18 "${irf}.irf+tlrc" -t19 "${irf}.irf+tlrc"
        end

        # get the peak of the HRF
        foreach irf ("`echo $functional| cut -c 23-27`_BV_style-happy")
            3dcalc -prefix ${irf}-AMP -expr "max(a,max(b,max(c,max(d,max(e,max(f,max(g,max(h,max(i,max(j,max(k,max(l,max(m,max(n,max(o,max(p,max(q,max(r,max(s,t)))))))))))))))))))" -a0 "${irf}.irf+tlrc" -b1 "${irf}.irf+tlrc" -c2 "${irf}.irf+tlrc" -d3 "${irf}.irf+tlrc" -e4 "${irf}.irf+tlrc" -f5 "${irf}.irf+tlrc" -g6 "${irf}.irf+tlrc" -h7 "${irf}.irf+tlrc" -i8 "${irf}.irf+tlrc" -j9 "${irf}.irf+tlrc" -k10 "${irf}.irf+tlrc" -l11 "${irf}.irf+tlrc" -m12 "${irf}.irf+tlrc" -n13 "${irf}.irf+tlrc" -o14 "${irf}.irf+tlrc" -p15 "${irf}.irf+tlrc" -q16 "${irf}.irf+tlrc" -r17 "${irf}.irf+tlrc" -s18 "${irf}.irf+tlrc" -t19 "${irf}.irf+tlrc"
        end
        cd ../../..
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
        -stim_file 1  "../../../../../../stim_file_for_original_data_happy.1D"   -stim_label 1 happy \
        -stim_minlag 1 0 \
        -stim_maxlag 1 19 \
        -stim_file 2  "../../../../../../stim_file_for_original_data_rest.1D"   -stim_label 2 rest \
        -stim_minlag 2 0 \
        -stim_maxlag 2 19 \
        -stim_file 3  "../../../../../../stim_file_for_original_data_count.1D"   -stim_label 3 count \
        -stim_minlag 3 0 \
        -stim_maxlag 3 19 \
        -censor "censor_file.1D" \
        -ortvec "dfile_rall.1D" motion \
        -iresp 1 $fdir/`echo $functional| cut -c 10-14`_standard-happy.irf \
        -polort 0 \
        -fout -rout -tout \
        -full_first \
        -bucket $fdir/`echo $functional| cut -c 10-14`_standard-happy.ssreg \
        -GOFORIT 2

        cd $fdir

        foreach irf ("`echo $functional| cut -c 10-14`_standard-happy")
            3dcalc -prefix ${irf}-AUC -expr "a+b+c+d+e+f+g+h+i+j+k+l+m+n+o+p+q+r+s+t" -a0 "${irf}.irf+tlrc" -b1 "${irf}.irf+tlrc" -c2 "${irf}.irf+tlrc" -d3 "${irf}.irf+tlrc" -e4 "${irf}.irf+tlrc" -f5 "${irf}.irf+tlrc" -g6 "${irf}.irf+tlrc" -h7 "${irf}.irf+tlrc" -i8 "${irf}.irf+tlrc" -j9 "${irf}.irf+tlrc" -k10 "${irf}.irf+tlrc" -l11 "${irf}.irf+tlrc" -m12 "${irf}.irf+tlrc" -n13 "${irf}.irf+tlrc" -o14 "${irf}.irf+tlrc" -p15 "${irf}.irf+tlrc" -q16 "${irf}.irf+tlrc" -r17 "${irf}.irf+tlrc" -s18 "${irf}.irf+tlrc" -t19 "${irf}.irf+tlrc"
        end

        # get the peak of the HRF
        foreach irf ("`echo $functional| cut -c 10-14`_standard-happy")
            3dcalc -prefix ${irf}-AMP -expr "max(a,max(b,max(c,max(d,max(e,max(f,max(g,max(h,max(i,max(j,max(k,max(l,max(m,max(n,max(o,max(p,max(q,max(r,max(s,t)))))))))))))))))))" -a0 "${irf}.irf+tlrc" -b1 "${irf}.irf+tlrc" -c2 "${irf}.irf+tlrc" -d3 "${irf}.irf+tlrc" -e4 "${irf}.irf+tlrc" -f5 "${irf}.irf+tlrc" -g6 "${irf}.irf+tlrc" -h7 "${irf}.irf+tlrc" -i8 "${irf}.irf+tlrc" -j9 "${irf}.irf+tlrc" -k10 "${irf}.irf+tlrc" -l11 "${irf}.irf+tlrc" -m12 "${irf}.irf+tlrc" -n13 "${irf}.irf+tlrc" -o14 "${irf}.irf+tlrc" -p15 "${irf}.irf+tlrc" -q16 "${irf}.irf+tlrc" -r17 "${irf}.irf+tlrc" -s18 "${irf}.irf+tlrc" -t19 "${irf}.irf+tlrc"
        end
        cd ../../..
    end
    cd ../..
end
