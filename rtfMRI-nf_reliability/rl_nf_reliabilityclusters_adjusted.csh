#!/bin/csh
foreach reliability ('0.40' '0.60' '0.70' '0.75')
 mkdir -p cluster_correction_adjusted/$reliability/ICC_maps_masked_by_ROI_and_threshold
 mkdir -p cluster_correction_adjusted/$reliability/masks_cluster_correction_by_map_and_ROI

 # for 18 patients in experimental group r=.4 yields p=.100025
 # for 18 patients in experimental group r=.6 yields p=.008479
 # for 18 patients in experimental group r=.7 yields p=.001219
 # for 18 patients in experimental group r=.75 yields p=.000338
 # for 16 patients in control group r=.4 yields p=.12475
 # for 16 patients in control group r=.6 yields p=.014007
 # for 16 patients in control group r=.7 yields p=.002535
 # for 16 patients in control group r=.75 yields p=.00082

 if ($reliability == 0.40) then
	 set threshold_expe='0.100025'
	 set threshold_ctrl='0.12475'
 else if ($reliability == 0.60) then
	 set threshold_expe='0.008479'
	 set threshold_ctrl='0.014007'
 else if ($reliability == 0.70) then
	 set threshold_expe='0.001219'
	 set threshold_ctrl='0.002535'
 else if ($reliability == 0.75) then
	 set threshold_expe='0.000338'
	 set threshold_ctrl='0.00082'
 else
	 echo 'Issue with reliability threshold'
 endif
  foreach meas ('amplitude' 'area_under_the_curve' 'canonical_amplitude' 'height' 'onset_delay' 'rise_decay_rate')
  #1
  set icc_map=icc_${meas}_active_data_BV_style_computed_with_matlab+tlrc
  3dcalc -prefix resid_icc_${meas}_active_data_BV_style_computed_with_matlab -expr "1-a" -a $icc_map
  set mask=Lamygdala
  touch 3dFWHMx.1D
  rm 3dFWHMx.1D
  set effect = resid_icc_${meas}_active_data_BV_style_computed_with_matlab
  3dFWHMx -acf -dset "${effect}+tlrc"  -mask ${mask}_resampled+tlrc  | grep -v "+" | tail -n 1 > ${effect}_${mask}_fwhmx_acf.txt 
  set athresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 1`
  set bthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 2`
  set cthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 3`
  set estthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 4`
  # radius ACF-r model-r gaussian_NEWmodel-r-r

  # get p-values based on our correlation thresholds
  3dClustSim -mask "${mask}_resampled+tlrc" -acf ${athresh} ${bthresh} ${cthresh} -pthr $threshold_expe  -athr .05  -iter 10000 -prefix cluster_correction_adjusted/$reliability/clustsim_${effect}_${mask}

  set numvox = `tail -1 cluster_correction_adjusted/$reliability/clustsim_${effect}_${mask}.NN3_2sided.1D | colex 2`
  echo USING $numvox VOXELS FOR CONTIGUITY THRESHOLD FOR $effect $mask at $reliability >> cluster_correction_adjusted/clustmaskthresholds.txt

  # Create dataset with only ICC values in ROI at some ICC threshold (thresholding)
  3dcalc -a $icc_map -b ${mask}_resampled+tlrc -expr "step(b)*step(a-$reliability)" -prefix cluster_correction_adjusted/$reliability/ICC_maps_masked_by_ROI_and_threshold/$mask/$icc_map
  # Create mask for cluster correction to apply to each map and thresholding
  # rmm: sqrt(1.75^2+1.755^2+1.75^2)=3.0340
  set vmul = `awk "BEGIN {print 1.75*1.75*1.75*$numvox}"`
  mkdir -p cluster_correction_adjusted/$reliability/masks_cluster_correction_by_map_and_ROI/${meas}/active/BV_style
  3dclust -savemask cluster_correction_adjusted/$reliability/masks_cluster_correction_by_map_and_ROI/${meas}/active/BV_style/cluster_correction_mask  3.034 $vmul   cluster_correction_adjusted/$reliability/ICC_maps_masked_by_ROI_and_threshold/$mask/$icc_map

#2
 set icc_map=icc_${meas}_active_data_standard_computed_with_matlab+tlrc
 3dcalc -prefix resid_icc_${meas}_active_data_standard_computed_with_matlab -expr "1-a" -a $icc_map
 touch 3dFWHMx.1D
 rm 3dFWHMx.1D
 set effect = resid_icc_${meas}_active_data_standard_computed_with_matlab
 3dFWHMx -acf -dset "${effect}+tlrc"  -mask ${mask}_resampled+tlrc  | grep -v "+" | tail -n 1 > ${effect}_${mask}_fwhmx_acf.txt 
 set athresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 1`
 set bthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 2`
 set cthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 3`
 set estthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 4`
 # radius ACF-r model-r gaussian_NEWmodel-r-r

  # get p-values based on our correlation thresholds
  3dClustSim -mask "${mask}_resampled+tlrc" -acf ${athresh} ${bthresh} ${cthresh} -pthr $threshold_expe  -athr .05  -iter 10000 -prefix cluster_correction_adjusted/$reliability/clustsim_${effect}_${mask}

  set numvox = `tail -1 cluster_correction_adjusted/$reliability/clustsim_${effect}_${mask}.NN3_2sided.1D | colex 2`
  echo USING $numvox VOXELS FOR CONTIGUITY THRESHOLD FOR $effect $mask at $reliability >> cluster_correction_adjusted/clustmaskthresholds.txt

  # Create dataset with only ICC values in ROI at some ICC threshold (thresholding)
  3dcalc -a $icc_map -b ${mask}_resampled+tlrc -expr "step(b)*step(a-$reliability)" -prefix cluster_correction_adjusted/$reliability/ICC_maps_masked_by_ROI_and_threshold/$mask/$icc_map
  # Create mask for cluster correction to apply to each map and thresholding
  # rmm: sqrt(1.75^2+1.755^2+1.75^2)=3.0340
  set vmul = `awk "BEGIN {print 1.75*1.75*1.75*$numvox}"`
  mkdir -p cluster_correction_adjusted/$reliability/masks_cluster_correction_by_map_and_ROI/${meas}/active/standard
  3dclust -savemask cluster_correction_adjusted/$reliability/masks_cluster_correction_by_map_and_ROI/${meas}/active/standard/cluster_correction_mask  3.034 $vmul   cluster_correction_adjusted/$reliability/ICC_maps_masked_by_ROI_and_threshold/$mask/$icc_map

  #3
  set icc_map=icc_${meas}_control_data_BV_style_computed_with_matlab+tlrc
  3dcalc -prefix resid_icc_${meas}_control_data_BV_style_computed_with_matlab -expr "1-a" -a $icc_map
  touch 3dFWHMx.1D
  rm 3dFWHMx.1D
  set effect = resid_icc_${meas}_control_data_BV_style_computed_with_matlab
  3dFWHMx -acf -dset "${effect}+tlrc"  -mask ${mask}_resampled+tlrc  | grep -v "+" | tail -n 1 > ${effect}_${mask}_fwhmx_acf.txt 
  set athresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 1`
  set bthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 2`
  set cthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 3`
  set estthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 4`
  # radius ACF-r model-r gaussian_NEWmodel-r-r

 # get p-values based on our correlation thresholds
  3dClustSim -mask "${mask}_resampled+tlrc" -acf ${athresh} ${bthresh} ${cthresh} -pthr $threshold_ctrl  -athr .05  -iter 10000 -prefix cluster_correction_adjusted/$reliability/clustsim_${effect}_${mask}

  set numvox = `tail -1 cluster_correction_adjusted/$reliability/clustsim_${effect}_${mask}.NN3_2sided.1D | colex 2`
  echo USING $numvox VOXELS FOR CONTIGUITY THRESHOLD FOR $effect $mask at $reliability >> cluster_correction_adjusted/clustmaskthresholds.txt

  # Create dataset with only ICC values in ROI at some ICC threshold (thresholding)
  3dcalc -a $icc_map -b ${mask}_resampled+tlrc -expr "step(b)*step(a-$reliability)" -prefix cluster_correction_adjusted/$reliability/ICC_maps_masked_by_ROI_and_threshold/$mask/$icc_map
  # Create mask for cluster correction to apply to each map and thresholding
  # rmm: sqrt(1.75^2+1.755^2+1.75^2)=3.0340
  set vmul = `awk "BEGIN {print 1.75*1.75*1.75*$numvox}"`
  mkdir -p cluster_correction_adjusted/$reliability/masks_cluster_correction_by_map_and_ROI/${meas}/control/BV_style
  3dclust -savemask cluster_correction_adjusted/$reliability/masks_cluster_correction_by_map_and_ROI/${meas}/control/BV_style/cluster_correction_mask  3.034 $vmul   cluster_correction_adjusted/$reliability/ICC_maps_masked_by_ROI_and_threshold/$mask/$icc_map

  #4
  set icc_map=icc_${meas}_control_data_standard_computed_with_matlab+tlrc
  3dcalc -prefix resid_icc_${meas}_control_data_standard_computed_with_matlab -expr "1-a" -a $icc_map
  touch 3dFWHMx.1D
  rm 3dFWHMx.1D
  set effect = resid_icc_${meas}_control_data_standard_computed_with_matlab
  3dFWHMx -acf -dset "${effect}+tlrc"  -mask ${mask}_resampled+tlrc  | grep -v "+" | tail -n 1 > ${effect}_${mask}_fwhmx_acf.txt 
  set athresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 1`
  set bthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 2`
  set cthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 3`
  set estthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 4`
  # radius ACF-r model-r gaussian_NEWmodel-r-r

  # get p-values based on our correlation thresholds
  3dClustSim -mask "${mask}_resampled+tlrc" -acf ${athresh} ${bthresh} ${cthresh} -pthr $threshold_ctrl  -athr .05  -iter 10000 -prefix cluster_correction_adjusted/$reliability/clustsim_${effect}_${mask}

  set numvox = `tail -1 cluster_correction_adjusted/$reliability/clustsim_${effect}_${mask}.NN3_2sided.1D | colex 2`
  echo USING $numvox VOXELS FOR CONTIGUITY THRESHOLD FOR $effect $mask at $reliability >> cluster_correction_adjusted/clustmaskthresholds.txt

  # Create dataset with only ICC values in ROI at some ICC threshold (thresholding)
  3dcalc -a $icc_map -b ${mask}_resampled+tlrc -expr "step(b)*step(a-$reliability)" -prefix cluster_correction_adjusted/$reliability/ICC_maps_masked_by_ROI_and_threshold/$mask/$icc_map
  # Create mask for cluster correction to apply to each map and thresholding
  # rmm: sqrt(1.75^2+1.755^2+1.75^2)=3.0340
  set vmul = `awk "BEGIN {print 1.75*1.75*1.75*$numvox}"`
  mkdir -p cluster_correction_adjusted/$reliability/masks_cluster_correction_by_map_and_ROI/${meas}/control/standard
  3dclust -savemask cluster_correction_adjusted/$reliability/masks_cluster_correction_by_map_and_ROI/${meas}/control/standard/cluster_correction_mask  3.034 $vmul   cluster_correction_adjusted/$reliability/ICC_maps_masked_by_ROI_and_threshold/$mask/$icc_map
 
 end
end
