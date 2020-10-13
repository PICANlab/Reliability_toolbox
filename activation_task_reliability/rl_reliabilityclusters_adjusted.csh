#!/bin/csh
foreach reliability ('0.40' '0.60' '0.70' '0.75')
 mkdir -p cluster_correction_adjusted/$reliability/ICC_maps_masked_by_ROI_and_threshold
 mkdir -p cluster_correction_adjusted/$reliability/masks_cluster_correction_by_map_and_ROI

 # for 91 participants in whole sample group r=.4 yields p=.100025
 # for 91 participants in whole sample group r=.6 yields p=.008479
 # for 91 participants in whole sample group r=.7 yields p=.001219
 # for 91 participants in whole sample group r=.75 yields p=.000338
 # for 57 patients in patients group r=.4 yields p=.12475
 # for 57 patients in patients group r=.6 yields p=.014007
 # for 57 patients in patients group r=.7 yields p=.002535
 # for 57 patients in patients group r=.75 yields p=.00082

 if ($reliability == 0.40) then
	 set threshold_all_sample='0.000086'
	 set threshold_patients='0.00205'
 else if ($reliability == 0.60) then
	 set threshold_all_sample='0.00001'
	 set threshold_patients='0.00001'
 else if ($reliability == 0.70) then
	 set threshold_all_sample='0.00001'
	 set threshold_patients='0.00001'
 else if ($reliability == 0.75) then
	 set threshold_all_sample='0.00001'
	 set threshold_patients='0.00001'
 else
	 echo 'Issue with reliability threshold'
 endif

foreach meas ('amplitude' 'area_under_the_curve' 'canonical_amplitude' 'height' 'onset_delay' 'rise_decay_rate')
# all_sample
 set sample = all_sample
 set icc_map=icc_${meas}_${sample}_computed_with_matlab+orig
 3dcalc -prefix resid_icc_${meas}_${sample}_computed_with_matlab -expr "1-a" -a $icc_map
# amygdala, DLPFC, rACC
 foreach mask ('Amygdala_MNIspace' 'DLPFC' 'rACC')
  touch 3dFWHMx.1D
  rm 3dFWHMx.1D
  set effect = resid_icc_${meas}_${sample}_computed_with_matlab
  3dFWHMx -acf -dset "${effect}+orig"  -mask binminicolin_${mask}_epan+orig  | grep -v "+" | tail -n 1 > ${effect}_${mask}_fwhmx_acf.txt 
  set athresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 1`
  set bthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 2`
  set cthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 3`
  set estthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 4`
  # radius ACF-r model-r gaussian_NEWmodel-r-r

 # get p-values based on our correlation thresholds
  3dClustSim -mask "binminicolin_${mask}_epan+orig" -acf ${athresh} ${bthresh} ${cthresh} -pthr $threshold_all_sample  -athr .05  -iter 10000 -prefix cluster_correction_adjusted/$reliability/clustsim_${effect}_${mask}

  set numvox = `tail -1 cluster_correction_adjusted/$reliability/clustsim_${effect}_${mask}.NN3_2sided.1D | colex 2`
  echo USING $numvox VOXELS FOR CONTIGUITY THRESHOLD FOR $effect $mask at $reliability >> cluster_correction_adjusted/clustmaskthresholds.txt

  # Create dataset with only ICC values in ROI at some ICC threshold (thresholding)
  3dcalc -a $icc_map -b "binminicolin_${mask}_epan+orig" -expr "step(b)*step(a-$reliability)" -prefix cluster_correction_adjusted/$reliability/ICC_maps_masked_by_ROI_and_threshold/$mask/$icc_map
  # Create mask for cluster correction to apply to each map and thresholding
  # rmm: sqrt(3.125^2+3.125^2+3.2^2)=5.4563
  # rmm: sqrt(1.75^2+1.755^2+1.75^2)=3.0340
  set vmul = `awk "BEGIN {print 3.125*3.125*3.2*$numvox}"`
  mkdir -p cluster_correction_adjusted/$reliability/masks_cluster_correction_by_map_and_ROI/${mask}/${meas}/${sample}
  3dclust -savemask cluster_correction_adjusted/$reliability/masks_cluster_correction_by_map_and_ROI/${mask}/${meas}/${sample}/cluster_correction_mask  5.4563 $vmul   cluster_correction_adjusted/$reliability/ICC_maps_masked_by_ROI_and_threshold/$mask/$icc_map
 end

# sgACC
 foreach mask ('sgACC_mask' 'sgACC_liberally_thresholded')
  touch 3dFWHMx.1D
  rm 3dFWHMx.1D
  set effect = resid_icc_${meas}_${sample}_computed_with_matlab
  3dFWHMx -acf -dset "${effect}+orig"  -mask ${mask}+orig  | grep -v "+" | tail -n 1 > ${effect}_${mask}_fwhmx_acf.txt 
  set athresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 1`
  set bthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 2`
  set cthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 3`
  set estthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 4`
  # radius ACF-r model-r gaussian_NEWmodel-r-r

 # get p-values based on our correlation thresholds
  3dClustSim -mask "${mask}+orig" -acf ${athresh} ${bthresh} ${cthresh} -pthr $threshold_all_sample  -athr .05  -iter 10000 -prefix cluster_correction_adjusted/$reliability/clustsim_${effect}_${mask}

  set numvox = `tail -1 cluster_correction_adjusted/$reliability/clustsim_${effect}_${mask}.NN3_2sided.1D | colex 2`
  echo USING $numvox VOXELS FOR CONTIGUITY THRESHOLD FOR $effect $mask at $reliability >> cluster_correction_adjusted/clustmaskthresholds.txt

  # Create dataset with only ICC values in ROI at some ICC threshold (thresholding)
  3dcalc -a $icc_map -b "${mask}_epan+orig" -expr "step(b)*step(a-$reliability)" -prefix cluster_correction_adjusted/$reliability/ICC_maps_masked_by_ROI_and_threshold/$mask/$icc_map
  # Create mask for cluster correction to apply to each map and thresholding
  # rmm: sqrt(3.125^2+3.125^2+3.2^2)=5.4563
  # rmm: sqrt(1.75^2+1.755^2+1.75^2)=3.0340
  set vmul = `awk "BEGIN {print 3.125*3.125*3.2*$numvox}"`
  mkdir -p cluster_correction_adjusted/$reliability/masks_cluster_correction_by_map_and_ROI/${mask}/${meas}/${sample}
  3dclust -savemask cluster_correction_adjusted/$reliability/masks_cluster_correction_by_map_and_ROI/${mask}/${meas}/${sample}/cluster_correction_mask  5.4563 $vmul   cluster_correction_adjusted/$reliability/ICC_maps_masked_by_ROI_and_threshold/$mask/$icc_map
 end


# patients
 set sample = patients
 set icc_map=icc_${meas}_${sample}_computed_with_matlab+orig
 3dcalc -prefix resid_icc_${meas}_${sample}_computed_with_matlab -expr "1-a" -a $icc_map
# amygdala, DLPFC, rACC
 foreach mask ('Amygdala_MNIspace' 'DLPFC' 'rACC') # add also the special ROI
  touch 3dFWHMx.1D
  rm 3dFWHMx.1D
  set effect = resid_icc_${meas}_${sample}_computed_with_matlab
  3dFWHMx -acf -dset "${effect}+orig"  -mask binminicolin_${mask}_epan+orig  | grep -v "+" | tail -n 1 > ${effect}_${mask}_fwhmx_acf.txt 
  set athresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 1`
  set bthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 2`
  set cthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 3`
  set estthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 4`
  # radius ACF-r model-r gaussian_NEWmodel-r-r

 # get p-values based on our correlation thresholds
  3dClustSim -mask "binminicolin_${mask}_epan+orig" -acf ${athresh} ${bthresh} ${cthresh} -pthr $threshold_patients  -athr .05  -iter 10000 -prefix cluster_correction_adjusted/$reliability/clustsim_${effect}_${mask}

  set numvox = `tail -1 cluster_correction_adjusted/$reliability/clustsim_${effect}_${mask}.NN3_2sided.1D | colex 2`
  echo USING $numvox VOXELS FOR CONTIGUITY THRESHOLD FOR $effect $mask at $reliability >> cluster_correction_adjusted/clustmaskthresholds.txt

  # Create dataset with only ICC values in ROI at some ICC threshold (thresholding)
  3dcalc -a $icc_map -b "binminicolin_${mask}_epan+orig" -expr "step(b)*step(a-$reliability)" -prefix cluster_correction_adjusted/$reliability/ICC_maps_masked_by_ROI_and_threshold/$mask/$icc_map
  # Create mask for cluster correction to apply to each map and thresholding
  # rmm: sqrt(3.125^2+3.125^2+3.2^2)=5.4563
  # rmm: sqrt(1.75^2+1.755^2+1.75^2)=3.0340
  set vmul = `awk "BEGIN {print 3.125*3.125*3.2*$numvox}"`
  mkdir -p cluster_correction_adjusted/$reliability/masks_cluster_correction_by_map_and_ROI/${mask}/${meas}/${sample}
  3dclust -savemask cluster_correction_adjusted/$reliability/masks_cluster_correction_by_map_and_ROI/$mask/${meas}/${sample}/cluster_correction_mask  5.4563 $vmul   cluster_correction_adjusted/$reliability/ICC_maps_masked_by_ROI_and_threshold/$mask/$icc_map
 end

# sgACC
 foreach mask ('sgACC_mask' 'sgACC_liberally_thresholded')
  touch 3dFWHMx.1D
  rm 3dFWHMx.1D
  set effect = resid_icc_${meas}_${sample}_computed_with_matlab
  3dFWHMx -acf -dset "${effect}+orig"  -mask ${mask}+orig  | grep -v "+" | tail -n 1 > ${effect}_${mask}_fwhmx_acf.txt 
  set athresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 1`
  set bthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 2`
  set cthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 3`
  set estthresh = `cat ${effect}_${mask}_fwhmx_acf.txt | colex 4`
  # radius ACF-r model-r gaussian_NEWmodel-r-r

 # get p-values based on our correlation thresholds
  3dClustSim -mask "${mask}+orig" -acf ${athresh} ${bthresh} ${cthresh} -pthr $threshold_patients  -athr .05  -iter 10000 -prefix cluster_correction_adjusted/$reliability/clustsim_${effect}_${mask}

  set numvox = `tail -1 cluster_correction_adjusted/$reliability/clustsim_${effect}_${mask}.NN3_2sided.1D | colex 2`
  echo USING $numvox VOXELS FOR CONTIGUITY THRESHOLD FOR $effect $mask at $reliability >> cluster_correction_adjusted/clustmaskthresholds.txt

  # Create dataset with only ICC values in ROI at some ICC threshold (thresholding)
  3dcalc -a $icc_map -b "${mask}_epan+orig" -expr "step(b)*step(a-$reliability)" -prefix cluster_correction_adjusted/$reliability/ICC_maps_masked_by_ROI_and_threshold/$mask/$icc_map
  # Create mask for cluster correction to apply to each map and thresholding
  # rmm: sqrt(3.125^2+3.125^2+3.2^2)=5.4563
  # rmm: sqrt(1.75^2+1.755^2+1.75^2)=3.0340
  set vmul = `awk "BEGIN {print 3.125*3.125*3.2*$numvox}"`
  mkdir -p cluster_correction_adjusted/$reliability/masks_cluster_correction_by_map_and_ROI/${mask}/${meas}/${sample}
  3dclust -savemask cluster_correction_adjusted/$reliability/masks_cluster_correction_by_map_and_ROI/${mask}/${meas}/${sample}/cluster_correction_mask  5.4563 $vmul   cluster_correction_adjusted/$reliability/ICC_maps_masked_by_ROI_and_threshold/$mask/$icc_map
 end
end
end
