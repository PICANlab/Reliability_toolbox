function rl_nf_compute_ICCs()
% rl_compute_ICCs.m will compute ICCs and semi partial correlations map with and without covariates at the second level and store them in a corr structure saved in the current directory
% assumes we have made a binary mask the size of the functional data
% and that we have already run rl_cbtsretssregs_auc.csh, rl_cbtsretssregs_canonical_amplitude.csh, rl_compute_gamma_variates.m
% This script assumes we have the AFNI matlab library and others function of the reliability toolbox in our path
% This code has been described in Compere et al. (2020)

dsmask = 'Lamygdala_resampled+tlrc.BRIK';
% load the mask
[mask]=BrikLoad(dsmask);

data_master = struct('first_level_modelization',{'canonical_amplitude',...
                                                 'amplitude',...
                                                 'area_under_the_curve',...
                                                 'onset_delay',...
                                                 'rise_decay_rate',...
                                                 'height'},...
                     'look_for_active_data_BV_style',{'data_for_voxel_wise_reliability_right_data/Active/*/Preprocessing_BV_style/Visit_1/*.results/ssregs_canonical_amplitude/*_BV_style-happy_canonical.ssreg+tlrc.BRIK*',...
                                               'data_for_voxel_wise_reliability_right_data/Active/*/Preprocessing_BV_style/Visit_1/*.results/ssregs_auc/*_BV_style-happy-AMP+tlrc.BRIK*',...
                                               'data_for_voxel_wise_reliability_right_data/Active/*/Preprocessing_BV_style/Visit_1/*.results/ssregs_auc/*_BV_style-happy-AUC+tlrc.BRIK*',...
                                               'data_for_voxel_wise_reliability_right_data/Active/*/Preprocessing_BV_style/Visit_1/*.results/ssregs_gammavar_computed_with_matlab/gamfitb0_*_BV_style-happy+tlrc.BRIK',...
                                               'data_for_voxel_wise_reliability_right_data/Active/*/Preprocessing_BV_style/Visit_1/*.results/ssregs_gammavar_computed_with_matlab/gamfitb1_*_BV_style-happy+tlrc.BRIK',...
                                               'data_for_voxel_wise_reliability_right_data/Active/*/Preprocessing_BV_style/Visit_1/*.results/ssregs_gammavar_computed_with_matlab/gamfitb2_*_BV_style-happy+tlrc.BRIK'},...
                     'look_for_control_data_BV_style',{'data_for_voxel_wise_reliability_right_data/Control/*/Preprocessing_BV_style/Visit_1/*.results/ssregs_canonical_amplitude/*_BV_style-happy_canonical.ssreg+tlrc.BRIK*',...
                                                       'data_for_voxel_wise_reliability_right_data/Control/*/Preprocessing_BV_style/Visit_1/*.results/ssregs_auc/*_BV_style-happy-AMP+tlrc.BRIK*',...
                                                       'data_for_voxel_wise_reliability_right_data/Control/*/Preprocessing_BV_style/Visit_1/*.results/ssregs_auc/*_BV_style-happy-AUC+tlrc.BRIK*',...
                                                       'data_for_voxel_wise_reliability_right_data/Control/*/Preprocessing_BV_style/Visit_1/*.results/ssregs_gammavar_computed_with_matlab/gamfitb0_*_BV_style-happy+tlrc.BRIK',...
                                                       'data_for_voxel_wise_reliability_right_data/Control/*/Preprocessing_BV_style/Visit_1/*.results/ssregs_gammavar_computed_with_matlab/gamfitb1_*_BV_style-happy+tlrc.BRIK',...
                                                       'data_for_voxel_wise_reliability_right_data/Control/*/Preprocessing_BV_style/Visit_1/*.results/ssregs_gammavar_computed_with_matlab/gamfitb2_*_BV_style-happy+tlrc.BRIK'},...
                     'look_for_active_data_standard',{'data_for_voxel_wise_reliability_right_data/Active/*/Preprocessing_standard/Visit_1/*.results/ssregs_canonical_amplitude/*_standard-happy_canonical.ssreg+tlrc.BRIK*',...
                                               'data_for_voxel_wise_reliability_right_data/Active/*/Preprocessing_standard/Visit_1/*.results/ssregs_auc/*_standard-happy-AMP+tlrc.BRIK*',...
                                               'data_for_voxel_wise_reliability_right_data/Active/*/Preprocessing_standard/Visit_1/*.results/ssregs_auc/*_standard-happy-AUC+tlrc.BRIK*',...
                                               'data_for_voxel_wise_reliability_right_data/Active/*/Preprocessing_standard/Visit_1/*.results/ssregs_gammavar_computed_with_matlab/gamfitb0_*_standard-happy+tlrc.BRIK',...
                                               'data_for_voxel_wise_reliability_right_data/Active/*/Preprocessing_standard/Visit_1/*.results/ssregs_gammavar_computed_with_matlab/gamfitb1_*_standard-happy+tlrc.BRIK',...
                                               'data_for_voxel_wise_reliability_right_data/Active/*/Preprocessing_standard/Visit_1/*.results/ssregs_gammavar_computed_with_matlab/gamfitb2_*_standard-happy+tlrc.BRIK'},...
                     'look_for_control_data_standard',{'data_for_voxel_wise_reliability_right_data/Control/*/Preprocessing_standard/Visit_1/*.results/ssregs_canonical_amplitude/*_standard-happy_canonical.ssreg+tlrc.BRIK*',...
                                                       'data_for_voxel_wise_reliability_right_data/Control/*/Preprocessing_standard/Visit_1/*.results/ssregs_auc/*_standard-happy-AMP+tlrc.BRIK*',...
                                                       'data_for_voxel_wise_reliability_right_data/Control/*/Preprocessing_standard/Visit_1/*.results/ssregs_auc/*_standard-happy-AUC+tlrc.BRIK*',...
                                                       'data_for_voxel_wise_reliability_right_data/Control/*/Preprocessing_standard/Visit_1/*.results/ssregs_gammavar_computed_with_matlab/gamfitb0_*_standard-happy+tlrc.BRIK',...
                                                       'data_for_voxel_wise_reliability_right_data/Control/*/Preprocessing_standard/Visit_1/*.results/ssregs_gammavar_computed_with_matlab/gamfitb1_*_standard-happy+tlrc.BRIK',...
                                                       'data_for_voxel_wise_reliability_right_data/Control/*/Preprocessing_standard/Visit_1/*.results/ssregs_gammavar_computed_with_matlab/gamfitb2_*_standard-happy+tlrc.BRIK'},...
                     'bucket',{2,...
                               1,...
                               1,...
                               1,...
                               1,...
                               1},...
                     'filename',{'corr_canonical_amplitude',...
                                 'corr_amplitude',...
                                 'corr_area_under_the_curve',...
                                 'corr_onset_delay',...
                                 'corr_rise_decay_rate',...
                                 'corr_height'},...
                     'threshold',{-200,...
                                  0,...
                                  -37000,...
                                  -1000,...
                                  -100000000000000,...
                                  -1000});

[dat_num, dat_txt, dat_raw]=xlsread('K99_covariates.xlsx');

% get paths for data of interest
for fl=1:numel(data_master)
  clear pre post predat postdat pre_temp post_temp corr subnum ID covariates;
  active_data_BV_style=dir(data_master(fl).look_for_active_data_BV_style);
  control_data_BV_style=dir(data_master(fl).look_for_control_data_BV_style);
  active_data_standard=dir(data_master(fl).look_for_active_data_standard);
  control_data_standard=dir(data_master(fl).look_for_control_data_standard);
  
  % On data from the experimental group preprocessed with the BV style
  % pipeline
  % Create a matrix in which the order of the participants' data and their
 % covariates match
  subnum=1;
  for subj = 1:numel(active_data_BV_style)
      folder_post=active_data_BV_style(subj).folder;
      folder_post(133)='2'; 
      ID=folder_post(98:102); 
      fprintf('%d Processing %s\n',subj,active_data_BV_style(subj).name);
      pre_temp(:,:,:,:) = BrikLoad(fullfile(active_data_BV_style(subj).folder,active_data_BV_style(subj).name));
      pre(:,:,:,subnum) = pre_temp(:,:,:,data_master(fl).bucket).*mask;
      post_temp(:,:,:,:) = BrikLoad(fullfile(folder_post,active_data_BV_style(subj).name));
      post(:,:,:,subnum) = post_temp(:,:,:,data_master(fl).bucket).*mask;
      covariates(subnum,:)=dat_num(find(strcmp(ID,dat_txt(:,2)))-1,[1 2 3 4 5 6 8]); % taking off the 7 column only for participants in the experimental group but taking off the 9th column for all participants 
      subnum=subnum+1;
  end
  
  % In each voxels verify the quality of the data by comparing their values
  % to a threshold and making sure the values are different to 0, rescale
  % outliers and then compute icc and semi partial correlations with and
  % without covariates when enough data are available and get rid of all
  % weird outcomes from semipartial correlation computation
  [x, y, z, n]=size(pre);
  for dimx = 1:x
    for dimy = 1:y
      for dimz = 1:z
    predat=squeeze(pre(dimx,dimy,dimz,:));
	postdat=squeeze(post(dimx,dimy,dimz,:));
	gooddata=find((predat>data_master(fl).threshold) & (postdat>data_master(fl).threshold) & (abs(predat)>.001) & (abs(postdat)>.001));
    if length(gooddata)>3 
      predatuse=rescaleoutliers(predat(gooddata));
      postdatuse= rescaleoutliers(postdat(gooddata));
      corr.icc(dimx,dimy,dimz) = icc([predatuse postdatuse]);
      [mr]=mregs(rescaleoutliers(predat(gooddata)),rescaleoutliers(postdat(gooddata)));
      corr.s_p_r_without_covariates(dimx,dimy,dimz)=sqrt(mr.Rsq);
    else
       corr.icc(dimx,dimy,dimz)=nan;
       corr.s_p_r_without_covariates(dimx,dimy,dimz)=0;
    end
    if length(gooddata)>max(3,2+size(covariates,2))
        [mr]=mhreg(covariates(gooddata), rescaleoutliers(predat(gooddata)), rescaleoutliers(postdat(gooddata)));
        corr.s_p_r_with_covariates(dimx,dimy,dimz)=sqrt(mr.dr.Rsq);
    else
         corr.s_p_r_with_covariates(dimx,dimy,dimz)=0;
    end
    ind=corr.s_p_r_without_covariates(dimx,dimy,dimz);
	if (ind>1) | (isnan(ind)) | (ind<0) | (ind~=real(ind))
	  corr.s_p_r_without_covariates(dimx,dimy,dimz)=0;
	end
	ind=corr.s_p_r_with_covariates(dimx,dimy,dimz);
	if (ind>1) | (isnan(ind)) | (ind<0) | (ind~=real(ind))
	  corr.s_p_r_with_covariates(dimx,dimy,dimz)=0;
	end
      end
    end
  end
 
  
  save(sprintf('%s_active_data_BV_style',data_master(fl).filename),'corr')
  
  % Does exactly the same on the data from the experimental group
  % preprocessed with the BV style pipeline
  clear pre post predat postdat pre_temp post_temp corr subnum ID covariates;
  subnum=1;
  for subj = 1:numel(control_data_BV_style)
    folder_post=control_data_BV_style(subj).folder;
    folder_post(134)='2';  
    ID=folder_post(99:103); 
    fprintf('%d Processing %s\n',subj,control_data_BV_style(subj).name);
    pre_temp(:,:,:,:) = BrikLoad(fullfile(control_data_BV_style(subj).folder,control_data_BV_style(subj).name));
    pre(:,:,:,subnum) = pre_temp(:,:,:,data_master(fl).bucket).*mask;
    post_temp(:,:,:,:) = BrikLoad(fullfile(folder_post,control_data_BV_style(subj).name));
    post(:,:,:,subnum) = post_temp(:,:,:,data_master(fl).bucket).*mask;
    covariates(subnum,:)=dat_num(find(strcmp(ID,dat_txt(:,2)))-1,[1 2 3 4 5 6 7 8]); % taking off the 7 column only for participants in the experimental group but taking off the 9th column for all participants (sex)    
    subnum=subnum+1;
  end

  [x, y, z, n]=size(pre);
  for dimx = 1:x
      for dimy = 1:y
          for dimz = 1:z
              predat=squeeze(pre(dimx,dimy,dimz,:));
              postdat=squeeze(post(dimx,dimy,dimz,:));
              gooddata=find((predat>data_master(fl).threshold) & (postdat>data_master(fl).threshold) & (abs(predat)>.001) & (abs(postdat)>.001));
              if length(gooddata)>3
                  predatuse=rescaleoutliers(predat(gooddata));
                  postdatuse= rescaleoutliers(postdat(gooddata));
                  corr.icc(dimx,dimy,dimz) = icc([predatuse postdatuse]);
                  [mr]=mregs(rescaleoutliers(predat(gooddata)),rescaleoutliers(postdat(gooddata)));
                  corr.s_p_r_without_covariates(dimx,dimy,dimz)=sqrt(mr.Rsq);
              else
                  corr.icc(dimx,dimy,dimz)=nan;
                  corr.s_p_r_without_covariates(dimx,dimy,dimz)=0;
              end
              if length(gooddata)>max(3,2+size(covariates,2))
                  [mr]=mhreg(covariates(gooddata), rescaleoutliers(predat(gooddata)), rescaleoutliers(postdat(gooddata)));
                  corr.s_p_r_with_covariates(dimx,dimy,dimz)=sqrt(mr.dr.Rsq);
              else
                  corr.s_p_r_with_covariates(dimx,dimy,dimz)=0;
              end
              ind=corr.s_p_r_without_covariates(dimx,dimy,dimz);
              if (ind>1) | (isnan(ind)) | (ind<0) | (ind~=real(ind))
                  corr.s_p_r_without_covariates(dimx,dimy,dimz)=0;
              end
              ind=corr.s_p_r_with_covariates(dimx,dimy,dimz);
              if (ind>1) | (isnan(ind)) | (ind<0) | (ind~=real(ind))
                  corr.s_p_r_with_covariates(dimx,dimy,dimz)=0;
              end
          end
      end
  end
  
  save(sprintf('%s_control_data_BV_style',data_master(fl).filename),'corr')
  
  % Does exactly the same on the data from the experimental group
  % preprocessed with the standard pipeline
  clear pre post predat postdat pre_temp post_temp corr subnum ID covariates;
  subnum=1;
  for subj = 1:numel(active_data_standard)
      folder_post=active_data_standard(subj).folder;
      folder_post(133)='2'; 
      ID=folder_post(98:102); 
      fprintf('%d Processing %s\n',subj,active_data_standard(subj).name);
      pre_temp(:,:,:,:) = BrikLoad(fullfile(active_data_standard(subj).folder,active_data_standard(subj).name));
      pre(:,:,:,subnum) = pre_temp(:,:,:,data_master(fl).bucket).*mask;
      post_temp(:,:,:,:) = BrikLoad(fullfile(folder_post,active_data_standard(subj).name));
      post(:,:,:,subnum) = post_temp(:,:,:,data_master(fl).bucket).*mask;
      covariates(subnum,:)=dat_num(find(strcmp(ID,dat_txt(:,2)))-1,[1 2 3 4 5 6 8]); % taking off the 7 column only for participants in the experimental group but taking off the 9th column for all participants (sex) => Rajouter un if
      subnum=subnum+1;
  end
  [x, y, z, n]=size(pre);
  for dimx = 1:x
      for dimy = 1:y
          for dimz = 1:z
              predat=squeeze(pre(dimx,dimy,dimz,:));
              postdat=squeeze(post(dimx,dimy,dimz,:));
              gooddata=find((predat>data_master(fl).threshold) & (postdat>data_master(fl).threshold) & (abs(predat)>.001) & (abs(postdat)>.001));
              if length(gooddata)>3
                  predatuse=rescaleoutliers(predat(gooddata));
                  postdatuse= rescaleoutliers(postdat(gooddata));
                  corr.icc(dimx,dimy,dimz) = icc([predatuse postdatuse]);
                  [mr]=mregs(rescaleoutliers(predat(gooddata)),rescaleoutliers(postdat(gooddata)));
                  corr.s_p_r_without_covariates(dimx,dimy,dimz)=sqrt(mr.Rsq);
              else
                  corr.icc(dimx,dimy,dimz)=nan;
                  corr.s_p_r_without_covariates(dimx,dimy,dimz)=0;
              end
              if length(gooddata)>max(3,2+size(covariates,2))
                  [mr]=mhreg(covariates(gooddata), rescaleoutliers(predat(gooddata)), rescaleoutliers(postdat(gooddata)));
                  corr.s_p_r_with_covariates(dimx,dimy,dimz)=sqrt(mr.dr.Rsq);
              else
                  corr.s_p_r_with_covariates(dimx,dimy,dimz)=0;
              end
              ind=corr.s_p_r_without_covariates(dimx,dimy,dimz);
              if (ind>1) | (isnan(ind)) | (ind<0) | (ind~=real(ind))
                  corr.s_p_r_without_covariates(dimx,dimy,dimz)=0;
              end
              ind=corr.s_p_r_with_covariates(dimx,dimy,dimz);
              if (ind>1) | (isnan(ind)) | (ind<0) | (ind~=real(ind))
                  corr.s_p_r_with_covariates(dimx,dimy,dimz)=0;
              end
          end
      end
  end
  
  
  save(sprintf('%s_active_data_standard',data_master(fl).filename),'corr')
  
  % Does exactly the same on the data from the control group
  % preprocessed with the standard pipeline
  clear pre post predat postdat pre_temp post_temp corr subnum ID covariates;
  subnum=1;
  for subj = 1:numel(control_data_standard)
    folder_post=control_data_standard(subj).folder;
    folder_post(134)='2';  
    ID=folder_post(99:103); 
    fprintf('%d Processing %s\n',subj,control_data_standard(subj).name);
    pre_temp(:,:,:,:) = BrikLoad(fullfile(control_data_standard(subj).folder,control_data_standard(subj).name));
    pre(:,:,:,subnum) = pre_temp(:,:,:,data_master(fl).bucket).*mask;
    post_temp(:,:,:,:) = BrikLoad(fullfile(folder_post,control_data_standard(subj).name));
    post(:,:,:,subnum) = post_temp(:,:,:,data_master(fl).bucket).*mask;
    covariates(subnum,:)=dat_num(find(strcmp(ID,dat_txt(:,2)))-1,[1 2 3 4 5 6 7 8]); % taking off the 7 column only for participants in the experimental group but taking off the 9th column for all participants (sex)    
    subnum=subnum+1;
  end

  [x, y, z, n]=size(pre);
  for dimx = 1:x
      for dimy = 1:y
          for dimz = 1:z
              predat=squeeze(pre(dimx,dimy,dimz,:));
              postdat=squeeze(post(dimx,dimy,dimz,:));
              gooddata=find((predat>data_master(fl).threshold) & (postdat>data_master(fl).threshold) & (abs(predat)>.001) & (abs(postdat)>.001));
              if length(gooddata)>3
                  predatuse=rescaleoutliers(predat(gooddata));
                  postdatuse= rescaleoutliers(postdat(gooddata));
                  corr.icc(dimx,dimy,dimz) = icc([predatuse postdatuse]);
                  [mr]=mregs(rescaleoutliers(predat(gooddata)),rescaleoutliers(postdat(gooddata)));
                  corr.s_p_r_without_covariates(dimx,dimy,dimz)=sqrt(mr.Rsq);
              else
                  corr.icc(dimx,dimy,dimz)=nan;
                  corr.s_p_r_without_covariates(dimx,dimy,dimz)=0;
              end
              if length(gooddata)>max(3,2+size(covariates,2))
                  [mr]=mhreg(covariates(gooddata), rescaleoutliers(predat(gooddata)), rescaleoutliers(postdat(gooddata)));
                  corr.s_p_r_with_covariates(dimx,dimy,dimz)=sqrt(mr.dr.Rsq);
              else
                  corr.s_p_r_with_covariates(dimx,dimy,dimz)=0;
              end
              ind=corr.s_p_r_without_covariates(dimx,dimy,dimz);
              if (ind>1) | (isnan(ind)) | (ind<0) | (ind~=real(ind))
                  corr.s_p_r_without_covariates(dimx,dimy,dimz)=0;
              end
              ind=corr.s_p_r_with_covariates(dimx,dimy,dimz);
              if (ind>1) | (isnan(ind)) | (ind<0) | (ind~=real(ind))
                  corr.s_p_r_with_covariates(dimx,dimy,dimz)=0;
              end
          end
      end
  end
  
  save(sprintf('%s_control_data_standard',data_master(fl).filename),'corr')
end

