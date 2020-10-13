function rl_correlation_helper(fl)
% helper routine called by rl_compute_ICCs_eliminate_outliers() that
% compute reliability estimates on whole sample
% This script assumes we have the AFNI matlab library and others function of the reliability toolbox in our path
% This code has been described in Compere et al. (2020)

global data_master dsmask mask dat dat_pat;

% get paths for data of interest
  clear pre post predat postdat pre_temp post_temp corr subnum ID filter_covariate covariates;
  fprintf('Working with fl=%d\n',fl);
  controls=dir(data_master(fl).look_for_controls_data);
  patients_cohort1=dir(data_master(fl).look_for_patients_cohort1_data);
  patients_cohort2=dir(data_master(fl).look_for_patients_cohort2_data);
  all=[controls; patients_cohort1; patients_cohort2];
  patients=[patients_cohort1; patients_cohort2];
  
  
  
  subnum=1;
  for subj = 1:numel(all)
      name_post=all(subj).name;
      if fl<4
          name_post(5)=name_post(5)+2;
          ID=str2num(all(subj).name(5:8));
      else
          name_post(14)=name_post(14)+2;
          ID=str2num(all(subj).name(14:17));
      end
    fprintf('%d Looking for %s and %s\n',subj,all(subj).name, name_post);

 % Create a matrix in which the order of the participants' data and their
 % covariates match and a filter for when covariates are available or not
 % and a filter to do computations of reliability with covariates when
 % covariates are available
    if exist(fullfile(all(subj).folder,name_post))
      %fprintf('   Found\n');
      pre_temp(:,:,:,:) = BrikLoad(fullfile(all(subj).folder,all(subj).name));
      pre(:,:,:,subnum) = pre_temp(:,:,:,data_master(fl).bucket).*mask;
      post_temp(:,:,:,:) = BrikLoad(fullfile(all(subj).folder,name_post));
      post(:,:,:,subnum) = post_temp(:,:,:,data_master(fl).bucket).*mask;
      if length(unique(find(dat(:,1)==ID)))~=0 && length(unique(dat(find(dat(:,1)==ID),2:end)))~=1
        filter_covariate(subnum)=1;
        covariates(subnum,:)=dat(find(dat(:,1)==ID),:);
      else
        filter_covariate(subnum)=0;
        covariates(subnum,:)=[0 0 0 0 0 0 0 0 0];
      end
      subnum=subnum+1;
    else
      fprintf('  *** unidentified participant: %s\n',name_post);
     
    end
  end
  
  % In each voxels verify the quality of the data by comparing their values
  % to a threshold and making sure the values are different to 0, rescale
  % outliers and then compute icc and semi partial correlations with and
  % without covariates when enought data are available and get read of all
  % weird outcomes from semipartial correlation computation
  [x, y, z, n]=size(pre);
  for dimx = 1:x
    for dimy = 1:y
      for dimz = 1:z
	predat=squeeze(pre(dimx,dimy,dimz,:));
	postdat=squeeze(post(dimx,dimy,dimz,:));
	gooddata=find((predat>data_master(fl).threshold) & (postdat>data_master(fl).threshold) & (abs(predat)>.001) & (abs(postdat)>.001));
    curfilter_covariate=find((predat>data_master(fl).threshold) & (postdat>data_master(fl).threshold) & (abs(predat)>.001) & (abs(postdat)>.001)  & (filter_covariate'>0));
    if length(gooddata)>3 
      predatuse=rescaleoutliers(predat(gooddata));
      postdatuse= rescaleoutliers(postdat(gooddata));
      corr.icc(dimx,dimy,dimz) = icc([predatuse postdatuse]);
    else
       corr.icc(dimx,dimy,dimz)=nan;
    end
    if length(curfilter_covariate)>3
        [mr]=mregs(rescaleoutliers(predat(curfilter_covariate)),rescaleoutliers(postdat(curfilter_covariate)));
        corr.s_p_r_without_covariates(dimx,dimy,dimz)=sqrt(mr.Rsq);
    else
       corr.s_p_r_without_covariates(dimx,dimy,dimz)=0;
    end
	if length(curfilter_covariate)>max(3,2+size(covariates,2))
	  [mr]=mhreg(covariates(curfilter_covariate,2:end), rescaleoutliers(predat(curfilter_covariate)), rescaleoutliers(postdat(curfilter_covariate)));
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
  
  
  save(sprintf('%s_all_sample',data_master(fl).filename),'corr')
  
