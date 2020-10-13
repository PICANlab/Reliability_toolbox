function rl_nf_compute_gamma_variates()
% does single subject gamma variate fit voxelwise
% as implemented on the dataset from Young et al (2017) AJP
%   Note that this paper combined data from two groups, and as such there is code below
%   for integrating across the 2 groups
% the script assumes we have calculated IRFs via rl_nf_ssregs_auc.csh
% This script assumes we have the AFNI matlab library and others function of the reliability toolbox in our path
% This code has been described in Compere et al. (2020)

dsmask = '/media/youngklab/data/K99/reliability_K99_data/Lamygdala_resampled+tlrc';


data_master = struct('look_for_data',{'data_for_voxel_wise_reliability_right_data/*/*/*/Visit_*/*.results/ssregs_auc/*_*-happy.irf+tlrc.BRIK*'});
 
data=dir(data_master.look_for_data);
fdir='ssregs_gammavar_computed_with_matlab';

% Unzipping to get real name %
for irf=1:size(data,1)
    [data_to_rescale, data_to_rescale_info]=BrikLoad(fullfile(data(irf).folder,data(irf).name));
end

% We start by rewriting every irf file in this folder rescaling the data
% from the first subbrik

data=dir(data_master.look_for_data);

for irf=1:size(data,1)
    mkdir(fullfile(data(irf).folder,'..',fdir))
    cd(fullfile(data(irf).folder,'..',fdir))
    [data_to_rescale, data_to_rescale_info]=BrikLoad(fullfile(data(irf).folder,data(irf).name));
    for brik=1:size(data_to_rescale,4)
       data_rescaled(:,:,:,brik)=data_to_rescale(:,:,:,brik)-data_to_rescale(:,:,:,1);
    end
 
   data_to_rescale_info.RootName=strcat(data(irf).name(1:end-10),'_rescaled+tlrc');
   optOut.Scale = 1;
   optOut.Prefix = data_to_rescale_info.RootName;
   optOut.verbose = 0;
   WriteBrik(data_rescaled,data_to_rescale_info,optOut);
   voxgamfitwithheight({strcat(data_to_rescale_info.RootName,'.BRIK')}, dsmask);
   for briknum=0:3
        movefile(sprintf('gamfitb%d+orig.BRIK',briknum),sprintf('gamfitb%d_%s+tlrc.BRIK',briknum,data_to_rescale_info.RootName(1:end-18)));
        movefile(sprintf('gamfitb%d+orig.HEAD',briknum),sprintf('gamfitb%d_%s+tlrc.HEAD',briknum,data_to_rescale_info.RootName(1:end-18)));
   end  
end