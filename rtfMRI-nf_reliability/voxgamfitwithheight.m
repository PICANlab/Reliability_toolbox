function res=voxgamfitwithheight(dslist, dsmask)
% does a gamma variate fit with height at  every voxel
% returns a 4-subbrik 
% e.g., res=voxgamfitwithheight({'NegTransient_CensorNeutandFixCoeffs_1022_02+orig.BRIK'},'/data/Siegle/ABM/struct/refbrain/minicolin+orig.BRIK')
if nargin<1, dslist=[]; end
if nargin<2, dsmask = '/data/Siegle/ATLAS/minicolin+orig'; end


addpath ('/data/Siegle/matlab/afni_matlab');
addpath ('/data/Siegle/matlab/physiotools');


curpath=pwd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the mask
[mask,dmaskinfo]=BrikLoad(dsmask);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do gamma fits

fprintf('======================================\nFitting data\n');

for ct=1:length(dslist)
  % read the file
  fname=char(dslist(ct));
  fprintf('loading %s\n',fname);
  [spath,ofilename]=strippath(char(dslist(ct)));
  if ~isempty(spath)
    cd(spath);
  end
  [ds1,dinfo]=BrikLoad(ofilename);
  
  allgam=zeros(size(mask,1),size(mask,2),size(mask,3),4);
  
  for x=1:size(mask,1)
    fprintf('\n');
    for y=1:size(mask,2) 
      fprintf('.'); 
      for z=1:size(mask,3)
	if mask(x,y,z)
	  %fprintf('%d,%d,%d...',x,y,z);
	  gam=gamfitwithheight(squeeze(ds1(x,y,z,:)),0);
	  % parms are [onset_delay rise_decay_rate height] - have
          % these and the correlation of the gamma fit with the raw data
	  allgam(x,y,z,:)=[gam.parms gam.r];
	end
      end
    end
  end

  % write the data out to its path
  if isempty(spath)
    fname=fprintf('gamma_%s',char(dslist(ct)));
  else
    fname=fprintf('%sgamma_%s',spath,ofilename);
  end
  
  fprintf('Writing %s\n',char(fname));  
  if ~isempty(spath)
    cd(spath);
  end
  
  for briknum=1:4
    outinfo=dmaskinfo;
    outinfo.RootName=sprintf('gamfitb%d+orig',briknum-1);
    %  outinfo.DATASET_RANK(2)=4;
    %  outinfo.BRICK_TYPES=[3 3 3 3];
    %  outinfo.BRICK_STATS=[0 1 1 1 1];
    %  outinfo.BRICK_FLOAT_FACS=[0 0 0 0];
   optOut.Scale = 1;
   optOut.Prefix = outinfo.RootName(1:end-5);
   optOut.verbose = 0;
  
   allgam(find(isnan(allgam)))=0;
  
   WriteBrik(squeeze(allgam(:,:,:,briknum)),outinfo,optOut);
  end
  
  
  cd(curpath);
end

