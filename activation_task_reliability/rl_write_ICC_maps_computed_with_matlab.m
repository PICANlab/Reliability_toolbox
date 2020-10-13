function rl_write_ICC_maps_computed_with_matlab()
% Write maps of reliability estimates computed with rl_compute_ICCs.m in
% .BRIK and .HEAD format
% This script assumes we have the AFNI matlab library and others function of the reliability toolbox in our path
% This code has been described in Compere et al. (2020)

maps = struct('name',{'canonical_amplitude_all_sample',...
                      'canonical_amplitude_patients',...
                      'amplitude_all_sample',...
                      'amplitude_patients',...
                      'area_under_the_curve_all_sample',...
                      'area_under_the_curve_patients',...
                      'onset_delay_all_sample',...
                      'onset_delay_patients',...
                      'rise_decay_rate_all_sample',...
                      'rise_decay_rate_patients',...
                      'height_all_sample',...
                      'height_patients'});     

                  dsmask='binminicolin+orig';
                      [tmp,V,Info]=BrikLoad(dsmask);
                 
                  for fl=1:numel(maps)
                      load(sprintf('corr_%s.mat',maps(fl).name));
                      % Write ICCs maps 
                      Opt.Scale = 1;
                      Opt.verbose = 0;
                      Info.RootName=sprintf('icc_%s_computed_with_matlab',maps(fl).name);
                      Opt.Prefix = sprintf('icc_%s_computed_with_matlab',maps(fl).name);
                      
                      [err, ErrMessage, Info] = WriteBrik (corr.icc, Info, Opt)
                    
                    % Write maps of semi partial correlations without
                    % covariates
                    Info.RootName=sprintf('s_p_corr_without_covariates_%s_computed_with_matlab',maps(fl).name);
                    Opt.Prefix = sprintf('s_p_corr_without_covariates_%s_computed_with_matlab',maps(fl).name);
                      
                    [err, ErrMessage, Info] = WriteBrik (corr.s_p_r_without_covariates, Info, Opt);
                      
                    % Write maps of semi partial correlations with
                    % covariates
                    Info.RootName=sprintf('s_p_corr_with_covariates_%s_computed_with_matlab',maps(fl).name);
                    Opt.Prefix = sprintf('s_p_corr_with_covariates_%s_computed_with_matlab',maps(fl).name);
                      
                    [err, ErrMessage, Info] = WriteBrik (corr.s_p_r_with_covariates, Info, Opt);
                  end

