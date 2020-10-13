function rl_nf_write_ICC_maps_computed_with_matlab()
% rl_nf_write_ICC_maps_computed_with_matlab.m will write maps computed
% by rl_nf_compute_ICCs under a format of .BRIK & .HEAD used by afni
% This script assumes we have the AFNI matlab library and others function of the reliability toolbox in our path
% This code has been described in Compere et al. (2020)


maps = struct('name',{'amplitude_active_data_BV_style',...
                      'amplitude_active_data_standard',...
                      'amplitude_control_data_BV_style',...
                      'amplitude_control_data_standard',...
                      'area_under_the_curve_active_data_BV_style',...
                      'area_under_the_curve_active_data_standard',...
                      'area_under_the_curve_control_data_BV_style',...
                      'area_under_the_curve_control_data_standard',...
                      'canonical_amplitude_active_data_BV_style',...
                      'canonical_amplitude_active_data_standard',...
                      'canonical_amplitude_control_data_BV_style',...
                      'canonical_amplitude_control_data_standard',...
                      'height_active_data_BV_style',...
                      'height_active_data_standard',...
                      'height_control_data_BV_style',...
                      'height_control_data_standard',...
                      'onset_delay_active_data_BV_style',...
                      'onset_delay_active_data_standard',...
                      'onset_delay_control_data_BV_style',...
                      'onset_delay_control_data_standard',...
                      'rise_decay_rate_active_data_BV_style',...
                      'rise_decay_rate_active_data_standard',...
                      'rise_decay_rate_control_data_BV_style',...
                      'rise_decay_rate_control_data_standard'});     

                  dsmask='Lamygdala_resampled+tlrc';
                      [tmp,V,Info]=BrikLoad(dsmask);
                  
                  for fl=1:numel(maps)
                      load(sprintf('corr_%s.mat',maps(fl).name));
                      % Write ICC mapts
                      Opt.Scale = 1;
                      Opt.verbose = 0;
                      Info.RootName=sprintf('icc_%s_computed_with_matlab+tlrc',maps(fl).name);
                      Opt.Prefix = sprintf('icc_%s_computed_with_matlab+tlrc',maps(fl).name);
                      Info.BRICK_TYPES=3;
                      
                      [err, ErrMessage, Info] = WriteBrik (corr.icc, Info, Opt);
                      
                      % Write semi partial correlations with no covariates
                      % maps
                      Info.RootName=sprintf('s_p_corr_without_covariates_%s_computed_with_matlab+tlrc',maps(fl).name);
                      Opt.Prefix = sprintf('s_p_corr_without_covariates_%s_computed_with_matlab+tlrc',maps(fl).name);
                      
                      [err, ErrMessage, Info] = WriteBrik (corr.s_p_r_without_covariates, Info, Opt);
                      
                      % Write semi partial correlations with covariates
                      % maps
                      Info.RootName=sprintf('s_p_corr_with_covariates_%s_computed_with_matlab+tlrc',maps(fl).name);
                      Opt.Prefix = sprintf('s_p_corr_with_covariates_%s_computed_with_matlab+tlrc',maps(fl).name);
                      
                      [err, ErrMessage, Info] = WriteBrik (corr.s_p_r_with_covariates, Info, Opt);
                  end