function rl_compute_ICCs()
% main computations for the reliability paper (ICCs and semi-partial
% correlations with and without covariates)
% This script assumes we have the AFNI matlab library and others function of the reliability toolbox in our path
% This code has been described in Compere et al. (2020)


global data_master dsmask mask dat dat_pat;

data_master = struct('first_level_modelization',{'canonical_amplitude',...
                                                 'amplitude',...
                                                 'area_under_the_curve',...
                                                 'onset_delay',...
                                                 'rise_decay_rate',...
                                                 'height'},...
                     'look_for_controls_data',{'ssanal/sret/ssregs_canonical_amplitude/sret2*-negtrials.ssreg+orig.BRIK',...
                                               'ssanal/sret/ssregs_auc/sret2*AMP+orig.BRIK',...
                                               'ssanal/sret/ssregs_auc/sret2*AUC+orig.BRIK',...
                                               'ssanal/sret/ssregs_gammavar_computed_with_matlab/gamfitb0_sret2*-negtrial+orig.BRIK',...
                                               'ssanal/sret/ssregs_gammavar_computed_with_matlab/gamfitb1_sret2*-negtrial+orig.BRIK',...
                                               'ssanal/sret/ssregs_gammavar_computed_with_matlab/gamfitb2_sret2*-negtrial+orig.BRIK'},...
                     'look_for_patients_cohort1_data',{'ssanal/sret/ssregs_canonical_amplitude/sret5*-negtrials.ssreg+orig.BRIK',...
                                                       'ssanal/sret/ssregs_auc/sret5*AMP+orig.BRIK',...
                                                       'ssanal/sret/ssregs_auc/sret5*AUC+orig.BRIK',...
                                                       'ssanal/sret/ssregs_gammavar_computed_with_matlab/gamfitb0_sret5*-negtrial+orig.BRIK',...
                                                       'ssanal/sret/ssregs_gammavar_computed_with_matlab/gamfitb1_sret5*-negtrial+orig.BRIK',...
                                                       'ssanal/sret/ssregs_gammavar_computed_with_matlab/gamfitb2_sret5*-negtrial+orig.BRIK'},...
                     'look_for_patients_cohort2_data',{'ssanal/sret/ssregs_canonical_amplitude/sret6*-negtrials.ssreg+orig.BRIK',...
                                                       'ssanal/sret/ssregs_auc/sret6*AMP+orig.BRIK',...
                                                       'ssanal/sret/ssregs_auc/sret6*AUC+orig.BRIK',...
                                                       'ssanal/sret/ssregs_gammavar_computed_with_matlab/gamfitb0_sret6*-negtrial+orig.BRIK',...
                                                       'ssanal/sret/ssregs_gammavar_computed_with_matlab/gamfitb1_sret6*-negtrial+orig.BRIK',...
                                                       'ssanal/sret/ssregs_gammavar_computed_with_matlab/gamfitb2_sret6*-negtrial+orig.BRIK'},...
                     'bucket',{3,...
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
                                  3700,...
                                  37000,...
                                  -1000,...
                                  -100000000000000,...
                                  -1000});
                              


dsmask = 'binminicolin+orig.BRIK';
% load the mask
[mask]=BrikLoad(dsmask);

%% load prerequisite files
dat=load('dataTable_all_sample_numeric_imputed_noheader.txt');
dat_pat=load('dataTable_patients_numeric_imputed_noheader.txt');

for fl=1:numel(data_master)
  rl_correlation_helper(fl);
end

for fl=1:numel(data_master)
  rl_correlation_helper2(fl);
end

