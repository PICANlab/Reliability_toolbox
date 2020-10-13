function rl_which_model_best()
% Compare semi partial correlations with and without covariates across
% different first model for each ROI independently
% This script assumes we have the AFNI matlab library and others function of the reliability toolbox in our path
% This code has been described in Compere et al. (2020)

ROIs={'binminicolin_Amygdala_MNIspace_epan+orig.BRIK',...
       'binminicolin_DLPFC_epan+orig.BRIK',...
       'binminicolin_rACC_epan+orig.BRIK',...
       'sgACC_mask+orig.BRIK',...
       'sg_ACC_liberally_thresholded+orig.BRIK'};
ROI_name={'Amygdala',...
       'DLPFC',...
       'rACC',...
       'sgACC conservatively thresholded',...
       'sgACC liberally thresholded'};
   
sr_maps_controls_and_patients={'s_p_corr_without_covariates_canonical_amplitude_all_sample_computed_with_matlab+orig.BRIK',...
                               's_p_corr_with_covariates_canonical_amplitude_all_sample_computed_with_matlab+orig.BRIK',...
                               's_p_corr_without_covariates_amplitude_all_sample_computed_with_matlab+orig.BRIK',...
                               's_p_corr_with_covariates_amplitude_all_sample_computed_with_matlab+orig.BRIK',...
                               's_p_corr_without_covariates_area_under_the_curve_all_sample_computed_with_matlab+orig.BRIK',...
                               's_p_corr_with_covariates_area_under_the_curve_all_sample_computed_with_matlab+orig.BRIK',...
                               's_p_corr_without_covariates_onset_delay_all_sample_computed_with_matlab+orig.BRIK',...
                               's_p_corr_with_covariates_onset_delay_all_sample_computed_with_matlab+orig.BRIK',...
                               's_p_corr_without_covariates_rise_decay_rate_all_sample_computed_with_matlab+orig.BRIK',...
                               's_p_corr_with_covariates_rise_decay_rate_all_sample_computed_with_matlab+orig.BRIK',...
                               's_p_corr_without_covariates_height_all_sample_computed_with_matlab+orig.BRIK',...
                               's_p_corr_with_covariates_height_all_sample_computed_with_matlab+orig.BRIK'};
sr_maps_patients={'s_p_corr_without_covariates_canonical_amplitude_patients_computed_with_matlab+orig.BRIK',...
                  's_p_corr_with_covariates_canonical_amplitude_patients_computed_with_matlab+orig.BRIK',...
                  's_p_corr_without_covariates_amplitude_patients_computed_with_matlab+orig.BRIK',...
                  's_p_corr_with_covariates_amplitude_patients_computed_with_matlab+orig.BRIK',...
                  's_p_corr_without_covariates_area_under_the_curve_patients_computed_with_matlab+orig.BRIK',...
                  's_p_corr_with_covariates_area_under_the_curve_patients_computed_with_matlab+orig.BRIK',...
                  's_p_corr_without_covariates_onset_delay_patients_computed_with_matlab+orig.BRIK',...
                  's_p_corr_with_covariates_onset_delay_patients_computed_with_matlab+orig.BRIK',...
                  's_p_corr_without_covariates_rise_decay_rate_patients_computed_with_matlab+orig.BRIK',...
                  's_p_corr_with_covariates_rise_decay_rate_patients_computed_with_matlab+orig.BRIK',...
                  's_p_corr_without_covariates_height_patients_computed_with_matlab+orig.BRIK',...
                  's_p_corr_with_covariates_height_patients_computed_with_matlab+orig.BRIK'};
sr_name={'Canonical amplitude without covariates',...
         'Canonical amplitude with covariates',...
         'Amplitude without covariates',...
         'Amplitude with covariates',...
         'Area under the curve without covariates',...
         'Area under the curve with covariates',...
         'Onset delay without covariates',...
         'Onset delay with covariates',...
         'Rise decay without covariates',...
         'Rise decay with covariates',...
         'Height without covariates',...
         'Height with covariates'};

     % For the entire sample
for ROI=1:numel(ROIs)
    [mask]=BrikLoad(ROIs{ROI});
    mask(find(mask))=1;
    values_sr_in_ROI_all_models_controls_and_patients=zeros(numel(find(mask)),numel(sr_maps_controls_and_patients));
    for sr_map=1:numel(sr_maps_controls_and_patients)
        sr=BrikLoad(sr_maps_controls_and_patients{sr_map});
        values_sr_in_ROI=sr.*mask;
        values_sr_in_ROI_all_models_controls_and_patients(:,sr_map)=values_sr_in_ROI(find(mask));
        % Check if distribution of semi partial correlations is normal
        if any(~isnan(values_sr_in_ROI(find(mask))))
            non_normal_distribution_in_ROI_all_models_controls_and_patients(sr_map)=kstest(values_sr_in_ROI(find(mask))); % result h is 1 if the test rejects the null hypothesis that the data come from a standard normal distribution, 0 otherwise
        else
            non_normal_distribution_in_ROI_all_models_controls_and_patients(sr_map)=1;
        end
    end
    % Plot the distribution of the reliability estimates in different
    % models for comparaison and the p value of the Kruskal-Wallis test
    figure(ROI); 
    multihistimg(values_sr_in_ROI_all_models_controls_and_patients,100,sr_name,1);
    title(sprintf('Patients and controls: %s',ROI_name{ROI}));
    % Test
   if any(find(non_normal_distribution_in_ROI_all_models_controls_and_patients==1)) && ~any(find(non_normal_distribution_in_ROI_all_models_controls_and_patients==0)) 
       [p,tbl,stat]=kruskalwallis(values_sr_in_ROI_all_models_controls_and_patients,sr_name,'off');
       hold on;
       xlabel(sprintf('Kruskal-Wallis: p=%1.3f',p));
       hold off
       % If Kruskal-Wallis test is significant, do post hoc comparison
       if p<0.05
           figure(ROI+1);
           [c,m]=multcompare(stat);  
       end
   elseif any(find(non_normal_distribution_in_ROI_all_models_controls_and_patients==0)) && ~any(find(non_normal_distribution_in_ROI_all_models_controls_and_patients==1))
       fprintf('One-wayANOVA')
   else
       fprintf('mixed distribution\n')
   end
end


% Does the same for the sample of patients only
for ROI=1:numel(ROIs)
    [mask]=BrikLoad(ROIs{ROI});
    mask(find(mask))=1;
    values_sr_in_ROI_all_models_patients=zeros(numel(find(mask)),numel(sr_maps_patients));
    for sr_map=1:numel(sr_maps_patients)
        sr=BrikLoad(sr_maps_patients{sr_map});
        values_sr_in_ROI=sr.*mask;
        values_sr_in_ROI_all_models_patients(:,sr_map)=values_sr_in_ROI(find(mask));
        if any(~isnan(values_sr_in_ROI(find(mask))))
            non_normal_distribution_in_ROI_all_models_patients(sr_map)=kstest(values_sr_in_ROI(find(mask))); % result h is 1 if the test rejects the null hypothesis that the data come from a standard normal distribution, 0 otherwise
        else
            non_normal_distribution_in_ROI_all_models_patients(sr_map)=1;
        end
    end
    figure(ROI+5);
    multihistimg(values_sr_in_ROI_all_models_patients,100,sr_name,1);   
    title(sprintf('Patients: %s',ROI_name{ROI}));
   % Test
   if any(find(non_normal_distribution_in_ROI_all_models_patients==1)) && ~any(find(non_normal_distribution_in_ROI_all_models_patients==0)) 
       [p,tbl,stat]=kruskalwallis(values_sr_in_ROI_all_models_patients,sr_name,'off');
       hold on;
       xlabel(sprintf('Kruskal-Wallis: p=%1.3f',p));
       hold off
       if p<0.05
           figure(ROI+6)
           [c,m]=multcompare(stat);
       end    
   elseif any(find(non_normal_distribution_in_ROI_all_models_patients==0)) && ~any(find(non_normal_distribution_in_ROI_all_models_patients==1))
       fprintf('One-wayANOVA')
   else
       fprintf('mixed distribution\n')
   end
end