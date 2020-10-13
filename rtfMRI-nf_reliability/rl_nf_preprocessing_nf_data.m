function rl_nf_preprocessing_nf_data()
% On data preprocessed with the BV style pipeline, emulate the analysis
% done by Turbo BrainVoyager: concatenate all the happy blocks minus mean
% value in previous rest block and create stim_file for first level
% analysis to come
% This script assumes we have the AFNI matlab library and others function of the reliability toolbox in our path
% This code has been described in Compere et al. (2020)

% In active participants, on all functionals for the transfer runs
participants=dir('data_for_voxel_wise_reliability_right_data/Active');
cd(sprintf('%s',participants(1).folder))
for subj=3:size(participants,1)
    cd(fullfile(sprintf('%s',participants(subj).name),'Preprocessing_BV_style'))
    functionals=dir('*/*/all_runs.*.BRIK');
    % On all functionals in the tlrc space
    for func=1:size(functionals,1)
        [~, signal_to_preproc, Info, ~]=BrikLoad(fullfile(functionals(func).folder,functionals(func).name));
        [x, y, z, n]=size(signal_to_preproc);
        for dimx = 1:x
            for dimy = 1:y
                for dimz = 1:z
                    % take out dummy scans and concatenate the signal of each happy
                    % blocks minus previous rest block mean
                    preproc=squeeze((signal_to_preproc(dimx,dimy,dimz,:)));
                    postproc=[preproc(21:40,1)-mean(preproc(1:20,1));...
                        preproc(101:120,1)-mean(preproc(61:80,1));...
                        preproc(141:160,1)-mean(preproc(121:140,1));...
                        preproc(201:220,1)-mean(preproc(181:200,1))];
                    signal_preproc(dimx,dimy,dimz,:)=postproc;
                end
            end
        end
        % Write new functional preprocessed
        Opt.Scale = 1;
        Opt.verbose = 0;
        Info.BRICK_LABS='#0~#1~#2~#3~#4~#5~#6~#7~#8~#9~#10~#11~#12~#13~#14~#15~#16~#17~#18~#19~#20~#21~#22~#23~#24~#25~#26~#27~#28~#29~#30~#31~#32~#33~#34~#35~#36~#37~#38~#39~#40~#41~#42~#43~#44~#45~#46~#47~#48~#49~#50~#51~#52~#53~#54~#55~#56~#57~#58~#59~#60~#61~#62~#63~#64~#65~#66~#67~#68~#69~#70~#71~#72~#73~#74~#75~#76~#77~#78~#79';
        Info.BRICK_TYPES=ones(size(signal_preproc,4),1);
        Info.RootName=sprintf('%s/Preprocessed_%s',functionals(func).folder,functionals(func).name(1:end-5));
        Opt.Prefix = sprintf('%s/Preprocessed_%s',functionals(func).folder,functionals(func).name(1:end-5));
        Info.DATASET_RANK=[3 80 0 0 0 0 0 0];
        Info.BRICK_STATS=Info.BRICK_STATS(1:160);
        Info.BRICK_FLOAT_FACS=Info.BRICK_FLOAT_FACS(1:80);
        Info.TAXIS_NUMS(1)=80;
        [err, ErrMessage, Info] = WriteBrik (signal_preproc, Info, Opt);
    end
    cd ../..
end
cd ../..
% In control participants, on all functionals for the baseline runs
participants=dir('data_for_voxel_wise_reliability_right_data/Control');
cd(sprintf('%s',participants(1).folder))
for subj=3:size(participants,1)
    cd(fullfile(sprintf('%s',participants(subj).name),'Preprocessing_BV_style'))
    functionals=dir('*/*/all_runs.*.BRIK');
    % On all functionals in the tlrc space
    for func=1:size(functionals,1)
        [~, signal_to_preproc, Info, ~]=BrikLoad(fullfile(functionals(func).folder,functionals(func).name));
        [x, y, z, n]=size(signal_to_preproc);
        for dimx = 1:x
            for dimy = 1:y
                for dimz = 1:z
                    % take out dummy scans and concatenate the signal of each happy
                    % blocks minus previous rest block mean
                    preproc=squeeze((signal_to_preproc(dimx,dimy,dimz,:)));
                    postproc=[preproc(21:40,1)-mean(preproc(1:20,1));...
                        preproc(101:120,1)-mean(preproc(61:80,1));...
                        preproc(141:160,1)-mean(preproc(121:140,1));...
                        preproc(201:220,1)-mean(preproc(181:200,1))];
                    signal_preproc(dimx,dimy,dimz,:)=postproc;
                end
            end
        end
        % Write new functional preprocessed
        Opt.Scale = 1;
        Opt.verbose = 0;
        Info.BRICK_LABS='#0~#1~#2~#3~#4~#5~#6~#7~#8~#9~#10~#11~#12~#13~#14~#15~#16~#17~#18~#19~#20~#21~#22~#23~#24~#25~#26~#27~#28~#29~#30~#31~#32~#33~#34~#35~#36~#37~#38~#39~#40~#41~#42~#43~#44~#45~#46~#47~#48~#49~#50~#51~#52~#53~#54~#55~#56~#57~#58~#59~#60~#61~#62~#63~#64~#65~#66~#67~#68~#69~#70~#71~#72~#73~#74~#75~#76~#77~#78~#79';
        Info.BRICK_TYPES=ones(size(signal_preproc,4),1);
        Info.RootName=sprintf('%s/Preprocessed_%s',functionals(func).folder,functionals(func).name(1:end-5));
        Opt.Prefix = sprintf('%s/Preprocessed_%s',functionals(func).folder,functionals(func).name(1:end-5));
        Info.DATASET_RANK=[3 80 0 0 0 0 0 0];
        Info.BRICK_STATS=Info.BRICK_STATS(1:160);
        Info.BRICK_FLOAT_FACS=Info.BRICK_FLOAT_FACS(1:80);
        Info.TAXIS_NUMS(1)=80;
        [err, ErrMessage, Info] = WriteBrik (signal_preproc, Info, Opt);
    end
    cd ../..
end
cd ../..
% Create stim_file for preprocessed data
stim_file=zeros(1,80);
stim_file(1:20:80)=1;
dlmwrite('stim_file_for_preprocessed_data.1D',stim_file)
stim_file_rest=zeros(1,260);
stim_file_rest(1:60:260)=1;
dlmwrite('stim_file_for_original_data_rest.1D',stim_file_rest)
stim_file_happy=zeros(1,260);
stim_file_happy(21:60:260)=1;
dlmwrite('stim_file_for_original_data_happy.1D',stim_file_happy)
stim_file_count=zeros(1,260);
stim_file_count(41:60:260)=1;
dlmwrite('stim_file_for_original_data_count.1D',stim_file_count)
