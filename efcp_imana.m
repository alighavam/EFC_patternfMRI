function varargout = efcp_imana(what, varargin)
    % Use a different baseDir when using your local machine or the cbs
    % server. Add more directory if needed. Use single quotes ' and not
    % double quotes " because some spm function raise error with double
    % quotes
    if ismac
        baseDir = '/Volumes/Diedrichsen_data$/data/Chord_exp/EFC_patternfMRI';
    elseif isunix
        baseDir = '/cifs/diedrichsen/data/Chord_exp/EFC_patternfMRI';
    else
        disp('Running on Windows or another OS');
    end
    
    bidsDir = 'BIDS'; % Raw data post AutoBids conversion
    anatomicalDir = 'anatomicals';
    imagingRawDir = 'imaging_data_raw'; % Temporary directory for raw functional data
    imagingDir = 'imaging_data'; % Preprocesses functional data
    fmapDir = 'fieldmaps'; % Fieldmap dir after moving from BIDS and SPM make fieldmap
    
    pinfo = dload(fullfile(baseDir,'participants.tsv'));

    % handling input args:
    sn = [];
    rtm = 0;
    prefix = 'u';
    vararginoptions(varargin,{'sn', 'rtm', 'prefix'})
    if isempty(sn)
        error('BIDS:move_unzip_raw_func -> ''sn'' must be passed to this function.')
    end

    % get participant row from participant.tsv
    participant_row = getrow(pinfo, pinfo.sn == sn);
    
    % get participant id:
    participant_id = participant_row.participant_id{1};
    
    % get BIDS session names:
    for i = 1:participant_row.num_ses
        % removes extra double quotations from the name:
        ses{i} = erase(participant_row.(sprintf('name_ses%d',i)){1}, '"');
    end
    
    switch(what)
        case 'BIDS:move_unzip_raw_func'
            % Moves, unzips and renames raw functional (BOLD) images from 
            % BIDS directory. After you run this function you will find
            % nRuns Nifti files named <subj_id>_run_XX.nii in the 
            % <project_id>/imaging_data_raw/<subj_id>/ directory.
            
            for i = 1:length(ses)
                runs = spmj_dotstr2array(participant_row.(sprintf('run_ses%d',i)){1});
                func_tmp_name = sprintf('%s.nii.gz', participant_row.(sprintf('func_name_ses%d',i)){1});
                
                % loop on runs of sess:
                for run = runs
                    % add run number to the name of the functional file:
                    func_name = replace(func_tmp_name,'XX',sprintf('%.02d',run));
                    
                    % path to the subj func file:
                    func_file = fullfile(baseDir,bidsDir,sprintf('sub-%s',participant_id),sprintf('ses-%s',ses{i}),'func',func_name);
                    
                    % destination path:
                    output_folder = fullfile(baseDir,imagingRawDir,participant_id,sprintf('ses-%s',ses{i}));
                    output_file = fullfile(output_folder,[participant_id sprintf('_run_%.02d.nii.gz',run)]);
                    
                    if ~exist(output_folder,"dir")
                        mkdir(output_folder);
                    end
                    
                    % copy file to destination:
                    [status,msg] = copyfile(func_file,output_file);
                    if ~status  
                        error('FUNC:move_unzip_raw_func -> subj %d raw functional (BOLD) was not moved from BIDS to the destenation:\n%s',sn,msg)
                    end
                    
                    % unzip the .gz files to make usable for SPM:
                    gunzip(output_file);

                    % delete the compressed file:
                    delete(output_file);
                end
            end
    
        case 'BIDS:move_unzip_raw_fmap'
            % Moves, unzips and renames raw fmap images from BIDS
            % directory. After you run this function you will find two
            % files named <subj_id>_phase.nii and <subj_id>_magnitude.nii
            % in the <project_id>/fieldmaps/<subj_id>/ directory. The
            % <subj_id>_phase.nii contains phase information derived from
            % the MRI signal, reflecting local magnetic field
            % inhomogeneities caused by factors such as tissue
            % susceptibility differences (e.g., at air-tissue interfaces in
            % the nasal cavities or sinuses). This phase information can be
            % used to compute a fieldmap, which is essential for correcting
            % geometric distortions (unwarping) in other MRI sequences.

            for i = 1:length(ses)
                fmap_magnitude_name = participant_row.fmap_magnitude{1};
                fmap_phase_name = participant_row.fmap_phase{1};
                
                magnitude = dir(fullfile(baseDir, bidsDir, sprintf('sub-%s',participant_id), sprintf('ses-%s',ses{i}), 'fmap', ['*_' fmap_magnitude_name '.nii.gz']));
                magnitude_path = fullfile(magnitude.folder, magnitude.name);
                phase = dir(fullfile(baseDir, bidsDir, sprintf('sub-%s',participant_id), sprintf('ses-%s',ses{i}), 'fmap', ['*_' fmap_phase_name '.nii.gz']));
                phase_path = fullfile(phase.folder, phase.name);

                output_dir = fullfile(baseDir, fmapDir, participant_id, sprintf('ses-%s',ses{i}));
                if ~exist(output_dir,'dir')
                    mkdir(output_dir)
                end
                output_magnitude_path = fullfile(output_dir, sprintf('%s_magnitude.nii.gz',participant_id));
                output_phase_path = fullfile(output_dir, sprintf('%s_phase.nii.gz',participant_id));

                [status,msg] = copyfile(magnitude_path,output_magnitude_path);
                [status,msg] = copyfile(phase_path,output_phase_path);

                % unzip the .gz files to make usable for SPM:
                gunzip(output_magnitude_path);
                delete(output_magnitude_path);

                gunzip(output_phase_path);
                delete(output_phase_path);
            end
        
        case 'FUNC:fix_dtype'
            % After Khanlab gradcorrect, data is transformed from uint16 to
            % int16. This makes some positive values in some images to be
            % negative. This function fixes this in imaging_data_raw.
            
            % Get all .nii files in all subdirectories
            file_list = dir(fullfile(baseDir, imagingRawDir, participant_id, '**', '*.nii'));
            
            % Loop through each found file
            for i = 1:length(file_list)
                % Get the full path of the file
                file_path = fullfile(file_list(i).folder, file_list(i).name);
                
                % Perform the desired operation on the file
                fprintf('Processing file: %s\n', file_path);
                
                epi = niftiread(file_path);
                if ~isa(epi, 'uint16')
                    epi_info = niftiinfo(file_path);
                    epi_uint16 = typecast(epi(:), 'uint16');
                    epi_uint16 = reshape(epi_uint16, size(epi));
                    epi_info.Datatype = 'uint16';
                    niftiwrite(epi_uint16, file_path, epi_info);
                end
            end
            
        case 'FUNC:make_fmap'                
            % Differences in magnetic susceptibility between tissues (e.g.,
            % air-tissue or bone-tissue interfaces) can cause
            % inhomogeneities in the magnetic field. These inhomogeneities
            % result in spatial distortions along the phase-encoding
            % direction, which is the direction in which spatial location
            % is encoded using a phase gradient. To account for these
            % distortions, this step generates a Voxel Displacement Map
            % (VDM) for each run, saved as files named
            % vdm5_sc<subj_id>_phase_run_XX.nii in the fieldmap directory.
            % 
            % The VDM assigns a value in millimeters to each voxel,
            % indicating how far it should be shifted along the
            % phase-encoding direction to correct for the distortion. If
            % you open the VDM in FSLeyes, you will notice that the
            % distortion is particularly strong in the temporal lobe due to
            % proximity to the nasal cavities, where significant
            % differences in magnetic susceptibility occur.
            % 
            % In the fieldmap directory, you will also find the intermediate
            % files bmask<subj_id>_magnitude.nii and
            % fpm_sc<subj_id>_phase.nii that are used for VDM calculation
            % 
            % In the imaging_data_raw directory, you will find unwarped
            % functional volumes named u<subj_id>_run_XX.nii. These
            % correspond to the corrected first volume of each functional
            % run. Open them in FSL to inspect how the distortion was
            % corrected using the VDM (this step is for quality checking;
            % the actual unwarping is performed in a later step).
            % 
            % In addition, the imaging_raw_data directory contains the
            % intermediate file wfmag_<subj_id>_run_XX.nii that is
            % necessary to perform unwarping in eah run.
            
            for i = 1:length(ses)
                epi_path = fullfile(baseDir, imagingRawDir, participant_id, sprintf('ses-%s',ses{i}));
                
                epi_files = dir(fullfile(epi_path, [participant_id '_run_*.nii']));
                
                % params:
                et1 = 0.00408*1000;
                et2 = 0.0051*1000;
                tert = 0.000334996*90*1000;
                
                spmj_makefieldmap(fullfile(baseDir, fmapDir, participant_id, sprintf('ses-%s',ses{i})), ...
                                  sprintf('%s_magnitude.nii', participant_id),...
                                  sprintf('%s_phase.nii', participant_id),...
                                  'phase_encode', -1, ... % It's -1 (A>>P) or 1 (P>>A) and can be found in imaging sequence specifications
                                  'et1', et1, ...
                                  'et2', et2, ...
                                  'tert', tert, ...
                                  'func_dir',fullfile(baseDir, imagingRawDir, participant_id, sprintf('ses-%s',ses{i})),...
                                  'epi_files', {epi_files.name});
            end

        case 'FUNC:realign_unwarp'   
            for i = 1:length(ses)
                % Do spm_realign_unwarp
                epi_runs = dir(fullfile(baseDir, imagingRawDir, participant_id, sprintf('ses-%s',ses{i}), [participant_id '_run_*.nii']));
                fmap_runs = dir(fullfile(baseDir, fmapDir, participant_id, sprintf('ses-%s',ses{i}), ['vdm5_sc' participant_id '_phase_run_*.nii']));
                
                epi_list = {};
                fmap_list = {};
                for run = 1:length(epi_runs)
                    epi_list{end+1} = fullfile(epi_runs(run).folder, epi_runs(run).name);
                    fmap_list{end+1} = fullfile(fmap_runs(run).folder, fmap_runs(run).name);
                end
                
                spmj_realign_unwarp(epi_list, fmap_list, 'rtm', rtm);
            end
            
        case 'FUNC:inspect_realign'
            % looks for motion correction logs into imaging_data, needs to
            % be run after realigned images are moved there from
            % imaging_data_raw
            for i = 1:length(ses)
                rp_files = dir(fullfile(baseDir, imagingRawDir, participant_id, sprintf('ses-%s',ses{i}), ['rp_' participant_id '_run_*.txt']));
                rp_list = {}; % Initialize as an empty cell array
                for run = 1:length(rp_files)
                    rp_list{end+1} = fullfile(rp_files(run).folder, rp_files(run).name);
                end
                smpj_plot_mov_corr(rp_list)
            end

        case 'FUNC:move_realigned_images'          
            for i = 1:length(ses)
                % Move images created by realign(+unwarp) into imaging_data
                realigned_epi_files = dir(fullfile(baseDir, imagingRawDir, participant_id, sprintf('ses-%s', ses{i}), sprintf('%s%s_run_*.nii', prefix, participant_id)));
                rp_files = dir(fullfile(baseDir, imagingRawDir, participant_id, sprintf('ses-%s', ses{i}), sprintf('rp_%s_run_*.txt', participant_id)));
                mean_epi_files = dir(fullfile(baseDir, imagingRawDir, participant_id, sprintf('ses-%s', ses{i}), '*mean*.nii'));
                
                dest_dir = fullfile(baseDir, imagingDir, participant_id, sprintf('ses-%s', ses{i}));
                if ~exist(dest_dir,'dir')
                    mkdir(dest_dir)
                end
                
                % loop on runs of the session:
                for run = 1:length(realigned_epi_files)
                    % realigned (and unwarped) images names:
                    source = fullfile(realigned_epi_files(run).folder, realigned_epi_files(run).name);
                    
                    out_name = realigned_epi_files(run).name(length(prefix) + 1:end); % remove the prefix from realigned (and unwarped) file names
                    dest = fullfile(dest_dir, out_name);
                    % move to destination:
                    fprintf('moving and renaming imaging_data_raw/%s --> imaging_data/%s\n',realigned_epi_files(run).name,out_name)
                    [status,msg] = copyfile(source,dest);
                    if ~status  
                        error('FUNC:move_realigned_images -> %s',msg)
                    end
                    
                    % realign parameters names:
                    source = fullfile(rp_files(run).folder, rp_files(run).name);
                    dest = fullfile(dest_dir, rp_files(run).name);
                    % move to destination:
                    [status,msg] = copyfile(source,dest);
                    if ~status  
                        error('FUNC:move_realigned_images -> %s',msg)
                    end
                end
                
                % move the mean epi file:
                source = fullfile(mean_epi_files(1).folder, mean_epi_files(1).name);
                dest = fullfile(dest_dir, mean_epi_files(1).name);
                % move to destination:
                [status,msg] = copyfile(source, dest);
                if ~status  
                    error('BIDS:move_realigned_images -> %s',msg)
                end
            end
            
        case 'FUNC:meanimage_bias_correction'
            % EPI images often contain smooth artifacts caused by MRI
            % physics which make the intensity of signal from the same
            % tissue (e.g., grey matter, white matter) non-uniform. This
            % step perform bias correction and creates an image where the
            % signal from each tissue type is more uniform. This image is
            % then co-registered to the anatomical image. Bias correction
            % help make co-registration more accurate. If the realignment
            % was done with respect to the first volume of each run of each
            % session, the mean image will be calculated on the first run
            % of each session and will be called 'meanu*_run_01.nii'
            % ('mean' indicates the image is average of the volumes and 'u'
            % indicates it's unwarped). Therefore, we do the bias
            % correction on this file. But if you do the realignment to the
            % mean epi of every run, the generated mean file will be named
            % 'umeanepi_*' and we do the bias correction on this file. In
            % addition, this step generates five tissue probability maps
            % (c1-5) for grey matter, white matter, csf, bone and soft
            % tissue.
            for i = 1:length(ses)
                mean_epi = dir(fullfile(baseDir, imagingDir, participant_id, sprintf('ses-%s',ses{i}), 'mean*.nii'));
                P{1} = fullfile(mean_epi(1).folder, mean_epi(1).name);
                spmj_bias_correct(P);
            end
            
        case 'FUNC:coreg'                                                      
            % coregister rbumean image to anatomical image for each session

            % (1) Manually seed the functional/anatomical registration
            % - Open fsleyes
            % - Add anatomical image and b*mean*.nii (bias corrected mean) image to overlay
            % - click on the bias corrected mean image in the 'Overlay
            %   list' in the bottom left of the fsleyes window.
            %   list to highlight it.
            % - Open tools -> Nudge
            % - Manually adjust b*mean*.nii image to the anatomical by 
            %   changing the 6 paramters (tranlation xyz and rotation xyz) 
            %   and Do not change the scales! 
            % - When done, click apply and close the tool tab. Then to save
            %   the changes, click on the save icon next to the mean image 
            %   name in the 'Overlay list' and save the new image by adding
            %   'r' in the beginning of the name: rb*mean*.nii. If you don't
            %   set the format to be .nii, fsleyes automatically saves it as
            %   a .nii.gz so either set it or gunzip afterwards to make it
            %   compatible with SPM.
            
            % (2) Run automated co-registration to register bias-corrected meanimage to anatomical image
            for i = 1:length(ses)
                mean_epi = dir(fullfile(baseDir, imagingDir, participant_id, sprintf('ses-%s',ses{i}), ['bmean' prefix, participant_id '_run_*.nii']));
                
                J.source = {fullfile(mean_epi(1).folder, mean_epi(1).name)}; 
                J.ref = {fullfile(baseDir, anatomicalDir, participant_id, [participant_id, '_T1w','.nii'])};
                J.other = {''};
                J.eoptions.cost_fun = 'nmi';
                J.eoptions.sep = [4 2];
                J.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                J.eoptions.fwhm = [7 7];
                matlabbatch{1}.spm.spatial.coreg.estimate=J;
                spm_jobman('run',matlabbatch);
            end
                
        case 'FUNC:make_samealign'
            % align to registered bias corrected mean image of each session
            % (rb*mean*.nii). Alignment happens only by changing the
            % transform matrix in the header files of the functional 4D
            % .nii files to the transform matrix that aligns them to
            % anatomical. The reason that it works is: 1) in the
            % realignment (+unwarping) process, we have registered every
            % single volume of every single run to the first volume of the
            % first run of the session (rtm=0). 2) In the same step, for each
            % session, a mean functional image (meanepi*.nii or meanu*.nii
            % based on the rtm option) was generated. This mean image is
            % alread in the space of all the functional volumes. Later we
            % coregister this image to the anatomical space. Therefore, if
            % we change the transformation matrices of all the functional
            % volumes to the transform matrix of the coregistered image,
            % they will all tranform into the anatomical coordinates space.
            
            for i = 1:length(ses)
                epi_files = dir(fullfile(baseDir, imagingDir, participant_id, sprintf('ses-%s',ses{i}), [participant_id '_run_*.nii']));
                epi_list = {}; % Initialize as an empty cell array
                for run = 1:length(epi_files)
                    epi_list{end+1} = fullfile(epi_files(run).folder, epi_files(run).name);
                end
                
                % select the reference image:
                if rtm==0
                    mean_epi = dir(fullfile(baseDir, imagingDir, participant_id, sprintf('ses-%s',ses{i}), ['bmean' prefix participant_id '_run_*.nii']));
                    P{1} = fullfile(mean_epi(1).folder, mean_epi(1).name);
                else
                    P{1} = fullfile(baseDir,imagingDir,subj_id,['rb' prefix 'meanepi_' subj_id '.nii']);
                end
                
                % select images to be realigned:
                Q = {};
                for r = 1:length(epi_list)
                    for j = 1:pinfo.numTR
                        Q{end+1} = fullfile(sprintf('%s,%d', epi_list{r}, j));
                    end
                end
                
                spmj_makesamealign_nifti(char(P),char(Q));
            end
        
        case 'FUNC:make_maskImage'
            % Make mask images (noskull and gray_only) for 1st level glm
            for i = 1:length(ses)
                epi_files = dir(fullfile(baseDir, imagingDir, participant_id, sprintf('ses-%s',ses{i}), [participant_id '_run_*.nii']));
                epi_list = {}; % Initialize as an empty cell array
                for run = 1:length(epi_files)
                    epi_list{end+1} = fullfile(epi_files(run).folder, epi_files(run).name);
                end
                
                % bias corrected mean epi image:
                if rtm==0
                    mean_epi = dir(fullfile(baseDir, imagingDir, participant_id, sprintf('ses-%s',ses{i}), ['bmean' prefix participant_id '_run_*.nii']));
                    nam{1} = fullfile(mean_epi(1).folder, mean_epi(1).name);
                else
                    nam{1} = fullfile(baseDir,imagingDir,subj_id,['rb' prefix 'meanepi_' subj_id '.nii']);
                end
                nam{2} = fullfile(baseDir, anatomicalDir, participant_id, ['c1', participant_id, '_T1w','.nii']);
                nam{3} = fullfile(baseDir, anatomicalDir, participant_id, ['c2', participant_id, '_T1w','.nii']);
                nam{4} = fullfile(baseDir, anatomicalDir, participant_id, ['c3', participant_id, '_T1w','.nii']);
                spm_imcalc(nam, fullfile(baseDir, imagingDir, participant_id, sprintf('ses-%s',ses{i}), 'rmask_noskull.nii'), 'i1>1 & (i2+i3+i4)>0.2')
                
                source = fullfile(baseDir, imagingDir, participant_id, sprintf('ses-%s',ses{i}), 'rmask_noskull.nii');
                dest = fullfile(baseDir, anatomicalDir, participant_id, sprintf('rmask_noskull_ses-%s.nii',ses{i}));
                movefile(source, dest);
                
                % gray matter mask for covariance estimation
                % ------------------------------------------
                nam={};
                % nam{1}  = fullfile(imagingDir,subj_id{sn}, 'sess1', ['rb' prefix 'meanepi_' subj_id{sn} '.nii']);
    
                % IS THIS CHANGE CORRECT??
                % nam{1}  = fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)),sprintf('sess%d',sess), ['rb' prefix 'meanepi_' char(pinfo.subj_id(pinfo.sn==sn)) '.nii']);
                % bias corrected mean epi image:
                if rtm==0
                    mean_epi = dir(fullfile(baseDir, imagingDir, participant_id, sprintf('ses-%s',ses{i}), ['bmean' prefix participant_id '_run_*.nii']));
                    nam{1} = fullfile(mean_epi(1).folder, mean_epi(1).name);
                else
                    nam{1} = fullfile(baseDir,imagingDir,subj_id,['rb' prefix 'meanepi_' subj_id '.nii']);
                end
                nam{2} = fullfile(baseDir, anatomicalDir, participant_id, ['c1', participant_id, '_T1w','.nii']);
                spm_imcalc(nam, fullfile(baseDir, imagingDir, participant_id, sprintf('ses-%s',ses{i}), 'rmask_gray.nii'), 'i1>1 & i2>0.4')
                
                source = fullfile(baseDir, imagingDir, participant_id, sprintf('ses-%s',ses{i}), 'rmask_gray.nii');
                dest = fullfile(baseDir, anatomicalDir, participant_id, sprintf('rmask_gray_ses-%s.nii',ses{i}));
                movefile(source,dest);
            end
    
        case '==============CASES FOR SINGLE SESSION PREP=============='
        case 'FUNC:make_one_ses'
            cnt = 1;
            for i = 1:length(ses)
                epi_path = fullfile(baseDir, imagingRawDir, participant_id, sprintf('ses-%s',ses{i}));
                epi_files = dir(fullfile(epi_path, [participant_id '_run_*.nii']));

                fmap_path = fullfile(baseDir, fmapDir, participant_id, sprintf('ses-%s',ses{i}));
                fmap_files = dir(fullfile(fmap_path, ['vdm5_sc*_phase_run_*.nii']));

                % copy to main subject folder:
                for run = 1:length(epi_files)
                    output_epi_file = fullfile(baseDir, imagingRawDir, participant_id, sprintf('%s_run_%02d.nii',participant_id,cnt));
                    copyfile(fullfile(epi_files(run).folder,epi_files(run).name), output_epi_file);
                    fprintf('copied %s to %s\n',epi_files(run).name, sprintf('%s_run_%02d.nii',participant_id,cnt))

                    output_fmap_file = fullfile(baseDir, fmapDir, participant_id, sprintf('vdm5_sc%s_phase_run_%02d.nii',participant_id,cnt));
                    copyfile(fullfile(fmap_files(run).folder,fmap_files(run).name), output_fmap_file);
                    fprintf('copied %s to %s\n',fmap_files(run).name, sprintf('vdm5_sc%s_phase_run_%02d.nii',participant_id,cnt))

                    cnt = cnt+1;
                end
            end
        
        case 'FUNC:realign_unwarp_one_ses'
            % for i = 1:length(ses)
            %     % Do spm_realign_unwarp
            %     epi_runs = dir(fullfile(baseDir, imagingRawDir, participant_id, sprintf('ses-%s',ses{i}), [participant_id '_run_*.nii']));
            %     fmap_runs = dir(fullfile(baseDir, fmapDir, participant_id, sprintf('ses-%s',ses{i}), ['vdm5_sc' participant_id '_phase_run_*.nii']));
            % 
            %     epi_list = {};
            %     fmap_list = {};
            %     for run = 1:length(epi_runs)
            %         epi_list{end+1} = fullfile(epi_runs(run).folder, epi_runs(run).name);
            %         fmap_list{end+1} = fullfile(fmap_runs(run).folder, fmap_runs(run).name);
            %     end
            % 
            %     spmj_realign_unwarp(epi_list, fmap_list, 'rtm', rtm);
            % end
            % Do spm_realign_unwarp
            epi_runs = dir(fullfile(baseDir, imagingRawDir, participant_id, [participant_id '_run_*.nii']));
            fmap_runs = dir(fullfile(baseDir, fmapDir, participant_id, ['vdm5_sc' participant_id '_phase_run_*.nii']));
            
            epi_list = {};
            fmap_list = {};
            for run = 1:length(epi_runs)
                epi_list{end+1} = fullfile(epi_runs(run).folder, epi_runs(run).name);
                fmap_list{end+1} = fullfile(fmap_runs(run).folder, fmap_runs(run).name);
            end

            spmj_realign_unwarp(epi_list, fmap_list, 'rtm', rtm);
        
        case 'FUNC:move_realigned_images_ss'
            % Move images created by realign(+unwarp) into imaging_data
            realigned_epi_files = dir(fullfile(baseDir, imagingRawDir, participant_id, sprintf('%s%s_run_*.nii', prefix, participant_id)));
            rp_files = dir(fullfile(baseDir, imagingRawDir, participant_id, sprintf('rp_%s_run_*.txt', participant_id)));
            mean_epi_files = dir(fullfile(baseDir, imagingRawDir, participant_id, '*mean*.nii'));
            
            dest_dir = fullfile(baseDir, imagingDir, participant_id);
            if ~exist(dest_dir,'dir')
                mkdir(dest_dir)
            end
            
            % loop on runs of the session:
            for run = 1:length(realigned_epi_files)
                % realigned (and unwarped) images names:
                source = fullfile(realigned_epi_files(run).folder, realigned_epi_files(run).name);
                
                out_name = realigned_epi_files(run).name(length(prefix) + 1:end); % remove the prefix from realigned (and unwarped) file names
                dest = fullfile(dest_dir, out_name);
                % move to destination:
                fprintf('moving and renaming imaging_data_raw/%s --> imaging_data/%s\n',realigned_epi_files(run).name,out_name)
                [status,msg] = copyfile(source,dest);
                if ~status  
                    error('FUNC:move_realigned_images -> %s',msg)
                end
                
                % realign parameters names:
                source = fullfile(rp_files(run).folder, rp_files(run).name);
                dest = fullfile(dest_dir, rp_files(run).name);
                % move to destination:
                [status,msg] = copyfile(source,dest);
                if ~status  
                    error('FUNC:move_realigned_images -> %s',msg)
                end
            end
            
            % move the mean epi file:
            source = fullfile(mean_epi_files(1).folder, mean_epi_files(1).name);
            dest = fullfile(dest_dir, mean_epi_files(1).name);
            % move to destination:
            [status,msg] = copyfile(source, dest);
            if ~status  
                error('BIDS:move_realigned_images -> %s',msg)
            end
        
        case 'FUNC:meanimage_bias_correction_ss'
            % EPI images often contain smooth artifacts caused by MRI
            % physics which make the intensity of signal from the same
            % tissue (e.g., grey matter, white matter) non-uniform. This
            % step perform bias correction and creates an image where the
            % signal from each tissue type is more uniform. This image is
            % then co-registered to the anatomical image. Bias correction
            % help make co-registration more accurate. If the realignment
            % was done with respect to the first volume of each run of each
            % session, the mean image will be calculated on the first run
            % of each session and will be called 'meanu*_run_01.nii'
            % ('mean' indicates the image is average of the volumes and 'u'
            % indicates it's unwarped). Therefore, we do the bias
            % correction on this file. But if you do the realignment to the
            % mean epi of every run, the generated mean file will be named
            % 'umeanepi_*' and we do the bias correction on this file. In
            % addition, this step generates five tissue probability maps
            % (c1-5) for grey matter, white matter, csf, bone and soft
            % tissue.
            mean_epi = dir(fullfile(baseDir, imagingDir, participant_id, 'mean*.nii'));
            P{1} = fullfile(mean_epi.folder, mean_epi.name);
            spmj_bias_correct(P);
        case 'FUNC:coreg_ss'
            % coregister rbumean image to anatomical image for each session
            
            % (1) Manually seed the functional/anatomical registration
            % - Open fsleyes
            % - Add anatomical image and b*mean*.nii (bias corrected mean) image to overlay
            % - click on the bias corrected mean image in the 'Overlay
            %   list' in the bottom left of the fsleyes window.
            %   list to highlight it.
            % - Open tools -> Nudge
            % - Manually adjust b*mean*.nii image to the anatomical by 
            %   changing the 6 paramters (tranlation xyz and rotation xyz) 
            %   and Do not change the scales! 
            % - When done, click apply and close the tool tab. Then to save
            %   the changes, click on the save icon next to the mean image 
            %   name in the 'Overlay list' and save the new image by adding
            %   'r' in the beginning of the name: rb*mean*.nii. If you don't
            %   set the format to be .nii, fsleyes automatically saves it as
            %   a .nii.gz so either set it or gunzip afterwards to make it
            %   compatible with SPM.
            
            % (2) Run automated co-registration to register bias-corrected meanimage to anatomical image
            mean_epi = dir(fullfile(baseDir, imagingDir, participant_id, ['bmean' prefix, participant_id '_run_*.nii']));
            
            J.source = {fullfile(mean_epi.folder, mean_epi.name)}; 
            J.ref = {fullfile(baseDir, anatomicalDir, participant_id, [participant_id, '_T1w','.nii'])};
            J.other = {''};
            J.eoptions.cost_fun = 'nmi';
            J.eoptions.sep = [4 2];
            J.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            J.eoptions.fwhm = [7 7];
            matlabbatch{1}.spm.spatial.coreg.estimate=J;
            spm_jobman('run',matlabbatch);

        case 'FUNC:make_samealign_ss'
            % align to registered bias corrected mean image of each session
            % (rb*mean*.nii). Alignment happens only by changing the
            % transform matrix in the header files of the functional 4D
            % .nii files to the transform matrix that aligns them to
            % anatomical. The reason that it works is: 1) in the
            % realignment (+unwarping) process, we have registered every
            % single volume of every single run to the first volume of the
            % first run of the session (rtm=0). 2) In the same step, for each
            % session, a mean functional image (meanepi*.nii or meanu*.nii
            % based on the rtm option) was generated. This mean image is
            % alread in the space of all the functional volumes. Later we
            % coregister this image to the anatomical space. Therefore, if
            % we change the transformation matrices of all the functional
            % volumes to the transform matrix of the coregistered image,
            
            epi_files = dir(fullfile(baseDir, imagingDir, participant_id, [participant_id '_run_*.nii']));
            epi_list = {}; % Initialize as an empty cell array
            for run = 1:length(epi_files)
                epi_list{end+1} = fullfile(epi_files(run).folder, epi_files(run).name);
            end
            
            % select the reference image:
            mean_epi = dir(fullfile(baseDir, imagingDir, participant_id, ['bmean' prefix participant_id '_run_*.nii']));
            P{1} = fullfile(mean_epi.folder, mean_epi.name);
            
            % select images to be realigned:
            Q = {};
            for r = 1:length(epi_list)
                for j = 1:pinfo.numTR
                    Q{end+1} = fullfile(sprintf('%s,%d', epi_list{r}, j));
                end
            end
            
            spmj_makesamealign_nifti(char(P),char(Q));
        
        case 'FUNC:make_maskImage_one_ses'
            % Make mask images (noskull and gray_only) for 1st level glm
            % bias corrected mean epi image:
            mean_epi = dir(fullfile(baseDir, imagingDir, participant_id, ['bmean' prefix participant_id '_run_01.nii']));
            nam{1} = fullfile(mean_epi.folder, mean_epi.name);
            
            nam{2} = fullfile(baseDir, anatomicalDir, participant_id, ['c1', participant_id, '_T1w','.nii']);
            nam{3} = fullfile(baseDir, anatomicalDir, participant_id, ['c2', participant_id, '_T1w','.nii']);
            nam{4} = fullfile(baseDir, anatomicalDir, participant_id, ['c3', participant_id, '_T1w','.nii']);
            spm_imcalc(nam, fullfile(baseDir, imagingDir, participant_id, 'rmask_noskull.nii'), 'i1>1 & (i2+i3+i4)>0.2')
            
            source = fullfile(baseDir, imagingDir, participant_id, 'rmask_noskull.nii');
            dest = fullfile(baseDir, anatomicalDir, participant_id, 'rmask_noskull.nii');
            movefile(source, dest);
            
            % gray matter mask for covariance estimation
            % ------------------------------------------
            nam={};
            % nam{1}  = fullfile(imagingDir,subj_id{sn}, 'sess1', ['rb' prefix 'meanepi_' subj_id{sn} '.nii']);

            % IS THIS CHANGE CORRECT??
            % nam{1}  = fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)),sprintf('sess%d',sess), ['rb' prefix 'meanepi_' char(pinfo.subj_id(pinfo.sn==sn)) '.nii']);
            % bias corrected mean epi image:
            mean_epi = dir(fullfile(baseDir, imagingDir, participant_id, ['bmean' prefix participant_id '_run_*.nii']));
            nam{1} = fullfile(mean_epi(1).folder, mean_epi(1).name);

            nam{2} = fullfile(baseDir, anatomicalDir, participant_id, ['c1', participant_id, '_T1w','.nii']);
            spm_imcalc(nam, fullfile(baseDir, imagingDir, participant_id, 'rmask_gray.nii'), 'i1>1 & i2>0.4')
            
            source = fullfile(baseDir, imagingDir, participant_id, 'rmask_gray.nii');
            dest = fullfile(baseDir, anatomicalDir, participant_id, 'rmask_gray.nii');
            movefile(source,dest);

    end 

end
