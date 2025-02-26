function varargout = efcp_anat(what, varargin)
    if ismac
        baseDir = '/Volumes/Diedrichsen_data$/data/Chord_exp/EFC_patternfMRI';
        baseDir = '/Users/alighavampour/Desktop/Projects/EFC_patternfMRI/data';
    elseif isunix
        baseDir = '/cifs/diedrichsen/data/Chord_exp/EFC_patternfMRI';
    else
        disp('Running on Windows or another OS');
    end
    
    bidsDir = 'BIDS'; % Raw data post AutoBids conversion
    anatomicalDir = 'anatomicals'; % anatomical files (individual space)
    freesurferDir = 'surfaceFreesurfer'; % freesurfer reconall output
    surfacewbDir = 'surfaceWB'; % fs32k template 
    suitDir = 'suit'; % SUIT 2.0 outputs
    regDir = 'ROI';
    wbDir = 'surfaceWB';
    SPMhome = spm('dir');
    
    pinfo = dload(fullfile(baseDir,'participants.tsv'));

    switch(what)
        case 'BIDS:move_unzip_raw_anat'
            % Moves, unzips and renames raw anatomical from 
            % BIDS directory. After you run this function you will find
            % an anatomical Nifti file named <subj_id>_T1w.nii in the 
            % <project_id>/anatomicals/<subj_id>/ directory.
            % This function assumes that the anatomical is in the first
            % session in the BIDS dir
            
            % handling input args:
            sn = [];
            vararginoptions(varargin,{'sn'})
            if isempty(sn)
                error('BIDS:move_unzip_raw_anat -> ''sn'' must be passed to this function.')
            end

            % get participant row from participant.tsv
            subj_row = getrow(pinfo, pinfo.sn == sn);
            
            % get subj_id
            subj_id = subj_row.participant_id{1};

            % get anatomical session
            anat_ses = subj_row.anat_ses;

            % get the anatomical name
            anat_name = subj_row.anat_name{1};
            
            % anatomical file
            anat_full_path = fullfile(baseDir,bidsDir,sprintf('sub-%s',subj_id),sprintf('ses-%.2d',anat_ses),'anat',sprintf('sub-%s_ses-%.2d_%s.nii.gz', subj_id, anat_ses, anat_name));
            
            % Define output directory
            output_folder = fullfile(baseDir,anatomicalDir, subj_id);
            if ~exist(output_folder, 'dir')
               mkdir(output_folder)
            end
            output_file = fullfile(output_folder,sprintf('%s_T1w_raw.nii.gz', subj_id));
            
            % copy file to destination:
            copyfile(anat_full_path, output_file);
            
            % unzip the .gz files to make usable for SPM:
            gunzip(output_file);
            
            % delete the compressed file:
            delete(output_file);
        
        
        case 'ANAT:center_ac' 
            % recenter to AC (manually retrieve coordinates)
            % Before running this step you need to manually fill in the AC
            % coordinates in the participants.tsv file
            % run spm display to get the AC coordinates

            sn=[];
            vararginoptions(varargin,{'sn'})
            if isempty(sn)
                error('ANAT:center_ac -> ''sn'' must be passed to this function.')
            end
            
            % get subj row from participants.tsv
            subj_row=getrow(pinfo,pinfo.sn== sn);
            subj_id = subj_row.participant_id{1};

            % Get the anat of subject
            subj_anat_img = fullfile(baseDir,anatomicalDir, subj_id, sprintf('%s_T1w_LPI.nii', subj_id));

            % get location of ac
            locACx = subj_row.locACx;
            locACy = subj_row.locACy;
            locACz = subj_row.locACz;

            loc_AC = [locACx locACy locACz];
            loc_AC = loc_AC';

            % Recenter
            V               = spm_vol(subj_anat_img);
            dat             = spm_read_vols(V);
            % Solve the equation A @ loc_AC + Translation = 0
            A = V.mat(1:3,1:3);
            Trans = -A * loc_AC;
            V.mat(1:3,4) = Trans;

            % Modify filename
            new_filename = fullfile(baseDir,anatomicalDir, subj_id, sprintf('%s_T1w.nii', subj_id));
            V.fname = new_filename;
            spm_write_vol(V,dat);

        case 'ANAT:segment' 
            % segment the anatomical image
            % check results when done
            sn=[];
            vararginoptions(varargin,{'sn'})
            if isempty(sn)
                error('ANAT:segment -> ''sn'' must be passed to this function.')
            end
            
            subj_row = getrow(pinfo,pinfo.sn == sn);
            subj_id = subj_row.participant_id{1};
            
            subj_anat = fullfile(baseDir, anatomicalDir, subj_id, sprintf('%s_T1w.nii', subj_id));
            
            J.channel.vols     = {subj_anat};
            J.channel.biasreg  = 0.001;
            J.channel.biasfwhm = 60;
            J.channel.write    = [1 0];
            J.tissue(1).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,1')};
            J.tissue(1).ngaus  = 1;
            J.tissue(1).native = [1 0];
            J.tissue(1).warped = [0 0];
            J.tissue(2).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,2')};
            J.tissue(2).ngaus  = 1;
            J.tissue(2).native = [1 0];
            J.tissue(2).warped = [0 0];
            J.tissue(3).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,3')};
            J.tissue(3).ngaus  = 2;
            J.tissue(3).native = [1 0];
            J.tissue(3).warped = [0 0];
            J.tissue(4).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,4')};
            J.tissue(4).ngaus  = 3;
            J.tissue(4).native = [1 0];
            J.tissue(4).warped = [0 0];
            J.tissue(5).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,5')};
            J.tissue(5).ngaus  = 4;
            J.tissue(5).native = [1 0];
            J.tissue(5).warped = [0 0];
            J.tissue(6).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,6')};
            J.tissue(6).ngaus  = 2;
            J.tissue(6).native = [0 0];
            J.tissue(6).warped = [0 0];
    
            J.warp.mrf     = 1;
            J.warp.cleanup = 1;
            J.warp.reg     = [0 0.001 0.5 0.05 0.2];
            J.warp.affreg  = 'mni';
            J.warp.fwhm    = 0;
            J.warp.samp    = 3;
            J.warp.write   = [1 1];
            matlabbatch{1}.spm.spatial.preproc=J;
            spm_jobman('run',matlabbatch);

        case 'SURF:reconall' 
            % Freesurfer reconall routine
            % Calls recon-all, which performs, all of the
            % FreeSurfer cortical reconstruction process
            sn=[];
            vararginoptions(varargin,{'sn'})
            if isempty(sn)
                error('SURF:reconall -> ''sn'' must be passed to this function.')
            end
            
            subj_row = getrow(pinfo, pinfo.sn == sn);
            subj_id = subj_row.participant_id{1};   
        
            % recon all inputs
            fs_dir = fullfile(baseDir, freesurferDir);
            anatomical_dir = fullfile(baseDir, anatomicalDir);
            anatomical_name = sprintf('%s_T1w.nii', subj_id);
            
            % Get the directory of subjects anatomical;
            freesurfer_reconall(fs_dir, subj_id, ...
                fullfile(anatomical_dir, subj_id, anatomical_name));
            
        case 'SURF:fs2wb'          
            % Resampling subject from freesurfer fsaverage to fs_LR        
            res  = 32;          % resolution of the atlas. options are: 32, 164
            hemi = [1, 2];      % list of hemispheres
            
            sn=[];
            vararginoptions(varargin,{'sn'})
            if isempty(sn)
                error('SURF:fs2wb -> ''sn'' must be passed to this function.')
            end
    
            subj_row = getrow(pinfo,pinfo.sn == sn);
            subj_id = subj_row.participant_id{1};  
            
            % get the subject id folder name
            outDir = fullfile(baseDir, surfacewbDir); 
            fs_dir = fullfile(baseDir,freesurferDir);
            surf_resliceFS2WB(subj_id, fs_dir, outDir, 'hemisphere', hemi, 'resolution', sprintf('%dk', res))
        
        case 'ROI:define'    
            sn = [];
            glm = 1;
            ses = 1;
            atlas = 'ROI';
            
            vararginoptions(varargin,{'sn', 'atlas', 'glm', 'ses'});
            
            if isfolder('/Volumes/diedrichsen_data$/data/Atlas_templates/fs_LR_32')
                atlasDir = '/Volumes/diedrichsen_data$/data/Atlas_templates/fs_LR_32';
            elseif isfolder('/cifs/diedrichsen/data/Atlas_templates/fs_LR_32')
                atlasDir = '/cifs/diedrichsen/data/Atlas_templates/fs_LR_32';
            end
            atlasH = {sprintf('%s.32k.L.label.gii', atlas), sprintf('%s.32k.R.label.gii', atlas)};
            atlas_gii = {gifti(fullfile(atlasDir, atlasH{1})), gifti(fullfile(atlasDir, atlasH{1}))};
            
            % get participant row from participant.tsv
            participant_row = getrow(pinfo, pinfo.sn == sn);
            
            % get participant id:
            participant_id = participant_row.participant_id{1};

            Hem = {'L', 'R'};
            R = {};
            r = 1;
            for h = 1:length(Hem)
                for reg = 1:length(atlas_gii{h}.labels.name)
                    R{r}.white = fullfile(baseDir, wbDir, participant_id, [participant_id '.' Hem{h} '.white.32k.surf.gii']);
                    R{r}.pial = fullfile(baseDir, wbDir, participant_id, [participant_id '.' Hem{h} '.pial.32k.surf.gii']);
                    R{r}.image = fullfile(baseDir, sprintf('glm%d',glm), participant_id, sprintf('ses-%.2d',ses), 'mask.nii');
                    R{r}.linedef = [5 0 1];
                    key = atlas_gii{h}.labels.key(reg);
                    R{r}.location = find(atlas_gii{h}.cdata==key);
                    R{r}.hem = Hem{h};
                    R{r}.name = atlas_gii{h}.labels.name{reg};
                    R{r}.type = 'surf_nodes_wb';
                    r = r+1;
                end
            end

            R = region_calcregions(R, 'exclude', [2 3; 2 4; 2 5; 4 5; 8 9; 2 8;...
                11 12; 11 13; 11 14; 13 14; 17 18; 11 17], 'exclude_thres', .8);
            
            output_path = fullfile(baseDir, regDir, participant_id, sprintf('ses-%.2d',ses));
            if ~exist(output_path, 'dir')
                mkdir(output_path)
            end
            
            Vol = fullfile(baseDir, sprintf('glm%d',glm), participant_id, sprintf('ses-%.2d',ses), 'mask.nii');
            for r = 1:length(R)
                img = region_saveasimg(R{r}, Vol, 'name', fullfile(baseDir, regDir, participant_id, sprintf('ses-%.2d',ses), sprintf('%s.%s.%s_glm%d.nii', atlas, R{r}.hem, R{r}.name, glm)));
            end
            
            save(fullfile(output_path, sprintf('%s_%s_glm%d_region.mat', participant_id, atlas, glm)), 'R');
    end
end
    
