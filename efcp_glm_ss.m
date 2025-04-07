function varargout = efcp_glm_ss(what, varargin)

    % Use a different baseDir when using your local machine or the cbs
    % server. Add more directory if needed. Use single quotes ' and not
    % double quotes " because some spm function raise error with double
    % quotes
    if ismac
        baseDir = '/Users/alighavampour/Desktop/Projects/EFC_patternfMRI/data/efcp_fMRI';
    elseif isunix
        baseDir = '';
    else
        disp('Running on Windows or another OS');
    end
    
    sn = [];
    glm = [];
    type = 'spmT';
    atlas = 'ROI';
    derivs = [0, 0];
    hrf_params = [4.5 11 1 1 2 0 32];
    vararginoptions(varargin,{'sn', 'type', 'glm', 'hrf_params', 'atlas','derivs'})
    
    glmEstDir = 'glm';
    behavDir = 'behavioural';
    anatomicalDir = 'anatomicals';
    imagingDir = 'imaging_data';
    wbDir = 'surfaceWB';
    regDir = 'ROI';
    
    pinfo = dload(fullfile(baseDir,'participants.tsv'));

    % get participant row from participant.tsv
    participant_row = getrow(pinfo, pinfo.sn==sn);
    
    % get subj_id
    participant_id = participant_row.participant_id{1};
    
    % define subject runs:
    runs = spmj_dotstr2array(participant_row.glm_runs{1});
    switch what
        case 'GLM:make_glm1'
            
            dat_file = dir(fullfile(baseDir, behavDir, participant_id, 'efcp_*.dat'));
            D = dload(fullfile(dat_file.folder, dat_file.name));
            D = getrow(D, ismember(D.BN,runs));
            
            D.repetition = ones(length(D.TN), 1); % Initialize repetition column with 1

            for i = 1:length(D.TN)
                if i == 1
                    D.repetition(i) = 1;
                else
                    if D.chordID(i) == D.chordID(i-1)
                        D.repetition(i) = 2;
                    end
                end
            end
    
            events.BN = [];
            events.TN = [];
            events.onset = [];
            events.duration = [];
            events.eventtype = {};
            events.chordID = [];
            events.repetition = [];
            
            for rep = unique(D.repetition)'
                for chordID = unique(D.chordID)'
                    rows = (D.chordID == chordID) & (D.repetition == rep);
                    events.BN = [events.BN; D.BN(rows)];
                    events.TN = [events.TN; D.TN(rows)];
                    events.onset = [events.onset; D.startTimeReal(rows)];
                    events.duration = [events.duration; repmat(10, [sum(rows),1])];
                    events.repetition = [events.repetition; D.repetition(rows)];
                    events.chordID = [events.chordID; D.chordID(rows)];
                    events.eventtype = [events.eventtype; repmat({sprintf('chordID:%d,repetition:%d', chordID, rep)}, [sum(rows), 1])];
                end
            end
            events = struct2table(events);
            events.onset = events.onset / 1000;
            events.duration = events.duration / 1000;
            
            varargout{1} = events;

        case 'GLM:make_event'
            operation  = sprintf('GLM:make_glm%d', glm);
            
            events = efcp_glm_ss(operation, 'sn', sn);
            
            % export:
            output_folder = fullfile(baseDir, behavDir, participant_id);
            writetable(events, fullfile(output_folder, sprintf('glm%d_events.tsv', glm)), 'FileType', 'text', 'Delimiter','\t')
            
        case 'GLM:design'
            % Import globals from spm_defaults 
            global defaults; 
            if (isempty(defaults))
                spm_defaults;
            end
            
            currentDir = pwd;
            
            if isempty(sn)
                error('GLM:design -> ''sn'' must be passed to this function.')
            end
            
            if isempty(glm)
                error('GLM:design -> ''glm'' must be passed to this function.')
            end

            % Load data once, outside of session loop
            % D = dload(fullfile(baseDir,behavDir,subj_id, sprintf('smp2_%d.dat', sn)));
            events_file = sprintf('glm%d_events.tsv', glm);
            
            Dd = dload(fullfile(baseDir, behavDir, participant_id, events_file));
            
            regressors = unique(Dd.eventtype);
            nRegr = length(regressors);
            
            % init J
            J = [];
            T = [];
            J.dir = {fullfile(baseDir, sprintf('glm%d', glm), participant_id)};
            J.timing.units = 'secs';
            J.timing.RT = 1;
            
            % number of temporal bins in which the TR is divided,
            % defines the discrtization of the HRF inside each TR
            J.timing.fmri_t = 16;
            
            % slice number that corresponds to that acquired halfway in
            % each TR
            J.timing.fmri_t0 = 1;
            
            epi_files = dir(fullfile(baseDir, imagingDir, participant_id, [participant_id '_run_*.nii']));
            for run = runs
                % Setup scans for current session
                J.sess(run).scans = {fullfile(epi_files(run).folder, epi_files(run).name)};
                
                % Preallocate memory for conditions
                J.sess(run).cond = repmat(struct('name', '', 'onset', [], 'duration', []), nRegr, 1);
                
                for regr = 1:nRegr
                    % cue = Dd.cue(regr);
                    % stimFinger = Dd.stimFinger(regr);
                    rows = find(Dd.BN == run & strcmp(Dd.eventtype, regressors(regr)));
                    % cue_id = unique(Dd.cue_id(rows));
                    % stimFinger_id = unique(Dd.stimFinger_id(rows));
                    % epoch = unique(Dd.epoch(rows));
                    % instr = unique(Dd.instruction(rows));
                    
                    % Regressor name
                    J.sess(run).cond(regr).name = regressors{regr};
                    
                    % Define durationDuration(regr));
                    J.sess(run).cond(regr).duration = Dd.duration(rows); % needs to be in seconds
                    
                    % Define onset
                    J.sess(run).cond(regr).onset  = Dd.onset(rows);
                    
                    % Define time modulator
                    % Add a regressor that account for modulation of
                    % betas over time
                    J.sess(run).cond(regr).tmod = 0;
                    
                    % Orthogonalize parametric modulator
                    % Make the parametric modulator orthogonal to the
                    % main regressor
                    J.sess(run).cond(regr).orth = 0;
                    
                    % Define parametric modulators
                    % Add a parametric modulators, like force or
                    % reaction time. 
                    J.sess(run).cond(regr).pmod = struct('name', {}, 'param', {}, 'poly', {});

                    %
                    % filling in "reginfo"
                    TT.sn        = sn;
                    TT.run       = run;
                    TT.name      = regressors(regr);
                    % TT.cue       = cue_id;
                    % TT.epoch     = epoch;
                    % TT.stimFinger = stimFinger_id;
                    % TT.instr = instr;       
                    
                    T = addstruct(T, TT);
                end

                % Specify high pass filter
                J.sess(run).hpf = 64;
                
                % J.sess(run).multi
                % Purpose: Specifies multiple conditions for a session. Usage: It is used
                % to point to a file (.mat or .txt) that contains multiple conditions,
                % their onsets, durations, and names in a structured format. If you have a
                % complex design where specifying conditions manually within the script is
                % cumbersome, you can prepare this information in advance and just
                % reference the file here. Example Setting: J.sess(run).multi =
                % {'path/to/multiple_conditions_file.mat'}; If set to {' '}, it indicates
                % that you are not using an external file to specify multiple conditions,
                % and you will define conditions directly in the script (as seen with
                % J.sess(run).cond).
                J.sess(run).multi     = {''};                        

                % J.sess(run).regress
                % Purpose: Allows you to specify additional regressors that are not
                % explicitly modeled as part of the experimental design but may account for
                % observed variations in the BOLD signal. Usage: This could include
                % physiological measurements (like heart rate or respiration) or other
                % variables of interest. Each regressor has a name and a vector of values
                % corresponding to each scan/time point.
                J.sess(run).regress   = struct('name', {}, 'val', {});                        

                % J.sess(run).multi_reg Purpose: Specifies a file containing multiple
                % regressors that will be included in the model as covariates. Usage: This
                % is often used for motion correction, where the motion parameters
                % estimated during preprocessing are included as regressors to account for
                % motion-related artifacts in the BOLD signal. Example Setting:
                % J.sess(run).multi_reg = {'path/to/motion_parameters.txt'}; The file
                % should contain a matrix with as many columns as there are regressors and
                % as many rows as there are scans/time points. Each column represents a
                % different regressor (e.g., the six motion parameters from realignment),
                % and each row corresponds to the value of those regressors at each scan.
                J.sess(run).multi_reg = {''};
                
                % Specify factorial design
                J.fact             = struct('name', {}, 'levels', {});

                % Specify hrf parameters for convolution with
                % regressors
                J.bases.hrf.derivs = derivs;
                J.bases.hrf.params = hrf_params;  % positive and negative peak of HRF - set to [] if running wls (?)
                defaults.stats.fmri.hrf=J.bases.hrf.params; 
                
                % Specify the order of the Volterra series expansion 
                % for modeling nonlinear interactions in the BOLD response
                % *Example Usage*: Most analyses use 1, assuming a linear
                % relationship between neural activity and the BOLD
                % signal.
                J.volt = 1;

                % Specifies the method for global normalization, which
                % is a step to account for global differences in signal
                % intensity across the entire brain or between scans.
                J.global = 'None';

                % remove voxels involving non-neural tissue (e.g., skull)
                J.mask = {fullfile(baseDir, anatomicalDir, participant_id, 'rmask_noskull.nii')};
                
                % Set threshold for brightness threshold for masking 
                % If supplying explicit mask, set to 0  (default is 0.8)
                J.mthresh = 0.;

                % Create map where non-sphericity correction must be
                % applied
                J.cvi_mask = {fullfile(baseDir, anatomicalDir, participant_id,  'rmask_gray.nii')};

                % Method for non sphericity correction
                J.cvi = 'fast';
            end
            
            % remove empty rows (e.g., when skipping runs)
            J.sess = J.sess(~arrayfun(@(x) all(structfun(@isempty, x)), J.sess));
            
            if ~exist(J.dir{1},"dir")
                mkdir(J.dir{1});
            end
            
            dsave(fullfile(J.dir{1},'reginfo.tsv'), T);
            spm_rwls_run_fmri_spec(J);
            
            cd(currentDir)
            
            currentDir = pwd;
            
        case 'GLM:estimate' % estimate beta values

            currentDir = pwd;

            if isempty(sn)
                error('GLM:estimate -> ''sn'' must be passed to this function.')
            end

            if isempty(glm)
                error('GLM:estimate -> ''glm'' must be passed to this function.')
            end

%             fprintf('- Doing glm%d estimation for subj %s\n', glm, day_id, subj_id);
            subj_est_dir = fullfile(baseDir, sprintf('glm%d', glm), participant_id);                
            SPM = load(fullfile(subj_est_dir,'SPM.mat'));
            SPM.SPM.swd = subj_est_dir;

            iB = SPM.SPM.xX.iB;

            save(fullfile(subj_est_dir, "iB.mat"), "iB");

            spm_rwls_spm(SPM.SPM);

            cd(currentDir)
        
        case 'GLM:change_SPM.mat_format'
            % This is kinda crappy but, rsatoolbox has a SPM object which
            % makes loading and working with SPM.mat easy. But this object
            % uses scipy.io.loadmat() which does not work with v7.3
            % structure formats. Therefore, I will change v7.3 with matlab
            % to v7 here.
            subj_est_dir = fullfile(baseDir, sprintf('glm%d', glm), participant_id);                
            SPM = load(fullfile(subj_est_dir,'SPM.mat'));
            save(fullfile(subj_est_dir,'SPM_v7.mat'), '-struct', 'SPM', '-v7');
            
        case 'GLM:T_contrasts'
            
            currentDir = pwd;
            
            replace_xCon   = true;

            if isempty(sn)
                error('GLM:T_contrasts -> ''sn'' must be passed to this function.')
            end

            if isempty(glm)
                error('GLM:T_contrasts -> ''glm'' must be passed to this function.')
            end

            % get the subject id folder name
            fprintf('Contrasts for participant %s\n', participant_id)
            glm_dir = fullfile(baseDir, sprintf('glm%d', glm), participant_id);

            % load the SPM.mat file
            SPM = load(fullfile(glm_dir, 'SPM.mat')); SPM=SPM.SPM;
            
            if replace_xCon
                SPM  = rmfield(SPM,'xCon');
            end
            
            T    = dload(fullfile(glm_dir, 'reginfo.tsv'));
            T.name = cellstr(string(T.name));
            contrasts = unique(T.name);
            
            % Per condition contrasts:
            for c = 1:length(contrasts)
                contrast_name = contrasts{c};
                xcon = zeros(size(SPM.xX.X,2), 1);
                xcon(strcmp(T.name, contrast_name)) = 1;
                xcon = xcon / sum(xcon);
                if ~isfield(SPM, 'xCon')
                    SPM.xCon = spm_FcUtil('Set', contrast_name, 'T', 'c', xcon, SPM.xX.xKXs);
                    cname_idx = 1;
                elseif sum(strcmp(contrast_name, {SPM.xCon.name})) > 0
                    idx = find(strcmp(contrast_name, {SPM.xCon.name}));
                    SPM.xCon(idx) = spm_FcUtil('Set', contrast_name, 'T', 'c', xcon, SPM.xX.xKXs);
                    cname_idx = idx;
                else
                    SPM.xCon(end+1) = spm_FcUtil('Set', contrast_name, 'T', 'c', xcon, SPM.xX.xKXs);
                    cname_idx = length(SPM.xCon);
                end
                SPM = spm_contrasts(SPM,1:length(SPM.xCon));
                % SPM = rmfield(SPM,'xVi'); % 'xVi' take up a lot of space and slows down code!
                % save(fullfile(glm_dir, 'SPM_light.mat'), 'SPM')
                
                % rename contrast images and spmT images
                conName = {'con','spmT'};
                for n = 1:numel(conName)
                    oldName = fullfile(glm_dir, sprintf('%s_%2.4d.nii',conName{n},cname_idx));
                    newName = fullfile(glm_dir, sprintf('%s_%s.nii',conName{n},SPM.xCon(cname_idx).name));
                    movefile(oldName, newName);
                end
            end
            
            % Custom contrasts:
            % contrast based on num fingers:
            % contrasts = {'lhand','rhand','bi'};
            % for c = 1:length(contrasts)
            %     contrast_name = contrasts{c};
            %     xcon = zeros(size(SPM.xX.X,2), 1);
            %     xcon(contains(T.name, contrast_name)) = 1;
            %     xcon = xcon / sum(xcon);
            %     if ~isfield(SPM, 'xCon')
            %         SPM.xCon = spm_FcUtil('Set', contrast_name, 'T', 'c', xcon, SPM.xX.xKXs);
            %         cname_idx = 1;
            %     elseif sum(strcmp(contrast_name, {SPM.xCon.name})) > 0
            %         idx = find(strcmp(contrast_name, {SPM.xCon.name}));
            %         SPM.xCon(idx) = spm_FcUtil('Set', contrast_name, 'T', 'c', xcon, SPM.xX.xKXs);
            %         cname_idx = idx;
            %     else
            %         SPM.xCon(end+1) = spm_FcUtil('Set', contrast_name, 'T', 'c', xcon, SPM.xX.xKXs);
            %         cname_idx = length(SPM.xCon);
            %     end
            %     SPM = spm_contrasts(SPM,1:length(SPM.xCon));
            %     % SPM = rmfield(SPM,'xVi'); % 'xVi' take up a lot of space and slows down code!
            %     % save(fullfile(glm_dir, 'SPM_light.mat'), 'SPM')
            % 
            %     % rename contrast images and spmT images
            %     conName = {'con','spmT'};
            %     for n = 1:numel(conName)
            %         oldName = fullfile(glm_dir, sprintf('%s_%2.4d.nii',conName{n},cname_idx));
            %         newName = fullfile(glm_dir, sprintf('%s_%s.nii',conName{n},SPM.xCon(cname_idx).name));
            %         movefile(oldName, newName);
            %     end
            % end

            save('SPM.mat', 'SPM', '-v7.3');

            cd(currentDir)
            
        case 'GLM:all'
            spm_get_defaults('cmdline', true);  % Suppress GUI prompts, no request for overwirte
            
            % Check for and delete existing SPM.mat file
            % spm_file = fullfile(baseDir, [glmEstDir num2str(glm)], ['subj' num2str(sn)], 'SPM.mat');
            spm_file = fullfile(baseDir, [glmEstDir num2str(glm)], participant_id, 'SPM.mat');
            if exist(spm_file, 'file')
                delete(spm_file);
            end
            
            efcp_glm_ss('GLM:make_event', 'sn', sn, 'glm', glm)
            efcp_glm_ss('GLM:design', 'sn', sn, 'glm', glm, 'hrf_params', hrf_params, 'derivs', [0,0])
            efcp_glm_ss('GLM:estimate', 'sn', sn, 'glm', glm)
            efcp_glm_ss('GLM:T_contrasts', 'sn', sn, 'glm', glm)
            efcp_glm_ss('SURF:vol2surf', 'sn', sn, 'glm', glm, 'type', 'spmT')
            efcp_anat('ROI:define_ss', 'sn', sn, 'glm', glm)
            efcp_glm_ss('HRF:ROI_hrf_get', 'sn', sn, 'glm', glm)

        case 'SURF:vol2surf'
            currentDir = pwd;
            
            glmEstDir = [glmEstDir num2str(glm)];
            
            V = {};
            cols = {};
            if strcmp(type, 'spmT')
%                 filename = ['spmT_' id '.func.gii'];
                files = dir(fullfile(baseDir, glmEstDir, participant_id, 'spmT_*.nii'));
                for f = 1:length(files)
                    fprintf([files(f).name '\n'])
                    V{f} = fullfile(files(f).folder, files(f).name);
                    cols{f} = files(f).name;
                end
            elseif strcmp(type, 'beta')
                SPM = load(fullfile(baseDir, glmEstDir, participant_id, 'SPM.mat')); SPM=SPM.SPM;
                files = dir(fullfile(baseDir, glmEstDir, participant_id, 'beta_*.nii'));
                files = files(SPM.xX.iC);
                for f = 1:length(files)
                    fprintf([files(f).name '\n'])
                    V{f} = fullfile(files(f).folder, files(f).name);
                    cols{f} = files(f).name;
                end
            elseif strcmp(type, 'psc')
                files = dir(fullfile(baseDir, glmEstDir, participant_id, 'psc_*.nii'));
                for f = 1:length(files)
                    fprintf([files(f).name '\n'])
                    V{f} = fullfile(files(f).folder, files(f).name);
                    cols{f} = files(f).name;
                end
            elseif strcmp(type, 'con')
                files = dir(fullfile(baseDir, glmEstDir, participant_id, 'con_*.nii'));
                for f = 1:length(files)
                    fprintf([files(f).name '\n'])
                    V{f} = fullfile(files(f).folder, files(f).name);
                    cols{f} = files(f).name;
                end
            elseif strcmp(type, 'res')
                V{1} = fullfile(baseDir, glmEstDir, participant_id, 'ResMS.nii');
                cols{1} = 'ResMS';
            end

            hemLpial = fullfile(baseDir, wbDir, participant_id,  [participant_id '.L.pial.32k.surf.gii']);
            hemRpial = fullfile(baseDir, wbDir, participant_id, [participant_id '.R.pial.32k.surf.gii']);
            hemLwhite = fullfile(baseDir, wbDir, participant_id, [participant_id '.L.white.32k.surf.gii']);
            hemRwhite = fullfile(baseDir, wbDir, participant_id, [participant_id '.R.white.32k.surf.gii']);
            
            hemLpial = gifti(hemLpial);
            hemRpial = gifti(hemRpial);
            hemLwhite = gifti(hemLwhite);
            hemRwhite = gifti(hemRwhite);

            c1L = hemLpial.vertices;
            c2L = hemLwhite.vertices;
            c1R = hemRpial.vertices;
            c2R = hemRwhite.vertices;

            GL = surf_vol2surf(c1L,c2L,V,'anatomicalStruct','CortexLeft', 'exclude_thres', 0.9, 'faces', hemLpial.faces);
            GL = surf_makeFuncGifti(GL.cdata,'anatomicalStruct', 'CortexLeft', 'columnNames', cols);
    
            save(GL, fullfile(baseDir, wbDir, participant_id, [glmEstDir '.'  type '.L.func.gii']))
    
            GR = surf_vol2surf(c1R,c2R,V,'anatomicalStruct','CortexRight', 'exclude_thres', 0.9, 'faces', hemRpial.faces);
            GR = surf_makeFuncGifti(GR.cdata,'anatomicalStruct', 'CortexRight', 'columnNames', cols);

            save(GR, fullfile(baseDir, wbDir, participant_id, [glmEstDir '.' type '.R.func.gii']))
            
            cd(currentDir)
            
        case 'HRF:ROI_hrf_get' % Extract raw and estimated time series from ROIs
            currentDir = pwd;
            
            glmDir = fullfile(baseDir, [glmEstDir num2str(glm)]);
            time_series = [];
            
            fprintf('Extracting region time series for participant %s...\n', participant_id);
            
            % load SPM.mat
            cd(fullfile(glmDir, participant_id));
            SPM = load('SPM.mat'); SPM=SPM.SPM;
            
            % load ROI definition (R)
            R = load(fullfile(baseDir, regDir, participant_id, sprintf('%s_%s_glm%d_region.mat', participant_id, atlas, glm))); R=R.R;
            
            % extract time series data
            [y_raw, y_adj, y_hat, y_res, B] = region_getts(SPM,R,'stats','mean');
            time_series.y_raw = y_raw;
            time_series.y_adj = y_adj;
            time_series.y_hat = y_hat;
            time_series.y_res = y_res;
            save(fullfile(baseDir, regDir, participant_id, sprintf('time_series_glm%d.mat', glm)),'-struct','time_series','-v7');
            cd(currentDir)
            
            % Get onset structure, cut-out the trials of choice, and plot evoked
            % response
            figure;
            D = spmj_get_ons_struct(SPM);
            r = 2;
            pre = 8;
            post = 21;
            for i=1:size(D.block,1)
                D.y_adj(i,:)=cut(time_series.y_adj(:,2),pre,round(D.ons(i)),post,'padding','nan')';
                D.y_hat(i,:)=cut(time_series.y_hat(:,2),pre,round(D.ons(i)),post,'padding','nan')';
                D.y_res(i,:)=cut(time_series.y_res(:,2),pre,round(D.ons(i)),post,'padding','nan')';
            end
            
            T = D; %getrow(D,mod(D.num,2)==0); % Get the first onset for each double 
            traceplot([-pre:post],T.y_adj,'errorfcn','stderr'); % ,
            hold on;
            traceplot([-pre:post],T.y_hat,'linestyle',':',...
                    'linewidth',3); %
            drawline([-7,0,7,14,21],'dir','vert','linestyle',':');
            
            drawline([0],'dir','vert','linestyle','--');
            drawline([0],'dir','horz','linestyle','-');
            hold off;
            xlabel('TR');
            ylabel('activation');
    end
end

