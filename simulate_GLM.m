% SET THESE VALUES:
% PARAMS:
hrf_params = [6 16 1 1 6 0 32];

% NUMBER OF VOLUMES:
num_scans = 549;

% EVENT ONSET TIMES:
ons = (3:5:num_scans)';

% DURATION OF EACH EVENT:
event_duration = 10;
dur = event_duration*ones(size(ons))/1000;

% ====================================== 
global defaults; 
if (isempty(defaults)) 
    spm_defaults;
end

% Empty job structure:
J = [];

% onets of the events:
J.U.ons = ons;  % Ali's note: This is defined like this for the sake of this 
                % simulation. SPM's definition of events is much more complete
                % and detailed. Check spm_rwls_run_fmri_spec() for the details
                % of how U is made and how it looks like in real SPM. 
J.U.dur = dur;
J.nscan = num_scans;

J.timing.units = 'secs';
J.timing.RT = 1;

% number of temporal bins in which the TR is divided,
% defines the discrtization of the HRF inside each TR
J.timing.fmri_t = 16;

% slice number that corresponds to that acquired halfway in
% each TR
J.timing.fmri_t0 = 1;

% Specify hrf parameters for convolution with
% regressors
J.bases.hrf.derivs = [0,0];
J.bases.hrf.params = hrf_params;  % positive and negative peak of HRF
defaults.stats.fmri.hrf = J.bases.hrf.params; 

% Specify the order of the Volterra series expansion 
% for modeling nonlinear interactions in the BOLD response
% *Example Usage*: Most analyses use 1, assuming a linear
% relationship between neural activity and the BOLD
% signal.
J.volt = 1;


% RUNNING THE SIM:
SPM = spm_rwls_sim(J);
X = spm_fMRI_design_sim(SPM);

%% Simulation functions:
function SPM = spm_rwls_sim(job)
    % Ali's note: this function tries to match the 
    % spm_rwls_run_fmri_specs() function.

    SPM = []; % empty SPM structure. This is supposed to be similar to 
              % the SPM.mat
  
    % Ali's note: The SPM structure has many many values. The values that 
    % I am adding here are needed for this simulation.
    SPM.nscan = job.nscan;

    %-Timing parameters & Basis functions
    %==========================================================================
    %-Repeat time
    %--------------------------------------------------------------------------
    SPM.xY.RT = job.timing.RT;
    
    %-Basis function parameters
    %--------------------------------------------------------------------------
    SPM.xBF.UNITS = job.timing.units;
    SPM.xBF.T     = job.timing.fmri_t;
    SPM.xBF.T0    = job.timing.fmri_t0;
    
    %-Basis functions
    %--------------------------------------------------------------------------
    bf = char(fieldnames(job.bases));
    if strcmp(bf,'hrf') % if hrf was params were defiend in job:
        if all(job.bases.hrf.derivs == [0 0])
            SPM.xBF.name = 'hrf';
        elseif all(job.bases.hrf.derivs == [1 0])
            SPM.xBF.name = 'hrf (with time derivative)';
        elseif all(job.bases.hrf.derivs == [1 1])
            SPM.xBF.name = 'hrf (with time and dispersion derivatives)';
        else
            error('Unknown HRF derivative choices.');
        end
        if (isfield(job.bases.hrf,'params'))            % Hand over params for hrf: JD 2/2013
            SPM.xBF.params = job.bases.hrf.params;
        end
    else
        error("this simulation will not work without defining hrf params\n")
    end
    
    %-Model interactions (Volterra)
    %--------------------------------------------------------------------------
    SPM.xBF.Volterra = job.volt;

    %-Events
    %----------------------------------------------------------------------
    % Ali's note: This is defined like this for the sake of this 
    % simulation. SPM's definition of events is much more complete
    % and detailed. Check spm_rwls_run_fmri_spec() for the details
    % of how U is made and how it looks like in real SPM. 
    SPM.U = job.U;
end

function X = spm_fMRI_design_sim(SPM)
    % Ali's note: this function tries to match the 
    % spm_fMRI_design() function.
    fMRI_T     = SPM.xBF.T;
    fMRI_T0    = SPM.xBF.T0;

    %-Time units, dt = time bin {secs}
    %--------------------------------------------------------------------------
    SPM.xBF.dt     = SPM.xY.RT/SPM.xBF.T;
    
    %-Get basis functions
    %--------------------------------------------------------------------------
    SPM.xBF        = spm_get_bf(SPM.xBF);
    
    %-Create convolved stimulus functions or inputs
    %======================================================================
    %-Get inputs, neuronal causes or stimulus functions U
    %----------------------------------------------------------------------
    U = spm_get_ons_sim(SPM);

    %-Convolve stimulus functions with basis functions
    %----------------------------------------------------------------------
    [X,~,~] = spm_Volterra_sim(U, SPM.xBF.bf, SPM.xBF.Volterra);

    %-Number of scans for this session
    %----------------------------------------------------------------------
    k = SPM.nscan(s);
    
    %-Resample regressors at acquisition times (32 bin offset)
    %----------------------------------------------------------------------
    if ~isempty(X)
        X = X((0:(k - 1))*fMRI_T + fMRI_T0 + 32,:);
    end
end

function U = spm_get_ons_sim(SPM)
    % Ali's note: This is written based on spm_get_ons() function.
    TR = 1;

    k     = SPM.nscan;
    T     = SPM.xBF.T;
    dt    = SPM.xBF.dt;
    
    U = SPM.U;

    %-Onsets
    %----------------------------------------------------------------------
    ons   = U.ons(:);
    
    %-Durations
    %----------------------------------------------------------------------
    dur   = U.dur(:);
    
    u     = ons.^0;
    
    %-Create stimulus functions (32 bin offset)
    %======================================================================
    ton       = round(ons*TR/dt) + 33;               % onsets
    tof       = round(dur*TR/dt) + ton + 1;          % offset
    sf        = sparse((k*T + 128),size(u,2));
    ton       = max(ton,1);
    tof       = max(tof,1);
    for j = 1:length(ton)
        if size(sf,1) > ton(j)
            sf(ton(j),:) = sf(ton(j),:) + u(j,:);
        end
        if size(sf,1) > tof(j)
            sf(tof(j),:) = sf(tof(j),:) - u(j,:);
        end
    end
    sf        = cumsum(sf);                         % integrate
    sf        = sf(1:(k*T + 32),:);                 % stimulus
    
    %-Place in ouputs structure U
    %----------------------------------------------------------------------
    U.dt   = dt;         % - time bin {seconds}
    U.u    = sf;         % - stimulus function matrix
end

function [X,Xname,Fc] = spm_Volterra_sim(U,bf,V)
    % Ali's note: This is writtent based on spm_Volterra()
    % V option is not supported in this simulation. 
    
    %-1st order terms
    %==========================================================================
    X     = [];
    for k = 1:size(U.u,2)
    for p = 1:size(bf,2)
        x = U.u(:,k);
        d = 1:length(x);
        x = conv(full(x),bf(:,p));
        x = x(d);
        X = [X x];
    end
    end
    
    % Ali's note: These outputs are there in the orignal function plus a
    % few other things but we don't need them here in this simulation.
    Fc = [];
    Xname = [];
end
