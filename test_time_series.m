sn = 7;
T = load(['/Users/alighavampour/Desktop/Projects/bimanual_wrist/data/fMRI/ROI/s' num2str(sn,'%.2d') '/time_series_glm1.mat']);
spm_file = load(['/Users/alighavampour/Desktop/Projects/bimanual_wrist/data/fMRI/glm1/s' num2str(sn,'%.2d') '/SPM.mat']);
SPM = spm_file.SPM;

%%
runs = 1:size(T.y_raw,1)/740;
roi = 2;

for run = runs
    x = 1:740;
    y_raw = T.y_raw(SPM.Sess(run).row,roi);
    y_adj = T.y_adj(SPM.Sess(run).row,roi);
    y_hat = T.y_hat(SPM.Sess(run).row,roi);
    y_res = T.y_res(SPM.Sess(run).row,roi);

    figure;
    % plot(x, y_raw, 'k', 'linewidth', 1.5); 
    % plot(x, y_res, 'k', 'linewidth', 1.5); 
    hold on;
    plot(x, y_adj, 'k', 'linewidth', 1.5); 
    plot(x, y_hat, 'r', 'linewidth', 1.5);
    title('run %d',run)
end

%% Fit hrf params - redo GLM:
% baseDir = '/Users/alighavam/Desktop/Projects/bimanual_wrist/data/fMRI';
% baseDir = '/Users/ali/Desktop/Projects/bimanual_wrist/data';
baseDir = '/Users/alighavampour/Desktop/Projects/EFC_patternfMRI/data/efcp_fMRI';

sn = 101;
glm = 1;

% spm_file = load(['/Users/alighavampour/Desktop/Projects/bimanual_wrist/data/fMRI/glm' num2str(glm) '/s' num2str(sn,'%.2d') '/SPM.mat']);
spm_file = load(fullfile(baseDir,['glm' num2str(glm)], ['s' num2str(sn,'%.2d')], 'SPM.mat'));
SPM = spm_file.SPM;
% load ROI definition (R)
R = load(fullfile(baseDir, 'ROI', ['s' num2str(sn,'%.2d')], sprintf('s%s_ROI_glm%d_region.mat', num2str(sn,'%.2d'), glm))); 
R=R.R;

region_data = region_getdata(SPM.xY.VY,R);
% region_data = load('/Users/ali/Desktop/Projects/bimanual_wrist/data/region_data_s04.mat'); region_data=region_data.region_data;

%%

r = 2;
% [4.5 11 1 1 1 0 32]
hrf_params = [5 8 1 1 1.1 0 32];
pre = 8;
post = 22;
run = 1;

Yraw = region_data{r}; % Pick out data as a T x P matrix 

% Get the hemodynamic response in micro-time resolution
SPM.xBF.UNITS    = 'secs'; % units of the hrf
SPM.xBF.T        = 16; % microtime resolution: number of time samples per scan
SPM.xBF.dt       = SPM.xY.RT/SPM.xBF.T;  % Delta-t per sample 
SPM.xBF.T0       = 1; % first time bin
SPM.xBF.name     = 'fitted_hrf';
SPM.xBF.order    = 1;
SPM.xBF.Volterra = 1;  % volterra expansion order?
SPM.xBF.bf = spm_hrf(SPM.xBF.dt, hrf_params);
SPM.xBF.length = size(SPM.xBF.bf,1)*SPM.xBF.dt; % support in seconds 

% Reconvolve the design matrix with new HRF
SPM = spmj_glm_convolve(SPM);

% Restimate the betas and get predicted and residual response
[beta, Yhat, Yres] = spmj_glm_fit(SPM,Yraw);

% Diagnostic plots
figure;
t=[SPM.xBF.dt*SPM.xBF.T0:SPM.xBF.dt:SPM.xBF.length]; 
plot(t,SPM.xBF.bf);
drawline([0],'dir','horz','linestyle','-');
title(sprintf('params = %s',num2str(hrf_params)));

% Run, convolved design matrix, sum of regressors of interest
figure;
X=SPM.xX.X(SPM.Sess(run).row,SPM.Sess(run).col);
plot(sum(X,2));
title(sprintf('Run %d - overall response', run));

figure;
Yhat_run1 = mean(Yhat(SPM.Sess(run).row,:),2);
Yres_run1 = mean(Yres(SPM.Sess(run).row,:),2);
Yadj_run1 = Yres_run1+Yres_run1; 
t= SPM.Sess(run).row;
hold on;
plot(t, Yhat_run1, 'r', 'LineWidth', 1.5)
plot(t, Yadj_run1, 'k', 'LineWidth', 1)
title(sprintf('Run %d - overall response', run));

% Get onset structure, cut-out the trials of choice, and plot evoked
% response
figure;
D = spmj_get_ons_struct(SPM);
Yadj = Yres+Yhat; 
for i=1:size(D.block,1)
    D.y_adj(i,:)=cut(mean(Yadj,2),pre,round(D.ons(i)),post,'padding','nan')';
    D.y_hat(i,:)=cut(mean(Yhat,2),pre,round(D.ons(i)),post,'padding','nan')';
    D.y_res(i,:)=cut(mean(Yres,2),pre,round(D.ons(i)),post,'padding','nan')';
end

T = getrow(D,contains(D.eventname,'repetition:1')); % Get the first onset for each double 
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
















