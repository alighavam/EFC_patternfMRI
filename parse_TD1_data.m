% Ali Ghavampour, 2024
% alighavam79@gmail.com
% aghavamp@uwo.ca

% Making usable analysis files from the TenDigit 1 study. This data was
% used in Ejaz 2015. 

% Experiment:
% 6 healthy subjects
% Digits of the right and left hand were flexed to 2.3N. In each run, a
% finger was probed for 3 trials. Each trial they pressed the finger 5
% times in a paced manner. 
% That gives 10 conditions = 5-digit * 2-hands
% one fMRI session, 8 runs, each condition repeated three times -> 30
% trials per run. 

% fMRI analysis:
% GLM had 10 regressors (1 per condition) per run. 
% So, 10-regressor * 8-runs = 80 overall.

clear;
close all;
clc;
baseDir = '/Volumes/Diedrichsen_data$/data/FingerPattern/tendigit1';
regDir = fullfile(baseDir,'RegionOfInterest');
glmDir = fullfile(baseDir,'GLM_firstlevel_run');


glm_betas = []; % binary .mat dataframe

% ========== get betas:
% load glm betas:
glm = load(fullfile(regDir,'reg_data_8.mat'));

% concatenate the glm values for all sn, region, and make a single matrix:
beta = [];
ResMs = [];
xyz = [];
sn = [];
region = [];

beta = glm.beta;
ResMs = glm.ResMs;
xyz = glm.xyz;
sn = glm.SN;
region = glm.regNum;

glm_betas.beta = beta; % last 8 regressors are run intercepts
glm_betas.ResMs = ResMs;
glm_betas.xyz = xyz;
glm_betas.sn = sn;
glm_betas.region = region;

% get and add region names to glm_beta:
% glm_betas.region_name = strings(size(glm_betas.region));
% glm_betas.hemi = zeros(size(glm_betas.region));

% N = length(unique(sn));
% for i = 1:N
%     % load subj region info:
%     tmp = load(fullfile(regDir,sprintf('s%.2d_regions_8.mat',i))).R;
%     for j = 1:length(tmp)
%         name = tmp{j}.name;
%         roi = tmp{j}.roi;
%         hemi = tmp{j}.hemi;
%         glm_betas.region_name(glm_betas.region==roi,:) = string(name);
%         glm_betas.hemi(glm_betas.region==roi) = hemi;
%     end
% end
% glm_betas = rmfield(glm_betas,'region_name');

% ========== get regressor info:
glm_info = []; % flat dataframe containing info of betas
for i = 1:N
    % load subj region info:
    sn_info = load(fullfile(glmDir,sprintf('s%.2d',i),'SPM_info.mat'));
    sn_info = rmfield(sn_info,'regtype');
    sn_info = rmfield(sn_info,'regname');
    glm_info = addstruct(glm_info,sn_info,'row',1);
end
sn = glm_info.SN;
glm_info.sn = sn;
glm_info = rmfield(glm_info,'SN');

save(fullfile('analysis','ef_glm_betas.mat'),'glm_betas')
dsave(fullfile('analysis','ef_glm_info.tsv'),glm_info);



