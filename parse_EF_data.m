% Ali Ghavampour, 2024
% alighavam79@gmail.com
% aghavamp@uwo.ca

% Making usable analysis files from the extensionFlexion study, from the
% Arbuckle et al. 2020 JNeuro paper. 

% Experiment:
% 9 healthy subjects
% Digits of right hand were pressed towards Extension and Flexion in three
% different force levels. For flexion: [1.5N, 2N, 2.5N] ; for extension: [1N, 1.5N, 2N]
% This gives 30 conditions = 5-digits * 2-directions * 3-forcelvl
% one fMRI session, 8 runs, each condition repeated twice -> 60 trials per
% run.

% fMRI analysis:
% SPM12, 248 regressors: (30-condition regressors + 1-intercept) * 8-runs

clear;
close all;
clc;
baseDir = '/Volumes/Diedrichsen_data$/data/FingerPattern/extensionFlexion/fmri';
regDir = fullfile(baseDir,'RegionOfInterest');
glmDir = fullfile(baseDir,'glm1');


glm_betas = []; % binary .mat dataframe

% ========== get betas:
% load glm1 betas:
glm1 = load(fullfile(regDir,'glm1_reg_betas.mat'));

% concatenate the glm values for all sn, region, and make a single matrix:
beta = [];
ResMs = [];
xyz = [];
sn = [];
region = [];
for i = 1:length(glm1.beta)
    beta = [beta; glm1.beta{i}'];
    ResMs = [ResMs; glm1.resMS{i}'];
    xyz = [xyz; glm1.xyz{i}];
    sn = [sn ; repelem(glm1.SN(i), size(glm1.xyz{i},1))'];
    region = [region ; repelem(glm1.region(i), size(glm1.xyz{i},1))'];
end

glm_betas.beta = beta(:,1:240); % last 8 regressors are run intercepts
glm_betas.ResMs = ResMs;
glm_betas.xyz = xyz;
glm_betas.sn = sn;
glm_betas.region = region;

% get and add region names to glm_beta:
glm_betas.region_name = strings(size(glm_betas.region));
glm_betas.hemi = zeros(size(glm_betas.region));

N = length(unique(sn));
for i = 1:N
    % load subj region info:
    tmp = load(fullfile(regDir,sprintf('regions_s%.2d.mat',i))).R;
    for j = 1:length(tmp)
        name = tmp{j}.name;
        roi = tmp{j}.roi;
        hemi = tmp{j}.hemi;
        glm_betas.region_name(glm_betas.region==roi,:) = string(name);
        glm_betas.hemi(glm_betas.region==roi) = hemi;
    end
end
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



