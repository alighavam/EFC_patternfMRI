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


glm_beta = []; % binary .mat file
glm_info = []; % flat dataframe containing info of betas

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

glm_beta.beta = beta;
glm_beta.ResMs = ResMs;
glm_beta.xyz = xyz;
glm_beta.sn = sn;
glm_beta.region = region;

% get and add region names to glm_beta:
glm_beta.region_name = strings(size(glm_beta.region));
glm_beta.hemi = zeros(size(glm_beta.region));

N = length(unique(sn));
for i = 1:N
    % load subj region info:
    tmp = load(fullfile(regDir,sprintf('regions_s%.2d.mat',i))).R;
    for j = 1:length(tmp)
        name = tmp{j}.name;
        roi = tmp{j}.roi;
        hemi = tmp{j}.hemi;
        glm_beta.region_name(glm_beta.region==roi,:) = string(name);
        glm_beta.hemi(glm_beta.region==roi) = hemi;
    end
end







