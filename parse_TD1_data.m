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

% ========== get regressor info:
glm_info = []; % flat dataframe containing info of betas
run = [1:8]';
digit = [1:5,1:5]';
hand = [1,1,1,1,1,2,2,2,2,2]';

run_vec = repelem(run,length(digit));
digit_vec = repmat(digit,length(run),1);
hand_vec = repmat(hand,length(run),1);

glm_info.run = run_vec;
glm_info.digit = digit_vec;
glm_info.hand = hand_vec;

save(fullfile('analysis','hd1_glm_betas.mat'),'-struct','glm_betas')
dsave(fullfile('analysis','hd1_glm_info.tsv'),glm_info);


