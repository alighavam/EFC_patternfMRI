%% Tendigit1

close all;
clc;

subj_name={'s01','s02','s03','s04','s05','s06'};
run={'01','02','03','04','05','06','07','08'};
hem={'lh','rh'};
regname={'S1','M1','PMd','PMv','SMA','V12','AIP','OPJ'};
regType=[1:8 1:8]';
regSide=[1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2]';

TD = load('analysis/TD1/reg_data_run.mat');

