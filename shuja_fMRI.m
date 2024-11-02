% fMRI = load('/Volumes/Diedrichsen_data$/data/Chord_exp/ChordPatternDist/analysis/reg_data_1.mat');
% info = load('/Volumes/Diedrichsen_data$/data/Chord_exp/ChordPatternDist/analysis/SPM_info_1.mat');
fMRI = load('/Users/ali/Desktop/Projects/EFC_patternfMRI/reg_data_1.mat');
info = load('/Users/ali/Desktop/Projects/EFC_patternfMRI/SPM_info_1.mat');


% assuming that info is same across participants:
info = getrow(info,info.sn==1);

% info about the fMRI data:
hem         = {'lh','rh'};
hemName     = {'LeftHem','RightHem'};
regname     = {'S1','M1','PMd','PMv','SMA','V12','SPLa','SPLp'};
numregions  = 8;
regSide     = [1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2];
regType     = [1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8];

% extract beta + univariate prewhitening:
beta = fMRI.beta ./ sqrt(fMRI.ResMs);

% plot the brain:
idx = fMRI.sn==1;
idx_m1 = fMRI.sn==1 & fMRI.regNum==2;
figure; scatter3(fMRI.xyz(idx,1),fMRI.xyz(idx,2),fMRI.xyz(idx,3),5,'k','filled'); hold on;
scatter3(fMRI.xyz(idx_m1,1),fMRI.xyz(idx_m1,2),fMRI.xyz(idx_m1,3),5,'r','filled');

% Get the avg across participants, sess, runs:
sn = unique(fMRI.sn);
sess = unique(info.sess);
run = unique(info.run);
mean_M1 = zeros(length(sn),31);

% loop on subj:
for i = 1:length(sn)
    % loop on sess:
    for j = 1:length(sess)
        for k = 1:length(run)
            row = fMRI.regNum==2 & fMRI.sn==sn(i);
            col = info.sess==sess(j) & info.run==run(k);
            tmp_beta = mean(beta(row,col),1); 
            mean_M1(i,:) = mean_M1(i,:) + tmp_beta / length(run)/length(sess);
        end
    end
end

figure;
mean_subj = mean(mean_M1,1);
sem_subj = std(mean_M1,[],1)/sqrt(8);
plot(1:31,mean_subj,'k','linewidth',2); hold on;
errorbar(1:31,mean_subj,sem_subj,'LineStyle','none','color','r');
xlabel('chord')
ylabel('mean beta in M1')
xlim([0,32])


% figure;
% imagesc(beta')
% colormap('viridis')
% clim([0, 4])
% xlabel('voxels')
% ylabel('betas')

