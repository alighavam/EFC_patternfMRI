from utils.config import *
import pandas as pd
import numpy as np
from utils import dataset as ds
from utils import emg_handler
from utils import measures
import os
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy import signal

def trial_routine(row):
    '''
    This function performs all the necessary preprocessing for a single trial
    '''
    C = pd.DataFrame()

    # add row to the dataframe:
    C = pd.concat([C, row], ignore_index=True)
    C = C.rename(columns={'trialCorr': 'trial_correct'})

    return C

def subject_routine(subject, fs_force=500, fs_emg=2148.1481, bpf=[20, 500], lpf=30, debug=0):
    """
    This function is used to preprocess the data of a subject
    
    params:
        subject: int, the subject number
        fs_force: the sampling frequency of the force data
        fs_emg: the sampling frequency of the emg data
    """
    # empty dataframe to store the data:
    D = pd.DataFrame()
    df_mov = pd.DataFrame(columns=['sn', 'BN', 'TN', 'trial_correct', 'state', 'time','f1','f2','f3','f4','f5'])
    
    # Load the .dat file:
    dat_file_name = os.path.join(DATA_PATH, 'behavioural', f'subj{subject}', f'efcp_{subject}.dat') 
    dat = pd.read_csv(dat_file_name, sep='\t')

    #  ============== PROCESSING THE RAW BEHAVIOURAL DATA ==============
    oldblock = -1
    # loop on trials:
    for i in range(dat.shape[0]):
        if dat['BN'][i] != oldblock:
            print(f'Processing block {dat["BN"][i]}')
            # load the .mov file:
            ext = int(dat['BN'][i])
            mov = ds.movload(os.path.join(DATA_PATH, 'behavioural', f'subj{subject}', f'efcp_{subject}_{ext:02d}.mov'))
            oldblock = dat['BN'][i]
        # trial routine:
        C = trial_routine(dat.iloc[[i]])

        # append the trial to the dataframe:
        D = pd.concat([D, C], ignore_index=True)

        # smooth the forces with Gaussian filter:
        tmp_mov = mov[dat['TN'][i]-1]
        tmp_mov[:,3:] = gaussian_filter1d(tmp_mov[:,3:], sigma=1.0, axis=0)
        # add the mov trial in the move dataframe:
        tmp = pd.DataFrame({'sn': np.full_like(tmp_mov[:,0], subject), 'BN': np.full_like(tmp_mov[:,0],dat['BN'][i]), 'TN': np.full_like(tmp_mov[:,0],dat['TN'][i]), 'trial_correct': np.full_like(tmp_mov[:,0],dat['trialCorr'][i]), 
                            'state': tmp_mov[:,0], 'time': tmp_mov[:,2],
                            'f1': tmp_mov[:,13], 'f2': tmp_mov[:,14], 'f3': tmp_mov[:,15], 'f4': tmp_mov[:,16], 'f5': tmp_mov[:,17]})
        df_mov = pd.concat([df_mov, tmp], ignore_index=True)

    #  ============== PROCESSING THE RAW EMG DATA ==============
    # iterate through emg files and load:
    col_names = ['sn','BN','TN','trial_correct','state','time','e1', 'e2', 'e3', 'e4', 'e5', 'f1', 'f2', 'f3', 'f4', 'f5']
    df_emg = pd.DataFrame(columns=col_names)
    uniqueBN = np.unique(D.BN)
    for i, bn in enumerate(uniqueBN):
        # getting the name of the file:
        fname = os.path.join(DATA_PATH, 'emg', f'subj{subject}', f'emg_run{bn}.csv')

        # load and clean up the file:
        file = pd.read_csv(fname, header=None, delimiter=',', skiprows=5) 
        file.drop(index=1, inplace=True)
        file.reset_index(drop=True, inplace=True)
        # header columns names:
        header = file.iloc[0].values
        # the emg signals:
        data = file[1:]
        data = data.apply(pd.to_numeric, errors='coerce').astype(float)

        # wanted channel names:
        channel_names = ['Analog 1', 'ext_D1', 'ext_D2', 'ext_D3', 'ext_D4', 'ext_D5', 'flx_D1', 'flx_D2', 'flx_D3', 'flx_D4', 'flx_D5']

        # Extract the wanted signals from the csv file
        raw_emg = np.zeros((data.shape[0],2*len(channel_names)), dtype=np.float32)
        col = 0
        for name in channel_names:
            for i_col, col_name in enumerate(header):
                if name in col_name:
                    raw_emg[:,col] = data.iloc[0:, i_col].to_numpy()
                    col = col+1

        # find trial indices based on triggers:
        t_trig = raw_emg[:,0]
        trig = raw_emg[:,1]
        t_signal = raw_emg[:,2]
        threshold = 0.05
        trial_start_idx, trial_end_idx, riseTime, fallTime  = emg_handler.find_trigger_rise_edge(t_trig, trig, t_signal, 
                                                                                                threshold=threshold, debug=debug)
        # get the full emg signal:
        sig = raw_emg[:,3::2]
        sig = np.nan_to_num(sig, nan=0.0)

        # bandpass filter the whole signal:
        sos = signal.butter(2, bpf, btype='bandpass', fs=fs_emg, output='sos')  # using 'sos' to avoid numerical errors
        sig_filtered = signal.sosfiltfilt(sos, sig, axis=0)
        
        sos_lpf = signal.butter(2, lpf, btype='lowpass', fs=fs_emg, output='sos')
        # get the trials:
        for i_trial, idx in enumerate(trial_start_idx):
            print(i_trial)
            trial_sig = sig_filtered[idx:trial_end_idx[i_trial]]
            trial_time = (t_signal[idx:trial_end_idx[i_trial]] - t_signal[idx]) * 1000  # in ms

            # de-mean rectify the signal:
            trial_sig = np.abs(trial_sig - np.mean(trial_sig, axis=0))
            
            # low-pass filter the demeaned signal:
            trial_sig_filtered = signal.sosfiltfilt(sos_lpf, trial_sig, axis=0)

            
            tmp_sn = np.full_like(trial_sig_filtered[:,0], subject, dtype=int)
            tmp_BN = np.full_like(trial_sig_filtered[:,0], bn, dtype=int)
            tmp_TN = np.full_like(trial_sig_filtered[:,0], i_trial+1, dtype=int)
            tmp_trial_correct = np.full_like(trial_sig_filtered[:,0], D[(D['BN'] == bn) & (D['TN'] == i_trial+1)]['trial_correct'].values, dtype=int)
            
            # match mov with emg:
            mov_time = df_mov[(df_mov['BN'] == bn) & (df_mov['TN'] == i_trial+1)]['time'].values
            mov_states = df_mov[(df_mov['BN'] == bn) & (df_mov['TN'] == i_trial+1)]['state'].values
            unique_states = np.unique(mov_states)
            tmp_state = np.full_like(trial_sig_filtered[:,0], 0)
            for j, state in enumerate(unique_states):
                first_index = np.where(mov_states == state)[0][0]
                last_index = np.where(mov_states == state)[0][-1]

                t1 = mov_time[first_index]
                if j != len(unique_states)-1:
                    t2 = mov_time[last_index+1]
                else:
                    t2 = mov_time[last_index]

                idx_time = np.where((trial_time >= t1) & (trial_time < t2))[0]
                tmp_state[idx_time] = state
            
            # fill in the dataframe:
            tmp = pd.DataFrame({'sn': tmp_sn, 'BN': tmp_BN, 'TN': tmp_TN, 'trial_correct': tmp_trial_correct, 
                            'state': tmp_state, 'time': trial_time,
                            'e1': trial_sig_filtered[:,0], 'e2': trial_sig_filtered[:,1], 'e3': trial_sig_filtered[:,2], 'e4': trial_sig_filtered[:,3], 'e5': trial_sig_filtered[:,4], 
                            'f1': trial_sig_filtered[:,5], 'f2': trial_sig_filtered[:,6], 'f3': trial_sig_filtered[:,7], 'f4': trial_sig_filtered[:,8], 'f5': trial_sig_filtered[:,9]})
            df_emg = pd.concat([df_emg, tmp], ignore_index=True)

    # save the data frames:
    D.to_csv(os.path.join(ANALYSIS_PATH, f'efcp_{subject}.csv'), index=False)
    df_mov.to_csv(os.path.join(ANALYSIS_PATH, f'efcp_{subject}_mov.csv'), index=False)
    df_emg.to_csv(os.path.join(ANALYSIS_PATH, f'efcp_{subject}_emg.csv'), index=False)

    return D, df_mov, df_emg

def ppp(subjNum, smoothing_window=30, fs=500):
    datFileName = os.path.join(DATA_PATH, 'behavioural', f'subj{subjNum}', f'efcp_{subjNum}.dat')   # input .dat file
    outFileName = os.path.join(ANALYSIS_PATH, f'subj{subjNum}', f'efcp_{subjNum}.csv')   # input .dat file

    D = pd.read_table(datFileName)
    D = D.loc[:, ~D.columns.str.contains('^Unnamed')]

    # cleaning up the D dataframe:
    D.rename(columns={'subNum': 'sn'}, inplace=True)
    D.loc[D['RT'] != 10000, 'RT'] = D['RT'] - 1400
    D.rename(columns={'RT': 'ET'}, inplace=True)


    # loading mov files and appending each block to movList:
    oldBlock = -1
    movList = []
    for i in range(len(D.BN)):
        if (oldBlock != D.BN[i]):
            # load mov file
            movPath = scriptPath + '/data/' + subjName + '/efc1_' + subjName[-2:] + '_' + '{:02d}.mov'.format(D.BN[i])
            mov = dataLoader.movload(movPath)
            movList.extend(mov)
            # print(mov[0])
            oldBlock = D.BN[i]
            print(len(movList))
    
    # adding the mov data to the dataframe
    D['mov'] = movList

    # loading emg data:
    emgList = [] # list to contain all emg trials

    # iterate through emg files and load:
    uniqueBN = np.unique(D.BN)
    for i in range(len(uniqueBN)):
        # getting the name of the file:
        fname = scriptPath + '/data/' + subjName + '/efc1_EMG_' + subjName[-2:] + '_' + '{:02d}.csv'.format(uniqueBN[i])

        # loading emg and separating trials:
        emg_selected, fs = dataLoader.emgload(fname, riseThresh=0.5, fallThresh=0.5, debug=0)

        # down sampling the signals:
        emg_selected, fs = emgHandler.downsample_emg(emg_selected, fs, target_fs=1000, debug=0)

        # filtering the signals - bandpass:
        emg_selected = emgHandler.filter_emg(emg_selected, fs=fs, low=20, high=500, order=2, debug=0)

        # rectifying the signals:
        emg_selected = emgHandler.rectify_emg(emg_selected, debug=0)

        # adding emg data of trials to emgList:
        emgList.extend(emg_selected)

    # adding emg data to the dataframe:
    D['emg'] = emgList

    return D