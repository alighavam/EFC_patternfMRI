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

def subject_routine(subject, fs=500):
    """
    This function is used to preprocess the data of a subject
    
    params:
        subject: int, the subject number
        fs: int, the sampling frequency of the force data
    """
    # empty dataframe to store the data:
    df = pd.DataFrame()
    df_mov = pd.DataFrame(columns=['sn', 'BN', 'TN', 'trial_correct', 'state', 'time','f1','f2','f3','f4','f5'])

    # Load the .dat file:
    dat_file_name = os.path.join(DATA_PATH, 'behavioural', f'subj{subject}', f'efcp_{subject}.dat') 
    dat = pd.read_csv(dat_file_name, sep='\t')
    
    oldblock = -1
    # loop on trials:
    for i in range(dat.shape[0]):
        if dat['BN'][i] != oldblock:
            print(f'Processing block {dat['BN'][i]}')
            # load the .mov file:
            ext = int(dat['BN'][i])
            mov = ds.movload(os.path.join(DATA_PATH, 'behavioural', f'subj{subject}', f'efcp_{subject}.dat'))
            oldblock = dat['BN'][i]
        print(f'Processing trial {dat["TN"][i]}')
        # trial routine:
        C = trial_routine(dat.iloc[[i]])

        # append the trial to the dataframe:
        df = pd.concat([df, C], ignore_index=True)

        # smooth the forces with Gaussian filter:
        tmp_mov = mov[dat['TN'][i]-1]
        tmp_mov[:,3:] = gaussian_filter1d(tmp_mov[:,3:], sigma=1.0, axis=0)
        # add the mov trial in the move dataframe:
        tmp = pd.DataFrame({'sn': np.full_like(tmp_mov[:,0], subject), 'BN': np.full_like(tmp_mov[:,0],dat['BN'][i]), 'TN': np.full_like(tmp_mov[:,0],dat['TN'][i]), 'trial_correct': np.full_like(tmp_mov[:,0],dat['trialCorr'][i]), 
                            'state': tmp_mov[:,0], 'time': tmp_mov[:,2],
                            'f1': tmp_mov[:,13], 'f2': tmp_mov[:,14], 'f3': tmp_mov[:,15], 'f4': tmp_mov[:,16], 'f5': tmp_mov[:,17]})
        df_mov = pd.concat([df_mov, tmp], ignore_index=True)

    # sort the dataframes by day:
    df = df.sort_values(by='day', kind='mergesort')
    df_mov = df_mov.sort_values(by='day', kind='mergesort')

    # save the data frames:
    df.to_csv(os.path.join(ANALYSIS_PATH, f'efc2_{subject}.csv'), index=False)
    
    df_mov.reset_index(drop=True, inplace=True)
    df_mov.to_pickle(os.path.join(ANALYSIS_PATH, f'efc2_{subject}_mov.pkl'))
    
    return df, df_mov

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