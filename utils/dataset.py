import numpy as np
import pandas as pd
from utils import emg_handler

def movload(fname):
    '''
        Loads .mov files given the path of the file. 
        The .mov files have and arbitrary format hence the need for a custom function to parse the data correctly.

        Input:
        fname: path to the .mov file

        Returns:
        A: a list of numpy arrays. Each element of the list is a numpy array containing the data of a trial.
        To know what each column of the numpy array represents, you need to take a look at the experiment code on robotcode repo on
        diedrichsen lab github.
    '''
    A = []
    fid = open(fname, 'rt')
    if fid == -1:
        raise Exception('Could not open ' + fname)

    trial = 0
    for line in fid:
        if line[0] == 'T':
            # print('Trial: ', line.split()[1])
            a = int(line.split()[1])
            trial += 1
            if a != trial:
                print('Trials out of sequence')
                trial = a
            A.append([])
            A[trial-1] = np.empty((0,23))
        else:
            lineData = line.strip().split('\t')
            a = np.array([float(x) for x in lineData], ndmin=2)
            # print(a)
            A[trial-1] = np.vstack((A[trial-1],a))
            # A[trial-1].extend(a)

    fid.close()
    return A

def emgload(fname, channel_names, riseThresh=0.5, fallThresh=0.5, debug=0):
    '''
        Description: Loads and handles the DELSYS emg data from fname.csv file. 

        Inputs:
        fname: the .csv file name

        riseThresh, fallThresh, debug: refer to emgHandler.find_trigger_rise_edge function.

        Returns:
        emg_selected: python list with len <number of trials>. Each list element is emg data for each tiral. It is a numpy array with the format of (N by ChannelNum). 
        N is the len of that trial. ChannelNum is the number of emg channels.

        fs: sampling frequency of the emg data.
    '''
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
    
    return raw_emg

def within_subj_var(data, partition_vec, cond_vec, subj_vec, subtract_mean=True):
    '''
        Estimate the within subject and noise variance for each subject.
    
        Args:
            data: 2D numpy array of shape (N-regressors by P-voxels) nan must be removed.
            partition_vec: 1D numpy array of shape (N-regressors) with partition index
            cond_vec: 1D numpy array of shape (N-regressors) with condition index
            subj_vec: 1D numpy array of shape (P-voxels) with subject index
            subtract_mean: Subtract the mean of voxels across conditions within a run.

        Returns:
            v_s:  1D array containing the subject variance.
            v_se: 1D array containing subject + noise variance.
    '''

    # In case partition_vec was not contiguous, e.g.,: [1, 1, 2, 2, 1, 1, 2, 2] instead of [1,1,1,1,2,2,2,2].
    # First, make partition indices contiguous by sorting the rows:
    sorted_indices = np.argsort(partition_vec)
    data = data[sorted_indices]
    partition_vec = partition_vec[sorted_indices]
    cond_vec = cond_vec[sorted_indices]

    cond = np.unique(cond_vec)
    subj = np.unique(subj_vec)
    partition = np.unique(partition_vec)

    v_s = np.zeros(len(subj))
    v_se = np.zeros(len(subj))
    # loop on subj:
    for sn in subj:
        Y = data[:, subj_vec == sn]

        # subtract mean of voxels across conditions within each run:
        if subtract_mean:
            N, P = Y.shape

            # Reshape Y to separate each partition
            Y_reshaped = Y.reshape(partition.shape[0], cond.shape[0], P)

            # mean of voxels over conditions for each partition:
            partition_means = Y_reshaped.mean(axis=1, keepdims=True)

            # subtract the partition means from the original reshaped Y and rehsape back to original:
            Y = (Y_reshaped - partition_means).reshape(N, P)

            cov_Y = Y @ Y.T / Y.shape[1]

            # avg of the main diagonal:
            avg_main_diag = np.sum(np.diag(cov_Y))/(len(cond)*len(partition))
            
            # avg of the main off-diagonal:
            mask = np.kron(np.eye(len(partition)), np.ones((len(cond),len(cond))))
            mask = mask - np.eye(mask.shape[0])
            avg_main_off_diag = np.sum(cov_Y * mask)/(np.sum(mask))
            
            # within partition variance:
            v_se[sn-1] = avg_main_diag

            # avg across session diagonals:
            mask = np.kron(np.ones((len(partition), len(partition))), np.eye(len(cond)))
            mask = mask - np.eye(mask.shape[0])
            avg_across_diag = np.sum(cov_Y * mask)/(np.sum(mask))

            # avg across session off-diagonals:
            mask = np.kron(1-np.eye(len(partition)), np.ones((len(cond), len(cond))))
            mask = mask - np.kron(np.ones((len(partition), len(partition))), np.eye(len(cond))) + np.eye(mask.shape[0])
            avg_across_off_diag = np.sum(cov_Y * mask)/(np.sum(mask))

            # across partition variance:            
            v_s[sn-1] = avg_across_diag

        else:
            pass

    return v_s, v_se