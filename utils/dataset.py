import numpy as np
import pandas as pd
import emg_handler

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
            print('Trial: ', line.split()[1])
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
    # loading emg file:
    

    # loading emg file:
    column_names = [i for i in range(0, 7)] # hard coded number of columns in dataframe - ideally you should calculate maximum of number of columns for future use of the code
    data = pd.read_csv(fname, header=None, delimiter=',', names=column_names)   # loads .csv file into a dataframe
    fs = float(data[0][6].split()[0])   # getting the sampling rate of the data - hard coded to get fs. Also, assumed that all emg channels have the same fs.
    raw_emg = pd.DataFrame(data.iloc[7:], dtype=float)  # making a new dataframe containing only the raw emg data - hard coded row num
    raw_emg = raw_emg.reset_index(drop=True)
    
    # finding triggers:
    trig = raw_emg[0]   # trigger channel
    riseIdx, fallIdx = emg_handler.find_trigger_rise_edge(trig, fs, riseThresh, fallThresh, debug) # trigger detector function
    
    # slecting emg data of each trial based on triggers
    emg_selected = [] # empty list to contain emg data
    for i in range(len(riseIdx)):
        # starting index of trial i:
        idxStart = riseIdx[i]

        # ending index of trial i:
        idxEnd = fallIdx[i]     

        # selecting emg data. Channel 0 is trigger so it's not included:
        emg_tmp = raw_emg.iloc[idxStart:idxEnd,1:].to_numpy()  

        # append trial's emg data to emg_selected , the format of emg_tmp is (N by ChannelNum):
        emg_selected.append(emg_tmp)   
        # print("Trial {}:\n".format(i), len(emg_tmp[:,0])) # a sanity check
    
    # print(np.shape(emg_selected[99])) # sanity check
    return emg_selected, fs

def within_subj_var(data, partition_vec, cond_vec, subj_vec, subtract_mean=True):
    '''
        Estimate the within subject and noise variance for each subject.
    
        Args:
            data: 2D numpy array of shape (N-regressors by P-voxels)
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

            # mean of voxels across conditions for each partition:
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