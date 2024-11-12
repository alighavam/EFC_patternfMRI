import numpy as np

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