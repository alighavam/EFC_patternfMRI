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
            vs, ve: 1D numpy array of shape (P-voxels) with within subject and noise variance
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

    return cov_Y
        

