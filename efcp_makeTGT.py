import argparse
import random
import os
import numpy as np
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--subNum', default='101',
                        help='Participant ID (e.g., subj100, subj101, ...)')
    
    args = parser.parse_args()
    subNum = args.subNum

    # chord list:
    set1 = np.array([91999,99199,99919,99991,92999,99299,99929,99992,
                     91111,91112,91121,91211,92111,91122,91212,91221,92112,92121,92211,91222,92122,92212,92221,92222])
    
    set2 = np.array([91999,99199,99919,92999,99299,99929,
                     91199,91919,99119,91299,91929,92199,92919,99129,99219,92299,92929,99229,
                     91119,91129,91219,92119,91229,92129,92219,92229])

    # make set 1 trials:
    nRep = 5 # number of repetition of each chord
    nVisits = 2 # number of visits of each chord
    nRuns = 8

    set1_trials = np.repeat(set1, nVisits)
    np.random.shuffle(set1_trials)
    set1_trials = np.repeat(set1_trials, nRep)

    # split the trials into runs:
    set1_runs = np.split(set1_trials, nRuns)

    # make set 2 trials:
    nRep = 5 # number of repetition of each chord
    nVisits = 2 # number of visits of each chord

    set2_trials = np.repeat(set2, nVisits)
    np.random.shuffle(set2_trials)
    set2_trials = np.repeat(set2_trials, nRep)

    set2_runs = np.split(set2_trials[:240], 8)
    set2_runs.append(set2_trials[240:])

    # all runs:
    runs = set1_runs + set2_runs
    
    # params:
    planTime = 0  # time for planning, if 0, go-cue is presented immediately
    execMaxTime = 10000  # maximum allowed time for execution
    feedbackTime = 500  # time to present feedback

    # make tgt files:
    # target file columns:
    column_names = ['subNum', 'chordID', 'planTime', 'execMaxTime', 'feedbackTime', 'iti']
    for r, run in enumerate(runs):
        # building the dataframe:
        df = pd.DataFrame(columns=column_names)
        df['chordID'] = run
        df['subNum'] = np.full_like(run, subNum)
        df['planTime'] = np.full_like(run, planTime)
        df['execMaxTime'] = np.full_like(run, execMaxTime)
        df['feedbackTime'] = np.full_like(run, feedbackTime)
        df['iti'] = np.random.randint(500, 1500, len(run))

        # saving the tgt file:
        fname = f'efcp_s{subNum}_run{r+1}.tgt'
        df.to_csv(os.path.join('target', fname), sep='\t', index=False)
        print(f'{fname} saved!')
        
    

