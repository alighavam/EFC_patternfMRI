import argparse
import random
import os
import numpy as np
import pandas as pd

def makeTGT_emg_deprecated(subNum):
    '''
    make target files for EMG experiment
    '''
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

def makeTGT_training_deprecated(subNum):
    ''' 
    make target files for training experiment
    '''
    # chord list:
    set = np.array([91999,99199,99919,92999,99299,99929,
                     91199,91919,99119,91299,91929,92199,92919,99129,99219,92299,92929,99229,
                     91119,91129,91219,92119,91229,92129,92219,92229])

    nDays = 4
    nRep = 5 # number of repetition of each chord
    nVisits = 3 # number of visits of each chord

    # params:
    planTime = 0  # time for planning, if 0, go-cue is presented immediately
    execMaxTime = 10000  # maximum allowed time for execution
    holdTime = 1000  # hold time of the chord
    feedbackTime = 500  # time to present feedback
    iti_range = [500, 1000]
    
    for i in range(nDays):
        trials = np.repeat(set, nVisits)
        np.random.shuffle(trials)
        trials = np.repeat(trials, nRep)

        runs = np.split(trials[:360], 9)
        runs.append(trials[360:])
        
        # make tgt files:
        # target file columns:
        column_names = ['subNum', 'chordID', 'planTime', 'execMaxTime', 'holdTime', 'feedbackTime', 'iti']
        for r, run in enumerate(runs):
            # building the dataframe:
            df = pd.DataFrame(columns=column_names)
            df['chordID'] = run
            df['subNum'] = np.full_like(run, subNum)
            df['planTime'] = np.full_like(run, planTime)
            df['execMaxTime'] = np.full_like(run, execMaxTime)
            df['holdTime'] = np.full_like(run, holdTime)
            df['feedbackTime'] = np.full_like(run, feedbackTime)
            df['iti'] = np.random.randint(iti_range[0], iti_range[1], len(run))

            # saving the tgt file:
            fname = f'efcp_s{subNum}_train_day{i+1}_run{r+1:02d}.tgt'
            df.to_csv(os.path.join('target', fname), sep='\t', index=False)
            print(f'{fname} saved!')

def makeTGT_training_day1(subNum):
    ''' 
    EFCP DAY1: 13 blocks of 50 trials = 650 trials. exec_max_time = 5s, repetition = 5, hold time = 600ms
    feedback_time+iti jittered.
    hold_till_end = 0
    Estimated time:
    '''

    # chord list:
    set = np.array([91999,99199,99919,92999,99299,99929,
                     91199,91919,99119,91299,91929,92199,92919,99129,99219,92299,92929,99229,
                     91119,91129,91219,92119,91229,92129,92219,92229])

    nRep = 5 # number of repetition of each chord
    nVisits = 5 # number of visits of each chord

    # params:
    hold_till_end = 0
    planTime = 0  # time for planning, if 0, go-cue is presented immediately
    execMaxTime = 5000  # maximum allowed time for execution
    holdTime = 600  # hold time of the chord
    feedbackTime = 500  # time to present feedback
    iti_range = [500, 1000]
    
    trials = np.repeat(set, nVisits)
    np.random.shuffle(trials)
    trials = np.repeat(trials, nRep)

    runs = np.split(trials, 13)
    
    # make tgt files:
    # target file columns:
    column_names = ['subNum', 'chordID', 'planTime', 'execMaxTime', 'holdTime', 'hold_till_end', 'feedbackTime', 'iti']
    for r, run in enumerate(runs):
        # building the dataframe:
        df = pd.DataFrame(columns=column_names)
        df['chordID'] = run
        df['subNum'] = np.full_like(run, subNum)
        df['planTime'] = np.full_like(run, planTime)
        df['execMaxTime'] = np.full_like(run, execMaxTime)
        df['holdTime'] = np.full_like(run, holdTime)
        df['hold_till_end'] = np.full_like(run, hold_till_end)
        df['feedbackTime'] = np.full_like(run, feedbackTime)
        df['iti'] = np.random.randint(iti_range[0], iti_range[1], len(run))
        
        # saving the tgt file:
        fname = f'efcp_s{subNum}_day1_run{r+1:02d}.tgt'
        df.to_csv(os.path.join('target', fname), sep='\t', index=False)
        print(f'{fname} saved!')

def makeTGT_training_day2(subNum):
    ''' 
    EFCP DAY2: 12 blocks of 50 trials + 1 block of 24 = 624 = 26*24 trials. exec_max_time = 3s, hold time = 600ms, repetition = 2
    feedback_time+iti jittered. Chords repeat for 2 consecutive trials.
    hold_till_end = 0
    Estimated time: 
    '''

    # chord list:
    set = np.array([91999,99199,99919,92999,99299,99929,
                     91199,91919,99119,91299,91929,92199,92919,99129,99219,92299,92929,99229,
                     91119,91129,91219,92119,91229,92129,92219,92229])

    nRep = 2 # number of repetition of each chord
    nVisits = 12 # number of visits of each chord

    # params:
    hold_till_end = 0
    planTime = 0  # time for planning, if 0, go-cue is presented immediately
    execMaxTime = 3000  # maximum allowed time for execution
    holdTime = 600  # hold time of the chord
    feedbackTime = 500  # time to present feedback
    iti_range = [500, 1000]
    
    trials = np.repeat(set, nVisits)
    np.random.shuffle(trials)
    trials = np.repeat(trials, nRep)

    runs = np.split(trials[:600], 12)
    runs.append(trials[600:])
    
    # make tgt files:
    # target file columns:
    column_names = ['subNum', 'chordID', 'planTime', 'execMaxTime', 'holdTime', 'hold_till_end' ,'feedbackTime', 'iti']
    for r, run in enumerate(runs):
        # building the dataframe:
        df = pd.DataFrame(columns=column_names)
        df['chordID'] = run
        df['subNum'] = np.full_like(run, subNum)
        df['planTime'] = np.full_like(run, planTime)
        df['execMaxTime'] = np.full_like(run, execMaxTime)
        df['holdTime'] = np.full_like(run, holdTime)
        df['hold_till_end'] = np.full_like(run, hold_till_end)
        df['feedbackTime'] = np.full_like(run, feedbackTime)
        df['iti'] = np.random.randint(iti_range[0], iti_range[1], len(run))

        # saving the tgt file:
        fname = f'efcp_s{subNum}_day2_run{r+1:02d}.tgt'
        df.to_csv(os.path.join('target', fname), sep='\t', index=False)
        print(f'{fname} saved!')

def makeTGT_emg_day3(subNum):
    ''' 
    EFCP DAY 3: 10 blocks of 26 trials = 260 trials. exec_max_time = 3s fix, repetition = 2, hold time = 600ms,
    feedback_time+iti jittered.
    hold_till_end = 1
    Estimated time: 
    '''

    # chord list:
    set = np.array([91999,99199,99919,92999,99299,99929,
                     91199,91919,99119,91299,91929,92199,92919,99129,99219,92299,92929,99229,
                     91119,91129,91219,92119,91229,92129,92219,92229])

    nRep = 2 # number of repetition of each chord
    nVisits = 5 # number of visits of each chord

    # params:
    hold_till_end = 1
    planTime = 0  # time for planning, if 0, go-cue is presented immediately
    execMaxTime = 3000  # maximum allowed time for execution
    holdTime = 600  # hold time of the chord
    feedbackTime = 500  # time to present feedback
    iti_range = [500, 1000]
    
    trials = np.repeat(set, nVisits)
    np.random.shuffle(trials)
    trials = np.repeat(trials, nRep)

    runs = np.split(trials, 10)
    
    # make tgt files:
    # target file columns:
    column_names = ['subNum', 'chordID', 'planTime', 'execMaxTime', 'holdTime', 'hold_till_end', 'feedbackTime', 'iti']
    for r, run in enumerate(runs):
        # building the dataframe:
        df = pd.DataFrame(columns=column_names)
        df['chordID'] = run
        df['subNum'] = np.full_like(run, subNum)
        df['planTime'] = np.full_like(run, planTime)
        df['execMaxTime'] = np.full_like(run, execMaxTime)
        df['holdTime'] = np.full_like(run, holdTime)
        df['hold_till_end'] = np.full_like(run, hold_till_end)
        df['feedbackTime'] = np.full_like(run, feedbackTime)
        df['iti'] = np.random.randint(iti_range[0], iti_range[1], len(run))

        # saving the tgt file:
        fname = f'efcp_s{subNum}_day3_run{r+1:02d}.tgt'
        df.to_csv(os.path.join('target', fname), sep='\t', index=False)
        print(f'{fname} saved!')

def makeTGT_emg_day4_deprecated(subNum):
    ''' 
    EFCP DAY 4: 10 blocks of 50 trials + 1 block of 20 trials = 520 trials. exec_max_time = 3s fix, 
    feedback_time+iti=1s jittered. Chords repeat in batches of 2 trials. 
    Estimated Time: 25-mins electrode placement + 60-mins chord production + 30-mins natural = 120mins
    '''

    # chord list:
    set = np.array([91999,99199,99919,92999,99299,99929,
                     91199,91919,99119,91299,91929,92199,92919,99129,99219,92299,92929,99229,
                     91119,91129,91219,92119,91229,92129,92219,92229])

    nRep = 2 # number of repetition of each chord
    nVisits = 10 # number of visits of each chord

    # params:
    planTime = 0  # time for planning, if 0, go-cue is presented immediately
    execMaxTime = 3000  # maximum allowed time for execution
    holdTime = 1000  # hold time of the chord
    feedbackTime = 500  # time to present feedback
    iti_range = [500, 1000]
    
    trials = np.repeat(set, nVisits)
    np.random.shuffle(trials)
    trials = np.repeat(trials, nRep)

    runs = np.split(trials[:500], 10)
    runs.append(trials[500:])
    
    # make tgt files:
    # target file columns:
    column_names = ['subNum', 'chordID', 'planTime', 'execMaxTime', 'holdTime', 'feedbackTime', 'iti']
    for r, run in enumerate(runs):
        # building the dataframe:
        df = pd.DataFrame(columns=column_names)
        df['chordID'] = run
        df['subNum'] = np.full_like(run, subNum)
        df['planTime'] = np.full_like(run, planTime)
        df['execMaxTime'] = np.full_like(run, execMaxTime)
        df['holdTime'] = np.full_like(run, holdTime)
        df['feedbackTime'] = np.full_like(run, feedbackTime)
        df['iti'] = np.random.randint(iti_range[0], iti_range[1], len(run))

        # saving the tgt file:
        fname = f'efcp_s{subNum}_day4_run{r+1:02d}.tgt'
        df.to_csv(os.path.join('target', fname), sep='\t', index=False)
        print(f'{fname} saved!')

def makeTGT_fMRI(subNum):
    '''
        make target files for fMRI experiment
    '''
    # chord list:
    set = np.array([91999,99199,99919,92999,99299,99929,
                     91199,91919,99119,91299,91929,92199,92919,99129,99219,92299,92929,99229,
                     91119,91129,91219,92119,91229,92129,92219,92229])

    nDays = 2
    nRep = 2 # number of repetition of each chord
    nRuns = 8

    # params:
    # trial timings:
    planTime = 0  # time for planning, if 0, go-cue is presented immediately
    execMaxTime = 3000  # maximum allowed time for execution
    success_holdTime = 600  # hold time of the chord
    feedbackTime = 1000  # time to present feedback
    iti = 0

    # scanner related timings:
    start_acq = 3000    # time to start the acquisition. for signal to reach equilibrium
    buff = 1000         # buffer time between trials.
    rest = 15000        # rest interval durations.
    end_acq = 10500     # how long to wait after the last trial before stopping the acquisition. Because BOLD is slow.

    for i in range(nDays):
        for r in range(nRuns):
            np.random.shuffle(set)
            half1 = np.repeat(set, nRep)
            np.random.shuffle(set)
            half2 = np.repeat(set, nRep)
            trials = np.concatenate((half1, half2))
            
            # make tgt files:
            # target file columns:
            column_names = ['subNum', 'chordID', 'planTime', 'success_holdTime', 'execMaxTime', 'feedbackTime', 'iti', 'startTime', 'endTime']
            df = pd.DataFrame(columns=column_names)
            df['chordID'] = trials
            df['subNum'] = np.full_like(trials, subNum)
            df['planTime'] = np.full_like(trials, planTime)
            df['success_holdTime'] = np.full_like(trials, success_holdTime)
            df['execMaxTime'] = np.full_like(trials, execMaxTime)
            df['feedbackTime'] = np.full_like(trials, feedbackTime)
            df['iti'] = np.full_like(trials, iti, dtype=int)

            # make start times:
            df['startTime'] = start_acq + np.arange(len(trials))*(execMaxTime + feedbackTime + iti + buff)
            df.loc[26:, 'startTime'] += rest - (execMaxTime + feedbackTime + iti + buff)
            df.loc[52:, 'startTime'] += rest - (execMaxTime + feedbackTime + iti + buff)
            df.loc[78:, 'startTime'] += rest - (execMaxTime + feedbackTime + iti + buff)

            # make end times:
            df['endTime'] = np.zeros(len(trials), dtype=int)
            df.loc[df.index[-1], 'endTime'] += df.loc[df.index[-1], 'startTime'] + end_acq

            # saving the tgt file:
            fname = f'efcp_s{subNum}_scan_day{i+1}_run{r+1:02d}.tgt'
            df.to_csv(os.path.join('target', fname), sep='\t', index=False)
            print(f'{fname} saved!')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--subNum', default='999',
                        help='Participant ID (e.g., subj100, subj101, ...)')
    parser.add_argument('--experiment', default='fMRI_pilot1',
                        help='Type of experiment (e.g., train, emg, fMRI, etc.)')
    
    args = parser.parse_args()
    subNum = args.subNum
    experiment = args.experiment
    
    if experiment == 'fMRI_pilot1':
        makeTGT_training_day1(subNum)
        makeTGT_training_day2(subNum)
        makeTGT_emg_day3(subNum)
        makeTGT_fMRI(subNum)
    elif experiment == 'train':
        makeTGT_training_day1(subNum)
    elif experiment == 'fMRI':
        makeTGT_fMRI(subNum)
    else:
        raise ValueError('Unknown experiment type!')
        
    
