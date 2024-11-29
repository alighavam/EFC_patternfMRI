import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal

def find_trigger_rise_edge(t_trig, trig, t_signal, threshold=0.1, debug=0):    # detect rising and falling edges of the trigger signal
    '''
        Description: Detects triggers from the trigger channel of the data.

        Inputs:
        t_trig: time vector of the trigger signal.

        trig: trigger signal

        t_signal: time vector of the emg data. Other channels might have different sampling rate than the trigger channel.
                  For some reasons, this is the case in the DELSYS system! data is 2148.1481 Hz and trigger is 2222.2222 Hz.

        riseThresh: The threshold to detect the rising edge of the trigger. Note that trigger data is normalized to its maximum absolute value in this function.
        So threhshold values are from 0 to 1. The rising edge in my data is the start of a trial.

        fallThresh: Same as riseThresh but for falling edges. The fall edge denotes the end of a trial in my data.
        
        debug: debug mode. Plots the detected triggers for eye insepction. Also, prints some sanity checks in the console. Make sure to run the function with the debug mode 1 
        the first time that you are detecting triggers in a trigger channel. Basically run in debug mode every block of the emg data.

        Returns:
        trial_start_idx: Contains the indices corresponding to the start of the trials. Indices correspond to t_signal not t_trig.
        trial_end_idx: Contains the indices corresponding to the end of the trials. Indices correspond to t_signal.
        riseTime: Contains the time of the rising edge of the trigger. The time is with respect to the t_signal.
        fallTime: Contains the time of the falling edge of the trigger. The time is with respect to the t_signal.
    ''' 

    # normalizing the trigger signal to its abs maximum
    trig = trig/np.amax(np.absolute(trig))

    # detect rising edges:
    derivative = np.diff(trig)
    crossings = np.where((derivative[:-1] < threshold) & (derivative[1:] >= threshold))[0] + 1
    # Remove indices that are too close to each other
    min_distance = 200  # Minimum number of samples between detections
    filtered_crossings = [crossings[0]]
    for idx in crossings[1:]:
        if idx - filtered_crossings[-1] >= min_distance:
            filtered_crossings.append(idx)
    riseIdx = np.array(filtered_crossings)

    # detecting falling edges:
    derivative = -np.diff(trig)
    crossings = np.where((derivative[:-1] < threshold) & (derivative[1:] >= threshold))[0] + 1
    # Remove indices that are too close to each other
    min_distance = 200  # Minimum number of samples between detections
    filtered_crossings = [crossings[0]]
    for idx in crossings[1:]:
        if idx - filtered_crossings[-1] >= min_distance:
            filtered_crossings.append(idx)
    fallIdx = np.array(filtered_crossings)

    # debug mode to make sure of the trigger detection:
    if debug:
        # number of detected triggers:
        print("\n\n======== Trigger Detection Results: ======== \n")
        print("Num Rise Trigger = {:d}".format(len(riseIdx)))
        print("Num Fall Triggers = {:d}".format(len(fallIdx)))
        print("Two numbers should be equal to the number of trials.\n")

        # falling edge triggers should be always after rising edge:
        if len(fallIdx) == len(riseIdx):
            diffRiseFall = fallIdx - riseIdx
            numNegative = sum(diffRiseFall <= 0)
            print("\nNumber of non-positive fall-rise edges = {:d}".format(numNegative))
            print("This value should be 0.\n")

        # plotting trigger signal along with the detected triggers:
        plt.figure(figsize=(14,3))
        plt.plot(t_trig, trig, label='trig', zorder=1)
        plt.scatter(t_trig[riseIdx], [trig[index] for index in riseIdx], s=6, c=[1,0,0], label='Rising Edges')
        plt.scatter(t_trig[fallIdx], [trig[index] for index in fallIdx], s=6, c=[0,0,0], label='Falling Edges')
        # plt.xlim(4.3e+05,5e+05)
        plt.show()
    
    # find the trigger indices with respect to t_signal:
    trig_riseTime = t_trig[riseIdx]
    trig_fallTime = t_trig[fallIdx]

    trial_start_idx = np.zeros(len(trig_riseTime), dtype=int)
    trial_end_idx = np.zeros(len(trig_fallTime), dtype=int)
    riseTime = np.zeros(len(trig_riseTime))
    fallTime = np.zeros(len(trig_fallTime))
    for i in range(len(trig_riseTime)):
        # get rise index:
        tDiff = np.absolute(t_signal - trig_riseTime[i])
        trial_start_idx[i] = np.argmin(tDiff)
        riseTime[i] = t_signal[np.argmin(tDiff)]

        tDiff = np.absolute(t_signal - trig_fallTime[i])
        trial_end_idx[i] = np.argmin(tDiff)
        fallTime[i] = t_signal[np.argmin(tDiff)]

    return trial_start_idx, trial_end_idx, riseTime, fallTime 

def downsample_emg(emg, fs, target_fs=1000, debug=0):
    
    # resampled emg:
    emg_resampled = []

    # designing lowpass filter with cutoff frequency of target_fs/2:
    sos = signal.butter(2, int(target_fs/2), btype='lowpass', fs=fs, output='sos')

    # iterating through signals and downsampling with zero-phase anti-aliasing filter:
    for i in range(len(emg)):
        # selecting the emg signal of trial i:
        emg_trial = emg[i]

        # number of samples in the resampled signal:
        target_len = int(np.floor(len(emg_trial[:,0])*target_fs/fs))

        # making an empty array to contain the resampled signal:
        emg_trial_resampled = np.empty((target_len, np.shape(emg_trial)[1]))

        # iterating through emg channels:
        for ch in range(np.shape(emg_trial)[1]):
            # selecting the signal:
            sig = emg_trial[:,ch]

            # zero-phase low pass filtering the signal to avoid aliasing:
            sig = signal.sosfiltfilt(sos, sig)

            # downsampling the signal:
            sig_resampled = signal.resample(sig, target_len)

            # appending the resampled signal to the resampled trial:
            emg_trial_resampled[:,ch] = sig_resampled

            # plotting resampled against original signal:
            if debug and i==0 and ch==0:
                # time vector for plotting:
                t_orig = np.linspace(0,len(sig)/fs,len(sig))
                t_resampled = np.linspace(0,len(sig_resampled)/target_fs,len(sig_resampled))
                print("Original Signal Length = {:d}".format(len(sig)))
                print("Resampled Signal Length = {:d}".format(len(sig_resampled)))

                # plotting:
                plt.figure()
                plt.plot(t_orig, sig, label='Original Signal')
                plt.plot(t_resampled, sig_resampled, label='Resampled Signal')
                plt.legend()
                plt.show()
        
        # appending the resampled trial to the resampled emg:
        emg_resampled.append(emg_trial_resampled)

    return emg_resampled, target_fs

def filter_emg(emg, fs, low=20, high=500, order=2, debug=0):
    # filtered emg:
    emg_filtered = []

    # designing bandpass filter:
    sos = signal.butter(order, [low, high-1], btype='bandpass', fs=fs, output='sos')  # using 'sos' to avoid numerical errors

    # iterating through signals and filtering:
    for i in range(len(emg)):
        # selecting the emg signal of trial i:
        emg_trial = emg[i]

        # making an empty array to contain the filtered signal:
        emg_trial_filtered = np.empty((np.shape(emg_trial)[0], np.shape(emg_trial)[1]))

        # iterating through emg channels:
        for ch in range(np.shape(emg_trial)[1]):
            # selecting the signal:
            sig = emg_trial[:,ch]

            # filtering the signal:
            sig_filtered = signal.sosfiltfilt(sos, sig)

            # appending the filtered signal to the emg_trial_filtered:
            emg_trial_filtered[:,ch] = sig_filtered

            # plotting filtered against original signal:
            if debug and i==0 and ch==0:
                # time vector for plotting:
                t = np.linspace(0,len(sig)/fs,len(sig))

                # plotting:
                plt.figure()
                plt.plot(t, sig, label='Original Signal')
                plt.plot(t, sig_filtered, label='Filtered Signal')
                plt.legend()
                plt.show()
        
        # appending the resampled trial to the resampled emg:
        emg_filtered.append(emg_trial_filtered)

    return emg_filtered
    
def rectify_emg(emg, debug=0):
    # rectified emg:
    emg_rectified = []

    # iterating through signals and rectifying:
    for i in range(len(emg)):
        # selecting the emg signal of trial i:
        emg_trial = emg[i]

        # making an empty array to contain the rectified signal:
        emg_trial_rectified = np.empty((np.shape(emg_trial)[0], np.shape(emg_trial)[1]))

        # iterating through emg channels:
        for ch in range(np.shape(emg_trial)[1]):
            # selecting the signal:
            sig = emg_trial[:,ch]

            # rectifying the signal:
            sig_rectified = np.absolute(sig)

            # appending the rectified signal to the emg_trial_rectified:
            emg_trial_rectified[:,ch] = sig_rectified

            # plotting rectified against original signal:
            if debug and i==0 and ch==0:
                # time vector for plotting:
                t = np.linspace(0,len(sig),len(sig))

                # plotting:
                plt.figure()
                plt.plot(t, sig, label='Original Signal')
                plt.plot(t, sig_rectified, label='Rectified Signal')
                plt.legend()
                plt.show()
        
        # appending the rectified trial:
        emg_rectified.append(emg_trial_rectified)

    return emg_rectified
