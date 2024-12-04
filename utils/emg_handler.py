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
        trial_start_idx[i] = np.nanargmin(tDiff)
        riseTime[i] = t_signal[np.nanargmin(tDiff)]

        tDiff = np.absolute(t_signal - trig_fallTime[i])
        trial_end_idx[i] = np.nanargmin(tDiff)
        fallTime[i] = t_signal[np.nanargmin(tDiff)]

    return trial_start_idx, trial_end_idx, riseTime, fallTime 
