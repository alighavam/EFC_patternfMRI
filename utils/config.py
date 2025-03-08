import os
import platform

# Store the current working directory in a variable
WORKING_DIRECTORY = os.getcwd()

DATA_PATH = WORKING_DIRECTORY + '/data/efcp_EMG'
ANALYSIS_PATH = WORKING_DIRECTORY + '/analysis/'
FIGURE_PATH = WORKING_DIRECTORY + '/figures/'

# If on server:
if platform.system() == 'Linux':
    DATA_PATH = '/cifs/diedrichsen/data/...sth'
    ANALYSIS_PATH = '/cifs/diedrichsen/data/...sth'