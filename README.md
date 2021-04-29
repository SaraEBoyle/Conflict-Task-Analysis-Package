# Conflict-Task-Analysis-Package

Photometry analysis for platform mediated conflict task

Inputs:
1. Analog input file (.csv)- file containing timing data for syncing each bpod trial to photometry data. Ensures timing is stable throughout entire session.
2. Photometry data file (.csv)- raw data array (n x m) where n is number of data entries, m = 1 is time, and m = 2, 3, ... is fluorescence data from each fiber
3. Bpod file (.mat)- contains all event and state data from bpod events
4. Centroid file (.csv)- Tracking data obtained in real time using Bonsai

Pipeline:
