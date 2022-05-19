Sorryfor the delay, I've been busy with interviewing for data engineering positions and as Karen can probably tell you, I'm not the best at replying right now...

DISCLAIMER:
These analysis scripts will ONLY WORK FOR ABF FILES. These scripts were not tested on traces recorded using sutterpatch software (the new rig) and were only used when analyzing traces from the old rig.

This repo holds two python scripts used for extracting spike timings from electrophysiology data. 
SpikeTimingAnalyis.py and SpikeTimingAnalysis_Contactenated_Files.py will prompt the user for 4 inputs in the following order:

1. a folder containing the file from the induction recording
2. the first frequency of the induction   
3. the second frequency of presynaptic stimulation during the induction
4. the delay between the ephys software delivering the presynaptic signal and the detection of the presynaptic artifact by the recording software. This delay is different for each experiment due to the presynaptic stimulating electrode being a different distance from the recording electrode in each sample.

inputs 2-4 Should (?) be optional if all you are interested in finding is the postsynaptic spike timings and nothing else (presynaptic inputs, relationship of presynaptic inputs to postsynaptic spikes, etc). However, this has not been tested.

There are 2 versions of the analysis script because Jin collected the data in different formats for some of the experiments. The induction was either recorded in 1 file or over 5 files. 

If the induction was recorded over 5 files. It is necessary to concatenate the files together using the analysis software for the old rig and then put the contacenated trace file in a folder to be analyzed, then use "SpikeTimingAnalysis_Concatenated_Files.py"

If the induction was recorded over 1 files. Then put the file containing the induction traces into a folder to be analyzed, then use "SpikeTimingAnalysis.py"

The output of these scripts will be a csv file. If you are interested in the postsynaptic action potentials, they are recorded in the "peak_indices" column. The action potential times are relative to the start of each stimulation in a given trace of the induction period. Absolute times were not necessary as these scripts were solely meant to compare STDP pre-post pairings. If you want absolute times you will have to edit the script. 

I have also included a list of the cells I analyzed in "Cells_to_analyze.txt". This list contains the cell name, the induction protocol used, age of the pup, and epsp offset (input #4 for the scripts).

The main libraries that were of use for me were as follows:

pyabf - python interface for working with axon binary format files (ABF is the proprietary file format used by Axon in the ephys recording software)
Docs: https://swharden.com/pyabf/

efel - electrophysiology feature extraction library
Docs: https://efel.readthedocs.io/en/latest/

Best of luck.
