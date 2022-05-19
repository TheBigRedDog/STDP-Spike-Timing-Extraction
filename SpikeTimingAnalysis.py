# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 12:38:06 2021

@author: Cliff Moran
"""

import pyabf
import matplotlib.pyplot as plt
import numpy as np
import efel
import pandas as pd
import tkinter as tk
import os
from tkinter import filedialog
        
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return int(array[idx])

def plotting_accuracy_test(Sweep_Number, sweep_voltage, sweep_time, Stimulation_Number, current_injection_times, stim_data):
    xmin, xmax = current_injection_times[Stimulation_Number][0], current_injection_times[Stimulation_Number][1]
    # create a plot with time on the horizontal axis
    postsynaptic_indices = stim_data[0]['peak_indices'] 
    presynaptic_indices = stim_data[0]['pre_stim_indices']
    
    x = np.arange(450, 2549, 1)
    plt.figure(figsize=(8, 5))
    plt.plot(x, sweep_voltage[xmin-50:xmax+50], lw=.5, alpha=.75, color='#86bf91') # for plotting V(t)
    
    plt.margins(0, .5)
    plt.ylabel('Membrane Potential (mV)')
    plt.xlabel('Time (.1 ms)')
    plt.legend(loc="upper center", bbox_to_anchor=(0.5, 1.15))
 
    
    # now add the tags as vertical lines
    
    for index in presynaptic_indices:
       
        posX = index
        comment = 'Presynaptic'
        color = "steelblue"
        plt.axvline(posX, label=comment, color=color, linewidth=1, ls='--')
        
    for index in postsynaptic_indices:
      
        posX = index
        comment = 'Postsynaptic'
        color = "red" 
        plt.axvline(posX, label=comment, color=color, linewidth=1, ls='--')

        
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc="upper right")

    ax = plt.gca()
    
    ax.set_ylim([-40, 40]) # for plotting V(t)
    
    # display the plot
    plt.title("Sweep #%x, Stim #%x" % (Sweep_Number, Stimulation_Number))
    plt.show()        
    
    
    return 

def filter_APs(stim_results, stim_time, stim_voltage, sweep_num, stim_num):
    '''
    validates that recorded peak times are action potentials and not pre-synaptic artifacts by checking the second derivate of V(t).
    if d2(V(t)) > 7.5, the peak is not likely to be an action potential and instead is an artifact.
    '''
    d2V = np.diff(stim_voltage, n=2)
    stimulation_artifact_indices = stim_results[0]['stimulation_artifact_indices']
    peak_indices_to_filter = stim_results[0]['peak_indices']
    AP_half_widths = stim_results[0]['spike_width2']
    AP_widths = stim_results[0]['AP_width']
    peaks_to_delete = []
    
    '''
    Filter peaks to find peaks that are too close to one another. Action potentials
    have a refractory period that cannot be made smaller. If detected peaks are too close,
    apply filtering logic to determine which peak is false positive.
    '''
    for peak_index_num in range(len(peak_indices_to_filter)-1):
        prior_peak_index = peak_indices_to_filter[peak_index_num]
        current_peak_index = peak_indices_to_filter[peak_index_num+1]
        # filtering for peaks near stimulation artifacts
        nearest_artifact = int(find_nearest(stimulation_artifact_indices, current_peak_index))

        if prior_peak_index < 125:
            peaks_to_delete.append(prior_peak_index)
        
        # filtering for peaks near one another
        if abs(current_peak_index - prior_peak_index) < 75:
            nearest_artifact = int(find_nearest(stimulation_artifact_indices, current_peak_index))
            
            if current_peak_index-2 == nearest_artifact or current_peak_index == nearest_artifact or current_peak_index-3 == nearest_artifact or current_peak_index-11 == nearest_artifact :
                peaks_to_delete.append(current_peak_index)
            elif current_peak_index-2 == nearest_artifact or prior_peak_index-3 == nearest_artifact or prior_peak_index-11 == nearest_artifact:
                peaks_to_delete.append(prior_peak_index)
            elif current_peak_index > nearest_artifact > prior_peak_index:
                peaks_to_delete.append(current_peak_index)
            else:
                peaks_to_delete.append(prior_peak_index)
                
            
             
    '''Filter detected peaks based on 2nd derivative of V(t) to see if the
    peak is a presynaptic artifact. Presynaptic artifacts are approximatley 6 
    time units long and have similar characteristics. This filtering allows us to catch
    the vast majority of false positive action potentials while allowing few false
    negatives.'''
    for potential_AP_index_num in range(len(peak_indices_to_filter)):
        
        potential_AP_index = peak_indices_to_filter[potential_AP_index_num]
        d2v_start = peak_indices_to_filter[potential_AP_index_num]-3
        d2v_end = peak_indices_to_filter[potential_AP_index_num]+2
        d2v_check = d2V[d2v_start:d2v_end]
        
        # check for each of the three condtions for filtera
        try:
            if np.std(d2v_check) > 5 and (2 > AP_half_widths[potential_AP_index_num] > 6):
                incorrect_peak_index = potential_AP_index
                peaks_to_delete.append(incorrect_peak_index)
        except:
            if np.std(d2v_check) > 5 and (2 > AP_widths[potential_AP_index_num] > 6):
                incorrect_peak_index = potential_AP_index
                peaks_to_delete.append(incorrect_peak_index)
        try: 
            if np.std(d2v_check) > 5 and potential_AP_index-11 in stimulation_artifact_indices and (2 > AP_half_widths[potential_AP_index_num] > 6):
                incorrect_peak_index = potential_AP_index
                peaks_to_delete.append(incorrect_peak_index)
        except: 
            if np.std(d2v_check) > 5 and potential_AP_index-11 in stimulation_artifact_indices and (2 > AP_widths[potential_AP_index_num] > 6):
                incorrect_peak_index = potential_AP_index
                peaks_to_delete.append(incorrect_peak_index)
    
    # delete peaks detected to be false positives
    for incorrect_peak_index in peaks_to_delete:
        peak_indices_to_filter = np.delete(peak_indices_to_filter, np.where(peak_indices_to_filter == incorrect_peak_index))
    
    filtered_peaks = peak_indices_to_filter
    
    return filtered_peaks
    
def get_postynaptic_current_injection_times(DAC_data):
    '''
    Input: 1d numpy array of size

    Parameters
    ----------
    DAC_data : Tuple
        Tuple containing information re. DAC data from individual sweeps.
        [0] = time points at which the DAC signal (postsynaptic current injection) is not 0
        [1] = data type (float)

    Returns
    curr_times : List of Tuples 
        Each tuple contains the start and end indices for a given command stimulus. These start and end indices 
        indicate the indices in the sweep_times array at which the stimulus starts and ends, respectively, for
        a given stimulus)
    '''
    curr_times = []
    injection_times = np.nonzero(DAC_data)[0]
    total_curr_times = np.split(injection_times, 10)
    for curr_time in total_curr_times:
        start_time = curr_time[0]
        end_time = curr_time[-1]
        curr_times.append((start_time, end_time))
    return curr_times

def calculate_AP_frequencies_constant_current(stim_results):
    AP_indices = stim_results[0]['peak_indices']
    # if postsynaptic current injection is constant across stimulation
    AP_frequency = len(AP_indices)/0.2
    
    return np.array([AP_frequency])

def calculate_AP_frequencies_changing_current(stim_results):
    AP_indices = stim_results[0]['peak_indices']
    # calculate frequency of action potentials for 2 100 ms halves of each stim
    num_APs_first_100ms = len([peak for peak in AP_indices if 500 <= peak <1500])
    num_APs_second_100ms = len([peak for peak in AP_indices if 1500 <= peak <= 2500])
    
    AP_hz_first_100ms = num_APs_first_100ms/0.1
    AP_hz_second_100ms = num_APs_second_100ms/0.1
    
    return np.array([AP_hz_first_100ms]), np.array([AP_hz_second_100ms])

def re_arrange_dataframe(df):
    column_titles = ['sweep_num', 'stim_num', 'voltage_base', 
                     'AP_frequency','pre_pre_pairings', 'pre_post_pairings','pre_stim_indices','peak_indices']
    df = df.reindex(columns=column_titles)
    return df

def calculate_pre_pre_pairings(stim_results):
    '''
    Calculate the number of pre-pre pairings in a given stimulation
    '''
    stim_pre_pre_pairings = 0
    pre_indices = stim_results[0]['pre_stim_indices']  # get the indices representing the pre-synaptic inputs
    # get the indices representing the post-synaptic inputs
    post_indices = stim_results[0]['peak_indices']
    # calculate for initial timepoint
    
    # check for pre-pre pairings between a given pre-synaptic index and the pre-synaptic index in front of it
    for pre_index_num in range(len(pre_indices)-1):
        posts_between_pres = [post_index for post_index in post_indices if pre_indices[pre_index_num] < post_index < pre_indices[pre_index_num+1]]
        if not posts_between_pres:
            stim_pre_pre_pairings += 1
            
    return [stim_pre_pre_pairings]
        
    
def calculate_pre_post_offsets(stim_results):
    '''
    Caclulates the presynaptic-postsynpatic timing differences for a given stimulus interval in a trace.
    Also calculates the number of pre-pre pairings in a given trace

    '''
    pre_indices = stim_results[0]['pre_stim_indices']  # get the indices representing the pre-synaptic inputs
    # get the indices representing the post-synaptic inputs
    post_indices = stim_results[0]['peak_indices']
    pre_post_pairings = []
    # calculate pre-post pairings for first presynaptic index
    # to calcualte threshold for pre-post pairings
    pre_pre_pair_after = pre_indices[1] - pre_indices[0]
    pre_index = pre_indices[0]
    for post_index in post_indices:
        pre_post_pair = post_index - pre_index
        if pre_post_pair < pre_pre_pair_after:
            pre_post_pairings.append(pre_post_pair)
        else:
            break

    # calculate pre-post pairings for values between first and last
    for pre_num in range(1, len(pre_indices)-1):
        pre_index = pre_indices[pre_num]
        # pre-pre difference for pre index prior to current pre index
        pre_pre_pair_prior = pre_indices[pre_num-1] - pre_index
        # pre-pre difference for pre index after the current pre-index
        pre_pre_pair_after = pre_indices[pre_num+1] - pre_index
        for post_index in post_indices:
            pre_post_pair = post_index - pre_index
            if pre_pre_pair_prior < pre_post_pair < pre_pre_pair_after:
                pre_post_pairings.append(pre_post_pair)

    # calculate spike timing differences for last presynaptic index
    # to calcualte threshold for pre-post pairings
    pre_pre_pair_prior = pre_indices[-2] - pre_indices[-1]
    pre_index = pre_indices[-1]
    for post_index in post_indices:
        pre_post_pair = post_index - pre_index
        if pre_post_pair > pre_pre_pair_prior:
            pre_post_pairings.append(pre_post_pair)

    pre_post_pairings = np.array(pre_post_pairings)/10
    
    return pre_post_pairings


def extract_pre_post_timepoints(sweep_num, sweep_voltage, sweep_time, stim_num, injection_times, stim_frequencies, stim_EPSP_offset, stim_offset=500):
    # stim offset is 500 to put start of recording 50 ms before each stimulation (x interval for abf.sweepX is 0.1 ms, hence 500)

    # get the start and end indices for the stim number for this given trace
    stim_start, stim_end = injection_times[stim_num][0], injection_times[stim_num][1]
    # get the start and end time in ms for the stim (necessary for the event detection algorithm)
    stim_start_time, stim_end_time = sweep_time[stim_start], sweep_time[stim_end]
    
    # splice the original sweep data to only include the data for this stimulus, plus a 50 ms buffer prior for baseline calculation
    stim_voltage = sweep_voltage[stim_start-stim_offset:stim_end+stim_offset]
    stim_time = sweep_time[stim_start-stim_offset:stim_end+stim_offset]
    
    stim_trace = {}

    # Set the 'T' (=time) key of the trace (this needs to be in ms)
    stim_trace['T'] = stim_time

    # Set the 'V' (=voltage) key of the trace
    stim_trace['V'] = stim_voltage

    # Set the 'stim_start' (time at which a stimulus starts, in ms)
    # key of the trace
    # Warning: this need to be a list (with one element)
    stim_trace['stim_start'] = [stim_start_time]

    # Set the 'stim_end' (time at which a stimulus end) key of the trace
    # Warning: this need to be a list (with one element)
    stim_trace['stim_end'] = [stim_end_time]
    # Multiple traces can be passed to tshe eFEL at the same time, so the
    # argument should be a list
    traces = [stim_trace]

    #get average value of detected peaks to set threshold for AP detection
    AP_results = efel.getFeatureValues(traces,
                                         ['voltage_base', 'peak_indices', 'peak_voltage'])
    
    AP_Threshold = np.mean(AP_results[0]['peak_voltage'])-20

    efel.api.setThreshold(AP_Threshold)
    
    # Choose whichever features you are interested in extracting from the trace
    stim_results = efel.getFeatureValues(traces,
                                         ['voltage_base', 'peak_indices', 'min_AHP_indices', 'spike_width2', 'min_AHP_values', 'AHP_depth', 'AP_width'])
    
    
    
    
    stim_results[0]['sweep_num'] = np.array([sweep_num+1])
    stim_results[0]['stim_num'] = np.array([stim_num+1])
    
    # now extract presynaptic stimulation timepoints
    hz1, hz2 = stim_frequencies[0], stim_frequencies[1]
    # convert hz to index steps to acquire presynaptic stimulation indices
    hz1_step, hz2_step = 1/hz1 * 10000, 1/hz2 * 10000
    hz_1_indices = np.arange(500+stim_EPSP_offset, 1500+stim_EPSP_offset, hz1_step) #offset by delay between presynaptic stimulation and resulting EPSP
    hz_2_indices = np.arange(1500+stim_EPSP_offset, 2500+stim_EPSP_offset, hz2_step)
    pre_stim_indices = np.concatenate((hz_1_indices, hz_2_indices), axis=0)
    artifact_1_indices = np.arange(500, 1500, hz1_step) #offset by delay between presynaptic stimulation and resulting EPSP
    artifact_2_indices = np.arange(1500, 2500, hz2_step)
    stimulation_artifact_indices = np.concatenate((artifact_1_indices, artifact_2_indices), axis=0)
    stim_results[0]['stimulation_artifact_indices'] = stimulation_artifact_indices
    stim_results[0]['pre_stim_indices'] = pre_stim_indices
  
    
    
    # Filter out pressynaptic inputs mistaken as postsyanptic 
    filtered_peaks = filter_APs(stim_results, stim_time, stim_voltage, sweep_num, stim_num)
    stim_results[0]['peak_indices'] = filtered_peaks
    
        


    return stim_results


# set threshold for potential detection of action potentials to -15 mV
efel.api.setThreshold(-10)

# prompt user to 
root = tk.Tk()
root.withdraw()
data_path = filedialog.askdirectory() + "/"
data_files_list = os.listdir(data_path)
abf_files = [data_path + file_name for file_name in data_files_list if file_name.endswith(".abf")] # gather list of .abf files in folder to pull data from
cell_number = data_path.split("/")[-2]


# prompt user to enter the 2 frequencies for the induction stimulation
hz1 = tk.simpledialog.askinteger("Frequency 1", "Please enter the first frequency used during induction")
hz2 = tk.simpledialog.askinteger("Frequency 2", "Please enter the second frequency used during induction")
stim_EPSP_offset = tk.simpledialog.askfloat("Presynaptic Offset", "Please enter the presynaptic offset (can be as precise as 0.1 ms)") * 10
stim_frequencies = (hz1, hz2)
np.set_printoptions(threshold=np.inf)

stim_dfs = []
sweep_num = 0

for abf_file in abf_files:
    abf = pyabf.ABF(abf_file)
    
    for Sweep_Number in range(len(abf.sweepList)):
        abf.setSweep(Sweep_Number)
        sweep_time = abf.sweepX  # time data units = (s)
        sweep_time *= 1000           # convert time data to units (ms)
        sweep_voltage = abf.sweepY   # voltage data, units = (mV)
        sweep_command = abf.sweepC  # stimulus waveform data, units = (pA)
        # get the start and end points for all current injections in the sweep
        current_injection_times = get_postynaptic_current_injection_times(sweep_command)
        # finds total number of stimulations
        num_stims = len(current_injection_times)
        
        for Stimulation_Number in range(num_stims):
            stim_data = extract_pre_post_timepoints(
                Sweep_Number, sweep_voltage, sweep_time, Stimulation_Number, current_injection_times, stim_frequencies, stim_EPSP_offset)
            AP_frequency = calculate_AP_frequencies_constant_current(stim_data) # calculate the frequency of postsynaptic firing during current injection
            pre_post_pairings = calculate_pre_post_offsets(stim_data)           # calculate the timing differences between pre and postsynaptic events for each period of current injection
            pre_pre_pairings = calculate_pre_pre_pairings(stim_data)
            stim_data[0]['pre_post_pairings'], stim_data[0]['AP_frequency'] = pre_post_pairings, AP_frequency
            stim_data[0]['pre_pre_pairings'] = pre_pre_pairings
            hope = plotting_accuracy_test(Sweep_Number, sweep_voltage, sweep_time, Stimulation_Number, current_injection_times, stim_data)
            
            stim_df = pd.DataFrame.from_dict(data=stim_data[0], orient='index')
            stim_df = stim_df.T # transpose to put columns on top
            stim_dfs.append(stim_df)
        sweep_num += 1
    del abf
        
trace_df = pd.concat(stim_dfs, axis=0)
trace_df = re_arrange_dataframe(trace_df)
trace_df.to_csv('%s%s_STDP_Analysis.csv' % (data_path, cell_number), index=True,header=True, encoding='utf-8')

trace_df.hist(column='pre_post_pairings', bins=25, grid=False, figsize=(12,8), color='#86bf91', zorder=2, rwidth=0.9)

