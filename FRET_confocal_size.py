#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 18:30:37 2021

@author: Mathew
"""


import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
import os
from matplotlib.colors import LogNorm
pathlist=[]

# Where to store the overall file containing means etc. for each experiment.
path_root=r"/Users/Mathew/Documents/Current analysis/20230718_lysate/Timepoints/"

# # Foldert to analyse here:
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20230718_lysate/Sample10/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20230718_lysate/Sample11/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20230718_lysate/Sample12/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20230718_lysate/Sample13/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20230718_lysate/Sample14/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20230718_lysate/Sample15/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20230718_lysate/Sample16/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20230718_lysate/Sample17/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20230718_lysate/Sample18/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20230718_lysate/Sample19/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20230718_lysate/Sample20/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20230718_lysate/Sample21/")
# pathlist.append(r"/Users/Mathew/Documents/Current analysis/20230718_lysate/Sample22/")
# pathlist.append(r"/Users/Mathew/Documents/Current analysis/20230718_lysate/Sample23/")
# pathlist.append(r"/Users/Mathew/Documents/Current analysis/20230718_lysate/Sample24/")
# pathlist.append(r"/Users/Mathew/Documents/Current analysis/20230718_lysate/Sample25/")
# pathlist.append(r"/Users/Mathew/Documents/Current analysis/20230718_lysate/Sample26/")
# pathlist.append(r"/Users/Mathew/Documents/Current analysis/20230718_lysate/Sample27/")
# pathlist.append(r"/Users/Mathew/Documents/Current analysis/20230718_lysate/Sample28/")
# pathlist.append(r"/Users/Mathew/Documents/Current analysis/20230718_lysate/Sample29/")
 

file_stem="AS"          # This is the part of the filename that will be searched for in each folder.

# Thresholds and other parameters:
    
channelA_thresh=8      # Threshold for Channel A (Green).
channelB_thresh=8 # Threshold for Channel B (Red).
channelA_AF=0.28       # Autofluorescence
channelB_AF=0.65
xtalk=0.03            # Cross-talk from A to B
size_split_med=5
size_split_lar=150

def load_files(filename_contains,path):
    print(path)
    num=0
    channelA_sample=[]             # Where channel A data will be stored
    channelB_sample=[]             # Where channel B data will be stored
    for root, dirs, files in os.walk(path):
      for name in files:
              # print(name)
              if filename_contains in name:
                  if 'pdf' not in name:
                      resultsname = name
                      # print(name)
                      num+=1
                      a=0
                      with open(path+name) as csvDataFile:                                                # Opens the file as a CSV
                            csvReader = csv.reader(csvDataFile,delimiter='\t')                           # Assigns the loaded CSV file to csvReader. 
                            for row in csvReader:
                                channelA_sample.append(row[0])                                                     # For every row in in csvReader, the values are apended to green and red.         
                                channelB_sample.append(row[1])
                                a+=1
            
                            # print ("Loaded %s, which contains %s rows."%(resultsname,a))
    rows=len(channelA_sample)
    # print("Loaded %s files in total, with a total of %s rows"%(num,rows))
        
    channelA_arr_sample=np.asarray(channelA_sample,dtype=np.float32)                              # Converts these to numpy arrays for vector calcs.
    channelB_arr_sample=np.asarray(channelB_sample,dtype=np.float32)
    return channelA_arr_sample,channelB_arr_sample,num

def maxQ():
    q_vals = np.zeros(shape=(20,20))

    for A in range(20):
        for B in range(20):
            channelA_only_events=channelA_arr[(channelA_arr>A)]                       # Total A events
            channelB_only_events=channelB_arr[(channelB_arr>B)]                       # Total B events
            channelA_events=channelA_arr[np.logical_and(channelA_arr>A, channelB_arr>B)]  # A coincident events             
            
            
            # Now need to account for chance events:
            
            channelB_shuffle=channelB_arr.copy()
            np.random.shuffle(channelB_shuffle)
            
            channelA_chance=channelA_arr[np.logical_and(channelA_arr>A, channelB_shuffle>B)]    # These are the chance events    
            
            # Now need to calculate Q value:
            
            var_real_events=float(len(channelA_events))
            var_A_events=float(len(channelA_only_events))
            var_B_events=float(len(channelB_only_events))
            var_chance_events=float(len(channelA_chance))
            Q=float((var_real_events-var_chance_events)/(var_A_events+var_B_events-(var_real_events-var_chance_events)))
    
            q_vals[A][B]=Q
            

    maximum_Q=np.amax(q_vals)
    result=np.where(q_vals == np.amax(q_vals))
    ThresholdA,ThresholdB=result
    
    print('The maximum value of Q is %.3f, with a threshold of %s in channel A, and %s in channel B.'%(maximum_Q,str(ThresholdA),ThresholdB))
    
    
    contourplot = plt.contourf(q_vals,20,origin='lower')
    cbar = plt.colorbar(contourplot)
    plt.xlabel("Channel B Threshold")
    plt.ylabel("Channel A Threshold")
    cbar.ax.set_ylabel('Q')
    


Output_all = pd.DataFrame(columns=['Path','Number_of_files','Threshold_A','Threshold_B','Events_A','Events_B','Events_coincindent',
                                       'Events_chance','Q','Total_Intensity_mean','Total_Intensity_SD','Total_Intensity_med','Intensity_A_mean','Intensity_A_SD','Intensity_A_med','Intensity_B_mean','Intensity_B_SD','Intensity_B_med',
                                       'Sizes_mean','Sizes_SD','Sizes_med','A_ave','B_ave','small','medium','large','monomer_in_agg'])
All_FRET_Small=pd.DataFrame()
All_FRET_Medium=pd.DataFrame()
All_FRET_Large=pd.DataFrame()

for path in pathlist:
    # path=path.replace('/Ab/','/ThT/')
    print(path)
    channelA_arr,channelB_arr,num=load_files(file_stem,path)
    
    
    
    # Now need to account for autofluorescence and crosstalk etc. 
    
    channelB_arr=(channelB_arr-xtalk*channelA_arr)-channelB_AF
    channelA_arr=channelA_arr-channelA_AF
    
    #This part is for the thresholding:
    
    channelA_only_events=channelA_arr[(channelA_arr>channelA_thresh)]                       # Total A events
    channelB_only_events=channelB_arr[(channelB_arr>channelB_thresh)]                       # Total B events
    channelA_events=channelA_arr[np.logical_and(channelA_arr>channelA_thresh, channelB_arr>channelB_thresh)]  # A coincident events             
    channelB_events=channelB_arr[np.logical_and(channelA_arr>channelA_thresh, channelB_arr>channelB_thresh)]  # B coincident events
    
    
   
    
    channelA_only_minus_events=channelA_arr[np.logical_and(channelA_arr>channelA_thresh, channelB_arr<channelB_thresh)] # Non-oincidence events
    channelB_only_minus_events=channelB_arr[np.logical_and(channelA_arr<channelA_thresh, channelB_arr>channelB_thresh)] # Non-oincidence events
    
    channelA_brightness=channelA_only_minus_events.mean()
    channelB_brightness=channelB_only_minus_events.mean()
    
    
    channelA_mean=channelA_events.mean()
    channelA_SD=channelA_events.std()
    channelA_med=np.median(channelA_events)
    
    channelB_mean=channelB_events.mean()
    channelB_SD=channelB_events.std()
    channelB_med=np.median(channelB_events)
    # Now need to account for chance events:
    
    channelB_shuffle=channelB_arr.copy()
    np.random.shuffle(channelB_shuffle)
    
    channelA_chance=channelA_arr[np.logical_and(channelA_arr>channelA_thresh, channelB_shuffle>channelB_thresh)]    # These are the chance events    
    channelB_chance=channelB_shuffle[np.logical_and(channelA_arr>channelA_thresh, channelB_shuffle>channelB_thresh)]
    
    # Now need to calculate Q value:
    
    var_real_events=float(len(channelA_events))
    var_A_events=float(len(channelA_only_events))
    var_B_events=float(len(channelB_only_events))
    var_chance_events=float(len(channelA_chance))
    Q=float((var_real_events-var_chance_events)/(var_A_events+var_B_events-(var_real_events-var_chance_events)))
    
    print(('There were %s A events, %s B events, %s coincidence events, and %s chance events. Q = %f.')%(var_A_events,var_B_events,var_real_events,var_chance_events,Q))
    
    
    # Now want histograms etc. 
    
    ln_events=np.log(channelB_events/channelA_events)
    ln_chance=np.log(channelB_chance/channelA_chance)
    
    FRET_events=channelB_events/(channelB_events+channelA_events)
    
    FRET_hist=np.histogram(FRET_events, bins=20, range=(0,1), normed=None, weights=None, density=None)
    FRET_hist_norm=FRET_hist[0][:]/len(channelA_only_events)
    
    # FRET_all[path]=FRET_hist_norm
    # FRET_all.to_csv(path_root + '/' + 'all_FRET.csv', sep = '\t')
    
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = "12"
    plt.figure(figsize=(8, 6))
    plt.hist(FRET_events, bins = 20,range=[0,1], rwidth=0.9,ec='black',color='#ff0000',alpha=0.8,label="Real Events")  
    plt.xlabel('Proximity Ratio')
    plt.ylabel('Number of events')
    plt.savefig(path+'/'+'FRET.pdf')
    plt.show()
    
    channelAout=[]
    channelAout.append(var_A_events)
    channelBout=[]
    channelBout.append(var_B_events)
    Qout=[]
    Qout.append(Q)
  
    output_stats=pd.DataFrame()
    output_stats['Donor']=channelAout
    output_stats['Acceptor']=channelBout
    output_stats['Q']=Qout
    output_stats.to_csv(path + '/' + 'Stats.csv', sep = '\t')
    
    
    
    raw_output=pd.DataFrame()
    
    raw_output['Donor']=channelA_events 
    raw_output['Acceptor']=channelB_events
    raw_output['FRET']=FRET_events
    raw_output['Z']=ln_events   
    
    raw_output.to_csv(path + '/' + 'Raw.csv', sep = '\t')
    
    textstr='Q = %.3f'%Q
    
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = "12"
    plt.figure(figsize=(8, 6))
    plt.hist(FRET_events, bins = 20,range=[0,1], rwidth=0.9,ec='black',color='#ff0000',alpha=0.8)
   
    plt.text(0.05,0.90, textstr,transform=plt.gca().transAxes)
    plt.legend(loc='upper right')         
    plt.xlabel('Proximity Ratio')
    plt.ylabel('Number of events')
    plt.savefig(path+'/'+'FRET.pdf')
    plt.show()
    
    
    channelB_arr_inv=channelB_arr*(-1)
    
    plt.plot(channelA_arr,color='green')
    plt.plot(channelB_arr_inv,color='red')
    plt.xlabel('Bin number')
    plt.ylabel('Intensity (photons)')
    plt.xlim(0,500)
    plt.ylim(-50,50)
    plt.savefig(path+'/'+'example_trace.pdf')
    plt.show()
    
    
    
    total_intensity=channelB_events+channelA_events
    plt.hist(total_intensity, bins = 60,range=[0,1000], rwidth=0.9,ec='black',color='#ff0000',alpha=0.8,)
    plt.yscale('log')
    plt.xlabel('Total intensity (photons)')
    plt.ylabel('Number of events')
    plt.savefig(path+'/'+'Intensity.pdf')
    plt.show()
    
    total_mean=total_intensity.mean()
    total_SD=total_intensity.std()
    total_med=np.median(total_intensity)
    
    sizes=2*(channelB_events+channelA_events)/channelA_brightness
    plt.hist(sizes, bins = 20,range=[0,50], rwidth=0.9,ec='black',color='#ff0000',alpha=0.8,)
    plt.yscale('log')
    plt.xlabel('Approximate size (monomers)')
    plt.ylabel('Number of events')
    plt.savefig(path+'/'+'Sizes.pdf')
    plt.show()
    
    events_monomer=sum(sizes)
    
    
    
    sizes_small=sizes[(sizes<size_split_med)] 
    sizes_medium=sizes[np.logical_and(sizes<size_split_lar,sizes>size_split_med)]
    sizes_large=sizes[(sizes>size_split_lar)] 
    
    
    # These are plots as in all of aSyn papers:
    
    FRET_small=FRET_events[(sizes<size_split_med)] 
    FRET_medium=FRET_events[np.logical_and(sizes<size_split_lar,sizes>size_split_med)]
    FRET_large=FRET_events[(sizes>size_split_lar)] 
   
    
    FRET_small_hist=[]
    
  
    j=0
    for i in range(0,20):
        
        temp=FRET_small[np.logical_and(FRET_small<j+0.05,FRET_small>j)]
        
        counter=len(temp)
        FRET_small_hist.append(counter)
        j+=0.05
    
    All_FRET_Small[path]=FRET_small_hist
    All_FRET_Small.to_csv(path_root + '/' + 'All_FRET_Small.csv', sep = '\t')
    FRET_medium_hist=[]
    

    j=0
    for i in range(0,20):
        
        temp=FRET_medium[np.logical_and(FRET_medium<j+0.05,FRET_medium>j)]
        
        counter=len(temp)
        FRET_medium_hist.append(counter)
        j+=0.05
    
    All_FRET_Medium[path]=FRET_medium_hist
    All_FRET_Medium.to_csv(path_root + '/' + 'All_FRET_Medium.csv', sep = '\t')
    
    FRET_large_hist=[]
    j=0
    for i in range(0,20):
        
        temp=FRET_large[np.logical_and(FRET_large<j+0.05,FRET_large>j)]
        
        counter=len(temp)
        FRET_large_hist.append(counter)
        
        
        
        j+=0.05
    All_FRET_Large[path]=FRET_large_hist
    
    All_FRET_Large.to_csv(path_root + '/' + 'All_FRET_Large.csv', sep = '\t')
   
    plt.hist(FRET_small, bins = 20,range=[0,1], width=0.05,facecolor='white',edgecolor='#949494')
    plt.xlabel('FRET Efficiency')
    plt.ylabel('Number of events')
    plt.savefig(path+'/'+'FRET_Small.pdf')
    plt.show()
    
        
    plt.hist(FRET_medium, bins = 20,range=[0,1], width=0.05,facecolor='white',edgecolor='#949494')
    plt.xlabel('FRET Efficiency')
    plt.ylabel('Number of events')
    plt.savefig(path+'/'+'FRET_Small.pdf')
    plt.show()
    
        
    plt.hist(FRET_large, bins = 20,range=[0,1], width=0.05,facecolor='white',edgecolor='#949494')
    plt.xlabel('FRET Efficiency')
    plt.ylabel('Number of events')
    plt.savefig(path+'/'+'FRET_Small.pdf')
    plt.show()

    
    small_num=len(sizes_small)
    med_num=len(sizes_medium)
    large_num=len(sizes_large)
    
    sizes_mean=sizes.mean()
    sizes_SD=sizes.std()
    sizes_med=np.median(sizes)
    Output_all= Output_all.append({'Path':path,'Number_of_files':num,'Threshold_A':channelA_thresh,'Threshold_B':channelB_thresh,'Events_A':var_A_events,'Events_B':var_B_events,'Events_coincindent':var_real_events,'Q':Q,
                                       'Events_chance':var_chance_events,'Total_Intensity_mean':total_mean,'Total_Intensity_SD':total_SD,'Total_Intensity_med':total_med,'Intensity_A_mean':channelA_mean,'Intensity_A_SD':channelA_SD,'Intensity_A_med':channelA_med,'Intensity_B_mean':channelB_mean,'Intensity_B_SD':channelB_SD,'Intensity_B_med':channelB_med,'Sizes_mean':sizes_mean,'Sizes_SD':sizes_SD,'Sizes_med':sizes_med,'A_ave':channelA_brightness,'B_ave':channelB_brightness,
                                       'small':small_num,'medium':med_num,'large':large_num,'monomer_in_agg':events_monomer},ignore_index=True)

    Output_all.to_csv(path_root + '/' + 'all_metrics.csv', sep = '\t')
    
    output_sizes=pd.DataFrame()
    

    raw_output['FRET']=FRET_events
    raw_output['Size']=sizes  
    
    twod=np.zeros((20,100),dtype=float)
    
    for i,j in zip(FRET_events,sizes):
        FRET=round(i*20)
        # print(FRET)
        Size=round(j)
        if(Size<100):
            if(FRET<20):
                # print(FRET)
                # print(Size)
                twod[FRET,Size]+=1
    np.savetxt(path + '/' + '2D.csv',twod)
    raw_output.to_csv(path + '/' + 'FRET_and_size.csv', sep = '\t')

    x_bins = np.linspace(0,1,20)
    y_bins = np.linspace(0,50,50)
    plt.hist2d(FRET_events, sizes, bins=[x_bins,y_bins],cmap='jet')
    # plt.colors.LogNorm()
    plt.colorbar()
    plt.savefig(path+'/'+'2D_Sizes.pdf')
    


    
# runfile('/Users/Mathew/Documents/Edinburgh Code/Global Fitting/2Gaussian.py', wdir='/Users/Mathew/Documents/Edinburgh Code/Global Fitting')