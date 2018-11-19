# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:01:58 2018

@author: Kendall
"""


import xtvReader
import pandas as pd
import numpy as np

# the xtv file to read
xtv_file = "spent_fuel_pool.xtv"
# define which time-dependent values to look at
# COME BACK AND MAKE THIS MORE GENERAL
trace_vars = ['pn-500A01','tln-500A01']


with open(xtv_file,'rb') as xtvFileHandle:
    xtvObj = xtvReader.XtvFile(xtvFileHandle, verbose=True)
    # gets a single value at a single point in time
    ##value_dataChannel = xtvObj.getDataChannel(2.0,'pn-500A01')
    # get values from several points in time and give it more basic component info
    ##time_points = [0.0, 5.0]
    ##value_timeData = xtvObj.getTimeData(time_points, 500, 'pipe', 'pn', 1)
   
    # extrace one variable to know how many time steps need to be saved
    vector = xtvObj.getTimeVector(trace_vars[0])
    time,values = zip(*vector)
    nvars = len(trace_vars)
    ntime = len(vector)
    # initialize an array for writing to a final csv
    data = np.zeros([ntime,nvars+1])
    data[:,0] = time
    # now loop through the variables and store them in "data"
    for v,var in enumerate(trace_vars):
        vector = xtvObj.getTimeVector(trace_vars[v])
        time,values = zip(*vector)
        data[:,v+1] = values
    
    # done getting things from the xtv file
    xtvFileHandle.close()    
    
    # now write to csv
    df = pd.DataFrame(data)
    head = trace_vars
    head.insert(0,'time')
    df.to_csv("finalizedCodeOutput.csv",header=head,index=False)
  
