#!/c/Users/Kendall/Miniconda2/python


import numpy as np
import pandas as pd
import sys

class PhysicalSystem:
    
    def __init__(self, diam_tank):
        self.diam_tank  = diam_tank
        #self.crit_level = crit_level
        #self.diam_valve = diam_valve
        #self.grav       = grav

    def echo_system_params(self):
        print 'The tank diameter is: ',self.diam_tank
        #print 'The critical level is:',self.crit_level
        #print 'The valve diameter is:',self.diam_valve

    def calculate_TW(self):    
        self.TW = self.diam_tank * 100
        return self.TW




if __name__ == "__main__":
    
    # read in the input file
    inputfileID = sys.argv[1]
    #inputfileID = open(r'input_for_pythoncode.i','r')
    inputData = np.loadtxt(inputfileID)
    #inputfileID.close()
    #print inputData
    
    # define constant system parameters
    system = PhysicalSystem(inputData)
    system.echo_system_params()

    # do calculations
    TW = system.calculate_TW()
    #print TW
    
    # create a matrix or array of the output data to write to csv
    outputData = np.zeros([1,2])
    outputData[0,0] = inputData
    outputData[0,1] = TW
    
    # write to output csv
    HEADER = ['diam', 'time_window']
    dataframe = pd.DataFrame(outputData)
    dataframe.to_csv("finalizedCodeOutput.csv",header=HEADER,index=False)
    
    
    
