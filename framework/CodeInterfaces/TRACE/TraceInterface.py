# Copyright 2017 Battelle Energy Alliance, LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
  Created on November 08, 2018
  @created: Kendall  Boniface
"""

from __future__ import division, print_function, absolute_import
import warnings
warnings.simplefilter('default',DeprecationWarning)

from GenericCodeInterface import GenericCode
import numpy as np
import csvUtilities as csvU
import dynamicEventTreeUtilities as detU
import csv
import glob
import os
import copy
import re
import math
import sys
import xtvReader
import pandas as pd

class Trace(GenericCode):
  """
    Trace Code Interface
  """
  def generateCommand(self, inputFiles, executable, clargs=None, fargs=None, preExec=None):
    """
      See base class.  Collects all the clargs and the executable to produce the command-line call.
      Returns tuple of commands and base file name for run.
      Commands are a list of tuples, indicating parallel/serial and the execution command to use.
      @ In, inputFiles, list, List of input files (lenght of the list depends on the number of inputs have been added in the Step is running this code)
      @ In, executable, string, executable name with absolute path (e.g. /home/path_to_executable/code.exe)
      @ In, clargs, dict, optional, dictionary containing the command-line flags the user can specify in the input (e.g. under the node < Code >< clargstype =0 input0arg =0 i0extension =0 .inp0/ >< /Code >)
      @ In, fargs, dict, optional, a dictionary containing the axuiliary input file variables the user can specify in the input (e.g. under the node < Code >< clargstype =0 input0arg =0 aux0extension =0 .aux0/ >< /Code >)
      @ In, preExec, string, optional, a string the command that needs to be pre-executed before the actual command here defined
      @ Out, returnCommand, tuple, tuple containing the generated command. returnCommand[0] is the command to run the code (string), returnCommand[1] is the name of the output root
    """
    # WILL NEED TO CONSIDER RESTART FILES HERE "trcrst or whatevername.rst"
    # find the input file (check that one input is provided)
    # create output file root
    #outputfile = 'out~' + inputFiles[0].getBase()
    executeCommand = [('parallel', executable +' -p '+ inputFiles[0].getFilename())]
    returnCommand = executeCommand, inputFiles[0].getBase() #output file root is same as input file
    return returnCommand

  def _readMoreXML(self,xmlNode):
    """
      Function to read the portion of the xml input that belongs to this specialized class and initialize
      some members based on inputs. This can be overloaded in specialize code interface in order to
      read specific flags
      @ In, xmlNode, xml.etree.ElementTree.Element, Xml element node
      @ Out, None.
    """
    GenericCode._readMoreXML(self,xmlNode)
    # add things here you want to read in the RAVEN XML file (in addition to the standard Generic
    # Code interface stuff)

  def createNewInput(self,currentInputFiles,oriInputFiles,samplerType,**Kwargs):
    """
      This method is used to generate an input based on the information passed in.
      @ In, currentInputFiles, list,  list of current input files (input files from last this method call)
      @ In, oriInputFiles, list, list of the original input files
      @ In, samplerType, string, Sampler type (e.g. MonteCarlo, Adaptive, etc. see manual Samplers section)
      @ In, Kwargs, dictionary, kwarded dictionary of parameters. In this dictionary there is another dictionary called "SampledVars"
            where RAVEN stores the variables that got sampled (e.g. Kwargs['SampledVars'] => {'var1':10,'var2':40})
      @ Out, newInputFiles, list, list of newer input files, list of the new input files (modified and not)
    """
    if 'dynamiceventtree' in str(samplerType).lower():
      # here we will be doing the stuff for the DET
      raise RunTimeError('DET not yet here')
    return GenericCode.createNewInput(self,currentInputFiles,oriInputFiles,samplerType,**Kwargs)

#  def checkForOutputFailure(self,output,workingDir):
#    """
#      This method is called by the RAVEN code at the end of each run  if the return code is == 0.
#      This method needs to be implemented by the codes that, if the run fails, return a return code that is 0
#      This can happen in those codes that record the failure of the job (e.g. not converged, etc.) as normal termination (returncode == 0)
#      This method can be used, for example, to parse the outputfile looking for a special keyword that testifies that a particular job got failed
#      (e.g. in RELAP5 would be the keyword "********")
#      @ In, output, string, the Output name root
#      @ In, workingDir, string, current working dir
#      @ Out, failure, bool, True if the job is failed, False otherwise
#    """
#    failure = True
#    goodWord  = ["Transient terminated by end of time step cards","Transient terminated by trip"]
#    try:
#      outputToRead = open(os.path.join(workingDir,output+'.o'),"r")
#    except:
#      return failure
#    readLines = outputToRead.readlines()
#
#    for goodMsg in goodWord:
#      if any(goodMsg in x for x in readLines[-20:]):
#        failure = False
#    return failure
    
  def finalizeCodeOutput(self, command, output, workingDir):
    """
      finalizeCodeOutput checks TRACE csv files and looks for iEvents and
      continous variables we specified in < boolMaapOutputVariables> and
      contMaapOutputVairables> sections of RAVEN_INPUT.xml file. Both
      < boolMaapOutputVariables> and <contMaapOutputVairables> should be
      contained into csv MAAP csv file
      In case of DET sampler, if a new branching condition is met, the
      method writes the xml for creating the two new branches.
      @ In, command, string, the command used to run the just ended job
      @ In, output, string, the Output name root
      @ In, workingDir, string, current working dir
      @ Out, output, string, output csv file containing the variables of interest specified in the input
    """
    # the following is taken from the relap code interface
    #outfile = os.path.join(workingDir,output+'.o')
    #outputobj=relapdata.relapdata(outfile,self.outputDeck)
    #if outputobj.hasAtLeastMinorData():
    #  outputobj.writeCSV(os.path.join(workingDir,output+'.csv'))
    #else:
    #  raise IOError('Relap5 output file '+ command.split('-o')[0].split('-i')[-1].strip()+'.o' + ' does not contain any minor edits. It might be crashed!')
    
    xtv_file = os.path.join(workingDir,output+'.xtv')
    # define which time-dependent values to look at
    trace_vars = ['pn-500A01','tln-500A01'] # COME BACK AND MAKE THIS MORE GENERAL
    with open(xtv_file,'rb') as xtvFileHandle:
        xtvObj = xtvReader.XtvFile(xtvFileHandle, verbose=True)
        # extract one variable to know how many time steps need to be saved
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
        
        # write to csv
        df = pd.DataFrame(data)
        head = trace_vars
        head.insert(0,'time')
        df.to_csv(os.path.join(workingDir,'finalizedCodeOutput.csv'),header=head,index=False)
    
    return "finalizedCodeOutput"




    