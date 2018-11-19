import numpy as np

def evaluate(self):
  return np.log(self.var_given_OAmean) - (self.var_OAsigma*self.var_OAsigma)/2.0