import numpy as np
import tensorflow as tf
from math import e


class Complex_SystemModel(object):
    def __init__(self, scenario, k, doa= x_dire, freq_values = None):
        # self.A = A
        # self.A_T = A.conj().T  
        self.scenario = scenario                                    # Narrowband or Broadband
        # self.n = A.size()[0]                                      # Number of sensors in element
        # self.m = A.size()[1]                                      # Number of grid sizes    
        self.n = n                                           
        self.k = k                                                  # Number of sources
        self.Scenario_define(freq_values)                           # Define parameters    
        self.Create_array()                                         # Define array indicies 
             
    def Scenario_define(self, freq_values):
        if self.scenario == "Broadband":
            ## frequencies initialization ##
            self.min_freq = freq_values[0]                          # Define minimal frequency value  
            self.max_freq = freq_values[1]                          # Define maximal frequency value  
            self.f_rng = np.linspace(start= self.min_freq,
                             stop= self.max_freq,
                              num= self.max_freq - self.min_freq,
                               endpoint = False)                    # Frequency range of interest  
            self.f_sampling = 2 * (self.max_freq)                   # Define sampling rate as twice the maximal frequency
            self.time_axis = np.linspace(0, 1, self.f_sampling,
                             endpoint = False)                      # Define time axis
            ## Array initialization ##
            self.dist = 1 / (2 * self.max_freq)                     # distance between array elements
        else: 
            ## frequencies initialization ##
            self.min_freq = None
            self.max_freq = None
            self.f_rng = None
            self.fs = None
            ## Array initialization ##
            self.dist = 1 / 2                                       # distance between array elements

    def Create_array(self):
        self.array = np.linspace(0, self.n, self.n, endpoint = False)   # create array of sensors locations
    

    def __str__(self):
        print("System Model Summery:")
        for key,value in self.__dict__.items():
            print (key, " = " ,value)
        return "End of Model"
