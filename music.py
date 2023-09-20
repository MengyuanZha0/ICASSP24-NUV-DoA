from scipy.linalg.matfuncs import absolute
##############################################
########## MUSIC & Root-MUSIC ################
##############################################

from System_Model import *
from operator import itemgetter
import numpy as np
import scipy as sc
import scipy.signal
from useful_func import *
import matplotlib.pyplot as plt
import torch
from scipy.signal import find_peaks, peak_prominences

class Model_Based_methods(object):
    def __init__(self, System_model, vectorized = True):
        self.angels = np.linspace(-1 * np.pi / 2, np.pi / 2, 360, endpoint=False)                        # angle axis for represantation of the MUSIC spectrum
        self.system_model = System_model
        self.dist = System_model.dist
        self.k = System_model.k   #number of sources
        self.n = System_model.n   #number of sensors
        self.m = m
        self.dist = System_model.dist
        self.vectorized = vectorized

    #############################################################################

    def Classic_MUSIC(self, y, NUM_OF_SOURCES, SPS=False, sub_array_size=0):
        '''
        Implementation of the model-based MUSIC algorithm in Narrow-band scenario.
        
        Input:
        @ y = samples vector shape : Nxsn
        @ NUM_OF_SOURCES = known number of sources flag
        @ SPS = pre-processing Spatial smoothing algorithm flag
        @ sub_array_size = size of sub array for spatial smoothing
        
        Output:
        @ DOA_pred = the predicted DOA's
        @ Spectrum = the MUSIC spectrum
        @ M = number of estimated/given sources
        
        '''
        y = np.array(y)

        # if NUM_OF_SOURCES:                                                      # NUM_OF_SOURCES = TRUE : number of sources is given 
        #     k = self.k                                                
        # else:                                                                   # NUM_OF_SOURCES = False : M is given using  multiplicity of eigenvalues
        #     # clustring technique                                   
        #     pass
        if SPS:
            number_of_sensors = self.n
            number_of_sub_arrays = number_of_sensors - sub_array_size + 1
            
            ## Averaged covariance matrix
            R_y = np.zeros((sub_array_size, sub_array_size)) + 1j * np.zeros((sub_array_size, sub_array_size))
            for j in range(number_of_sub_arrays):
                y_sub = y[j:j + sub_array_size,:]
                R_y += np.cov(y_sub)
            R_y /= number_of_sub_arrays
            # print('1')
            # R_y = np.mean(np.array(Rx_sub_arrays), 0)
        else:
            R_y = np.cov(y)                                                     # Create covariance matrix from samples
            # print(R_y.shape)

        eigenvalues, eigenvectors = np.linalg.eig(R_y)                          # Find the eigenvalues and eigenvectors using EVD
        Un = eigenvectors[:, k:]                                                # Take only the eigenvectors associated with Noise subspace 
        # print("size of Un is{}".format(Un.shape))
        Spectrum,_= self.spectrum_calculation(Un)
        DOA_pred, _ = scipy.signal.find_peaks(Spectrum)                         # Find maximal values in spectrum
        # print(Spectrum)
        # print(type(Spectrum))
        # print(Spectrum.shape)
        DOA_pred = list(DOA_pred)
        DOA_pred.sort(key = lambda x: Spectrum[x], reverse = True)
        return DOA_pred, Spectrum, k                                                    # return estimated DOA

    #############################################################################

    def Classic_Root_MUSIC(self, y, NUM_OF_SOURCES, SPS=False, sub_array_size=0):
        '''
        Implementation of the model-based Root-MUSIC algorithm in Narrow-band scenario.
        
        Input:
        @ X = samples vector shape : Nx1
        @ NUM_OF_SOURCES = known number of sources flag
        @ SPS = pre-processing Spatial smoothing algorithm flag
        @ sub_array_size = size of sub array for spatial smoothing
        
        Output:
        @ DOA_pred = the predicted DoA's
        @ roots = the roots of true DoA's 
        @ M = number of estimated/given sources
        @ roots_angels_all = all the roots produced by the algorithm
        @ DOA_pred_all = all the angels produced by the algorithm 
      
        '''
        # print(y.shape)
        y = np.array(y)
        
        if NUM_OF_SOURCES:                                                              # NUM_OF_SOURCES = TRUE : number of sources is given
            k = self.k
        else:                                                                   # NUM_OF_SOURCES = False : M is given using  multiplicity of eigenvalues        
            # clustring technique
            pass
        
        if SPS:
            number_of_sensors = self.n
            number_of_sub_arrays = number_of_sensors - sub_array_size + 1
            
            ## Averaged covariance matrix
            R_y = np.zeros((sub_array_size, sub_array_size)) + 1j * np.zeros((sub_array_size, sub_array_size))
            for j in range(number_of_sub_arrays):
                y_sub = y[j:j + sub_array_size,:]
                R_y += np.cov(y_sub)
            R_y /= number_of_sub_arrays
            # print('1')
            # R_x = np.mean(np.array(Rx_sub_arrays), 0)
        else:
            R_y = np.cov(y)  
            # print(R_y.shape)                                                    # Create covariance matrix from samples
        
        eigenvalues, eigenvectors = np.linalg.eig(R_y)                          # Find the eigenvalues and eigenvectors using EVD
        Us, Un = eigenvectors[:, :k] , eigenvectors[:, k:]
        F = Un @ np.conj(Un).T                                                  # Set F as the matrix contains Information
        coeff = sum_of_diag(F)                                                  # Calculate the sum of the diagonals of F
        roots = list(find_roots(coeff))                                         
                                                                                # By setting the diagonals as the coefficients of the polynomial
                                                                                # Calculate its roots
        roots.sort(key = lambda x : abs(abs(x) - 1))                             
        roots_inside = [root for root in roots if ((abs(root) - 1) < 0)][:k]    # Take only roots which are inside unit circle
        
        roots_angels = np.angle(roots_inside)                                   # Calculate the phase component of the roots 
        DOA_pred = np.arcsin((1/(2 * np.pi * self.dist)) * roots_angels)        # Calculate the DOA out of the phase component
        DOA_pred = (180 / np.pi) * DOA_pred                                     # Convert from radians to degrees
        
        roots_angels_all = np.angle(roots)                                      # Calculate the phase component of the roots 
        DOA_pred_all = np.arcsin((1/(2 * np.pi * self.dist)) * roots_angels_all)                              # Calculate the DOA our of the phase component
        DOA_pred_all = (180 / np.pi) * DOA_pred_all                                     # Convert from radians to Deegres
        return DOA_pred, roots, k, DOA_pred_all, roots_angels_all

    #############################################################################

    def spectrum_calculation(self, Un, f=1, Array_form="ULA"):
        Spectrum_equation = []
        for angle in self.angels:
            a = self.system_model.SV_Creation(theta= angle, f= f, Array_form = Array_form)
            # print(a.shape)
            a = a[:Un.shape[0]]                                         # sub-array response for Spatial smoothing 
            Spectrum_equation.append(np.conj(a).T @ Un @ np.conj(Un).T @ a)
        Spectrum_equation = np.array(Spectrum_equation, dtype=np.complex)
        Spectrum = 1 / Spectrum_equation
        return Spectrum, Spectrum_equation
