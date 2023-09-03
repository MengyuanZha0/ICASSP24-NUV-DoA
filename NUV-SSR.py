import torch
import numpy as np
from tqdm import tqdm
from time import time
from scipy.signal import find_peaks, peak_prominences



class Complex_NUVEstimator_k:     # Estimator for k known
    def __init__(self, ss_model, least_squares=True, vectorized=True, paper_scale=False):
        self.ss_model = ss_model 
        self.least_squares = least_squares
        self.vectorized = vectorized

        self.loss_mse = torch.nn.MSELoss(reduction='sum')
        self.loss_l1 = torch.nn.L1Loss(reduction='sum')


    def predict(self, y, k, r, max_iterations, convergence_threshold):
        ss_model = self.ss_model

        A_H = torch.zeros(ss_model.m, ss_model.n, dtype = torch.cdouble).cpu()
        A_H = ss_model.A.conj().T


        delta = []                                           # delta records the difference of norm of q between each iteration it+1 and it

        ### 1. Initial Guess ###
        q = torch.zeros(max_iterations + 2, ss_model.m, dtype = torch.cdouble).cpu()
        q = 0.01 * torch.ones(max_iterations + 2, ss_model.m, dtype=torch.cdouble).cpu() 


        
        ### 2. EM Algorithm ###
        iterations = max_iterations

        for it in range(max_iterations):
          # 2a. Precision Matrix
          q[it] = q[it]

          torch.diag(torch.square(q[it]))

          W_inv = torch.zeros(n, n, dtype=torch.cdouble)
          W_inv =  A @ torch.diag(torch.square(q[it])) @ A_H 

          W_inv = W_inv + (r/l) * torch.eye(ss_model.n, dtype = torch.cdouble)          
          W = torch.zeros(n, n, dtype=torch.cdouble)
          W = torch.linalg.inv(W_inv)

          # 2b. Gaussian Posteriori Distribution
          mean = torch.zeros(m, dtype=torch.cdouble) 
          variance = torch.zeros(m, dtype=torch.cdouble)

          mean = abs(torch.diag(torch.square(q[it])) @ A_H @ W @ y )
          variance =  abs(torch.pow(q[it], 2) -  torch.diag(torch.pow(q[it], 4)) @ torch.diagonal(A_H @ W @ A) )
          q[it+1] = torch.sqrt(torch.square(mean) + variance)
            
          if torch.norm(q[it + 1] - q[it]) < convergence_threshold:   # stopping criteria
              q[-1] = q[it + 1]
              iterations = it
              break
          else:
              delta.append(torch.norm(q[it + 1] - q[it]))
            

          q[-1]=q[it+1]


        ### 3. MAP Estimator of the sparse signal###

        u = torch.zeros(l, ss_model.m, dtype = torch.cdouble)
        u = torch.diag(torch.square(q[-1])) @ A_H @ W @ y

        #### find peaks ####

        Spectrum = abs(u)
        x = np.linspace(0, m, m, endpoint=False)
        DOA_pred,_ = find_peaks(Spectrum)
        DOA_pred = list(DOA_pred)
        DOA_pred.sort(key = lambda x: Spectrum[x], reverse = True)
        DOA = []
        for i in DOA_pred[0:k]:
          DOA.append(i* (180 / m) - 90)
        DOA.sort()

        return [abs(q[-1]), abs(u), DOA, iterations, delta]


    def PMSE(self, pred, DOA):
      prmse_list = []
      for p in list(permutations(pred, len(pred))):
          p = np.array(p)
          DOA = np.array(DOA)
          prmse_val = mean_squared_error(pred, DOA)
          # error = (((p - DOA) * np.pi / 180) + np.pi / 2) % np.pi - np.pi / 2
          # prmse_val = (1 / len(p)) * (np.linalg.norm(error) ** 2)
          prmse_list.append(prmse_val)
      return np.min(prmse_list)


