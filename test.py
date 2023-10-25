##########################################################
################### simulating NUV-SSR ###################
##########################################################

#### generate dictionary matrix of size nxm ####

import numpy as np
import tensorflow as tf
from math import e

def ULA_narrow (n, m):
  A = torch.zeros(n, m, dtype=torch.cdouble).cpu()
  for i in range(n):
    for s in range(m):
      x = -1j * i * np.pi * np.sin(s * (np.pi)/m - np.pi/2)
      A[i,s] = e ** x
  return [A]
  
#### initialization####

m=360                                                         #number of grids      
[A] = ULA_narrow(n,m)                           #steering matrix given the true direction x_dire

ss_model = Complex_SystemModel(A, "Narrowband", k)
estimator = Complex_NUVEstimator_k(ss_model)

y_mean = torch.zeros(samples, ss_model.n, 1, dtype = torch.cdouble) # generate y_mean by averaging l snapshots for each sample
for i in range(samples):
  for j in range(ss_model.n):
    for s in range(l):
      y_mean[i, j, 0] += y_train[i, s, j, 0] / l

#### estimation ####

from sklearn.metrics import mean_squared_error
from itertools import permutations
samples_run = samples


r_t = [1e-0, 9e-1, 7e-1, 5e-1, 3e-1]                                  # searching the best performed tuning parameter

for r_tuning in r_t:
  print('======================================')
  print('r_tuning = {}'.format(r_tuning))
  x_pred = [0] * samples_run
  q_pred = [0] * samples_run
  theta = [0] * samples
  MSE = [0] * samples_run 
  for i in range(samples_run):
    [q_pred[i], x_pred[i], theta[i], iterations, deltas] = estimator.predict(y_mean[i, :, 0], k, r=r_tuning, max_iterations=10000,  convergence_threshold=4e-4)

    MSE[i] = estimator.PMSE([theta[i]], [x_dire[i]]) # mean square error for each sample
    
    # print('MSE for {}th sample = {}'.format(i, MSE[i]))
    # print('------------------------------------------')
  MSE_dB = 10 * (math.log10(np.mean(MSE)))
  print('averaged MSE in dB = {}'.format(MSE_dB))
  MSE_linear = math.sqrt(np.mean(MSE))
  print('averaged RMSE in linear = {}'.format(MSE_linear))
  print('--------------------------------------------')

SNR = 10*math.log10((x_var + mean_c) / r_2)
print('SNR = {}'.format(SNR))

#### plotting ####
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
# fig,axs = plt.subplots(179)


p = 0

fig,ax = plt.subplots(figsize=(5, 3))

x_c = x_pred[p].cpu()
y_c = np.linspace(-90, 90, len(x_c), endpoint=False)
  # y_c = np.linspace(-0.5, 0.5, len(x_c), endpoint=False)
ax.plot(y_c, x_c, marker = '*')

# plt.plot(0, x_c[50], "o")
ax.set_title('M = {}'.format(len(x_c) ))
plt.savefig('plot/Vanilla_m=360_snr_high.eps', format='eps')



##########################################################
################### simulating NUV-DoA ###################
##########################################################

#### generate dictionary matrices for each little windows ####
def ULA_refine (self, n, resol, m, b1, b2):                                               #[b1, b2] are the grid boundary for each individual windows
        angles = np.linspace( b1 * np.pi/180, b2 * np.pi/180, m, endpoint=True)
        A = torch.zeros(n, m, dtype=torch.cdouble)
   
        for i in range(n):
          for s in range(m):
             x = -1j * i * np.pi * np.sin(angles[s])
             A[i,s] = e ** x
        return [A]

#### initialization ####

ss_model = Complex_SystemModel("Narrowband", n, k)

estimator = Complex_NUVEstimator_k(ss_model)

print(y_train.shape)
y_mean = torch.zeros(samples, n, 1, dtype = torch.cdouble).to('cuda')
for i in range(samples):
  for j in range(n):
    for s in range(l):
      y_mean[i, j, 0] += y_train[i, s, j, 0] / l

from sklearn.metrics import mean_squared_error
from itertools import permutations



B1 = -89                                                                        # B1 and B2 are left and right boundary centers of all windows
B2 = 89
resol = 0.05                                                                    # resolution
centers = np.linspace( B1, B2, int((B2-B1)/resol+1), endpoint=True)             # all the center azimuths of all the windows
# print(centers)
center_size = len(centers)                                                      # number of centers
window_size = resol * 50                                                        # each window covers the interval [center - window_size, center + window_size]
# window_size = 0.5
print(window_size)
m = int((window_size * 2)/resol)+1                                              # grid size of each result vector of its center angle = number of steering vectors of each window, here m=101
print('m = {}'.format(m))
u = torch.zeros(len(centers), m, dtype = torch.float64)
A = torch.zeros(center_size, n, m, dtype = torch.cdouble).cpu()                 # dictionary matrices
A_H = torch.zeros(center_size, m, n, dtype = torch.cdouble).cpu()               # transpose of dictionary matrices
        # print('middle_index = {}'.format(middle_index))
for ct in range(len(centers)):
    print(centers[ct])
    [A[ct, :, :]] = ss_model.ULA_refine(n, resol, m, b1 = centers[ct] - window_size, b2 = centers[ct] + window_size)
    A_H[ct, :, :] = A[ct, :, :].conj().T


#### estimation ####
import gc
A = A.to('cuda')
A_H = A_H.to('cuda')
samples_run = samples
x_pred = [0] * samples_run
q_pred = [0] * samples
  # est = [0] * samples
MSE = [0] * samples_run
theta = [0] * samples_run
spec_pred = [0] * samples_run
spec_true = [0] * samples_run
for i in range(samples_run):                                                    #r=r_2 * 6000 when k=1
  print(i)
  print('true doa is {}'.format(x_dire[i]))
  r_tuning =  [20000, 15000, 10000, 300, 100, 90]
  MSE_tuning = [0] * len(r_tuning)
  for r in range(len(r_tuning)):
    print('r_tuning = {}'.format(r_tuning[r]))
    [x_pred[i], theta_tuning, _, iterations, deltas] = estimator.predict(y_mean[i, :, 0], resol = 0.05, k = k, r=r_tuning[r], A = A, A_H = A_H, m = m, max_iterations=3000,  convergence_threshold=1e-3)
    print('iterations = {}'.format(iterations))
    print('predicted doa is {}'.format(theta_tuning))
    MSE_tuning[r] = estimator.PMSE(theta_tuning, x_dire[i])
    print('MSE_tuning = {}'.format(MSE_tuning[r]))
    gc.collect()
    torch.cuda.empty_cache()

  MSE[i] = min(MSE_tuning)
  print('MSE[i] = {}'.format(MSE[i]))
  print('-----------------------------------------')
  gc.collect()
  torch.cuda.empty_cache()

MSE_dB = 10 * (math.log10(np.mean(MSE)))
print('averaged MSE in dB = {}'.format(MSE_dB))
MSE_linear = math.sqrt(np.mean(MSE))
print('averaged RMSE in linear = {}'.format(MSE_linear))
SNR = 10*math.log10((0.5 + 4) / r2)
print('SNR = {}'.format(SNR))

#### plotting ####

from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

p = 0

fig,ax = plt.subplots(figsize=(6, 3))

  # x_c = xx_pred[p]
  # axs = axs.ravel()
x_c = x_pred[p].cpu()
y_c = np.linspace(B1, B2, len(x_c), endpoint=False)
  # y_c = np.linspace(-0.5, 0.5, len(x_c), endpoint=False)
ax.plot(y_c, x_c, marker = '*')
# plt.plot(0, x_c[50], "o")
# ax.set_title('true doa = {}'.format(x_dire[p] ))
  # argrelextrema(np.array(x_c), np.greater)
print('true doa = {}'.format(x_dire[p]))
print('predicted doa = {}'.format(theta[p]))
print('σ = 3r2 = 30')



###############################################################################
############# Hierarchical Approach for Multi-source Estimation ###############
###############################################################################

k = 2

ss_model = Complex_SystemModel("Narrowband", n, k)

estimator = Complex_NUVEstimator_k(ss_model)
mm = 18000
# baseline = Model_Based_methods(ss_model)
vanilla = Vanilla_NUVEstimator(ss_model)           # we choose NUV-SSR as our baseline to cancel the interference from the observation 
array = np.linspace(0, n, n, endpoint=False)

print(y_train.shape)
y_mean = torch.zeros(samples, n, 1, dtype = torch.cdouble).to('cuda')
for i in range(samples):
  for j in range(n):
    for s in range(l):
      y_mean[i, j, 0] += y_train[i, s, j, 0] / l


from sklearn.metrics import mean_squared_error
from itertools import permutations
import gc

#### prediction based on NUV-SSR (vanilla) methods ####

doa_vanilla = [0] * samples
doa_pred = [[0] * k] * samples
resol = 0.05
interf_const = 2.0
MSE_vanilla_all = [0] * samples
MSE_NUV = [0] * samples

for i in range(samples):
  print(i)
  [_, _, doa_vanilla[i], iterations, deltas] = vanilla.predict(y_mean[i, :, 0], k, r=9, A = A_vanilla, max_iterations=9000,  convergence_threshold=4e-4)
  print('NUV vanilla prediction = {}'.format(doa_vanilla[i]))
  print('true direction = {}'.format(x_dire[i]))
  RMSE_vanilla = [0] * k
  B1 = [0] * k
  B2 = [0] * k
  for r in range(k):
    RMSE_vanilla[r] = math.sqrt(estimator.PMSE([doa_vanilla[i][r]], [x_dire[i][r]]))
    B1[r] = int(doa_vanilla[i][r] - 6 * RMSE_vanilla[r]) - 2
    B2[r] = int(doa_vanilla[i][r] + 6 * RMSE_vanilla[r]) + 2

  MSE_vanilla_all[i] = estimator.PMSE([doa_vanilla[i]], [x_dire[i]])
  print('MSE of vanilla NUV prediction = {}'.format(MSE_vanilla_all[i]))   
  print('B1 = {}'.format(B1))
  print('B2 = {}'.format(B2))

  for r in range(k):
    interf_vec = ss_model.ULA_intervec(n, 0.01, 1, angles = doa_vanilla[i][1-r]).to('cuda')
    y = y_mean[i, :, 0] - interf_const * interf_vec
    # print(y)
    [_, doa_nuv, _, iterations, deltas] = estimator.predict(y, resol = 0.05, k = 1, r=40, B1 = B1[r], B2 = B2[r], max_iterations=5000,  convergence_threshold=6e-4)
    doa_pred[i][r] = doa_nuv
    gc.collect()
    torch.cuda.empty_cache()
  print('NUV-EM prediction = {}'.format(doa_pred[i]))
  MSE_NUV[i] = estimator.PMSE([doa_pred[i]], [x_dire[i]])
  print(MSE_NUV[i])
  print('---------------------------------------------------')

MSE_NUV_dB = 10 * (math.log10(np.mean(MSE_NUV)))
print('NUV-predicted averaged MSE in dB = {}'.format(MSE_NUV_dB))
MSE_NUV_linear = math.sqrt(np.mean(MSE_NUV))
print('NUV-predicted averaged RMSE in linear = {}'.format(MSE_NUV_linear))
MSE_vanilla_dB = 10 * (math.log10(np.mean(MSE_vanilla_all)))
print('vanilla NUV averaged MSE in dB = {}'.format(MSE_vanilla_dB))
MSE_vanilla_linear = math.sqrt(np.mean(MSE_vanilla_all))
print('vanilla NUV averaged RMSE in linear = {}'.format(MSE_vanilla_linear))
