import numpy as np
import tensorflow as tf
from math import e

###########################################
###### generate the dictionary matrix ######
###########################################


def ULA_narrow (n, m):
  A = torch.zeros(n, m, dtype=torch.cdouble)
  for i in range(n):
    for s in range(m):
      x = -1j * i * np.pi * np.sin(s * (np.pi)/m - np.pi/2)
      # x = -1j * i * np.pi * np.sin((s-90) * (np.pi)/m)
      # x = -1j * i * np.pi * np.sin(s * (np.pi)/m )
      # x_t = torch.tensor([x.real, x.imag])
      # y = x.real + x.imag * 1j
      # R[i,s] = (e ** x).real
      # I[i,s] = (e ** x).imag
      A[i,s] = e ** x
  return [A]
