import numpy as np
import torch
from math import e

class DataGenerator:
    def __init__(self, k, gap, sample, n, x_var, mean_c, r2, l):
        self.k = k                                                              # number of sources
        self.gap = gap                                                          # minimum gap between each source directions
        self.sample = sample                                                  
        self.n = n                                                              # number of sensor arrays
        self.x_var = x_var                                                      # variance of non-coherent Gaussian source signal
        self.mean_c = mean_c                                                    # constant mean of non-coherent Gaussian source signal
        self.r2 = r2                                                            # variance of complex AWGN
        self.l = l                                                              # number of snapshots

    def doa_generator(self):
        """Generates k DOAs with a specified gap."""
        while True:
            doa = np.round(np.random.rand(self.k) * 150, decimals=2) - 75
            doa.sort()
            difference_between_angles = np.diff(doa)

            if np.all(difference_between_angles > self.gap) and np.all(difference_between_angles < (180 - self.gap)):
                return doa

    def steeringmatrix(self, doa):
        """Generates samples of steering matrices given DOA."""
        STM_A = torch.zeros(self.sample, self.n, self.k, dtype=torch.cdouble)

        for j in range(self.sample):
            sig_direc = np.deg2rad(doa[j])

            for i in range(self.n):
                for s in range(self.k):
                    x = -1j * i * np.pi * np.sin(sig_direc[s])
                    STM_A[j][i][s] = e ** x

        return STM_A

    def nonco_signal_generator(self):
        """Generates non-coherent source signals."""
        x = np.sqrt(self.x_var / 2) * (np.random.randn(1, self.k) + 1j * np.random.randn(1, self.k))
        return torch.tensor(x)


def generate_experiment_data(generator):
    """Experiment data generation using the DataGenerator methods."""
    x_dire = [generator.doa_generator() for _ in range(generator.sample)]
    
    STM_A = generator.steeringmatrix(x_dire)

    x_mean = generator.mean_c * torch.ones(generator.k)
    x_true = torch.zeros(generator.sample, l, generator.k, 1, dtype=torch.cdouble)
    
    for j in range(generator.sample):
        for t in range(l):
            x_true[j, t, :, 0] = generator.nonco_signal_generator() + x_mean

    y_train = torch.zeros(generator.sample, l, generator.n, 1, dtype=torch.cdouble)
    
    for j in range(generator.sample):
        for t in range(l):
            er1 = torch.tensor(np.random.normal(0, generator.r2 / 2, generator.n))
            er2 = torch.tensor(np.random.normal(0, generator.r2 / 2, generator.n))
            y_train[j, t, :, 0] = STM_A[j].matmul(x_true[j, t, :, 0]) + er1 + er2 * 1j
            
    return x_dire, x_true, y_train


# Initialization
k = 1
gap = 15
sample = 1
n = 16
x_var = 0.5
mean_c = 2
r2 = 1e-1
l = 100

# Generate data
generator = DataGenerator(k, gap, sample, n, x_var, mean_c, r2, l)
x_dire, x_true, y_train = generate_experiment_data(generator)
# print(x_dire)
