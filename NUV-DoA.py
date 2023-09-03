import torch
import numpy as np
from tqdm import tqdm
from time import time
from scipy.signal import find_peaks, peak_prominences
from scipy import signal

# Try Hermitian conjugate here
# A_H = R_T - I_T * 1j

class Complex_NUVEstimator_k:
    def __init__(self, ss_model, least_squares=True, vectorized=True, paper_scale=False):
        self.ss_model = ss_model  # groundtruth of the system y = Phi @ x + v, have some properties, e.g. dimension = m, mean = mu...
        self.least_squares = least_squares
        self.vectorized = vectorized
        # self.paper_scale = paper_scale
        # self.square_window = square_window

        self.loss_mse = torch.nn.MSELoss(reduction='sum')
        self.loss_l1 = torch.nn.L1Loss(reduction='sum')

    def predict(self, y, resol, k, r, A, A_H, m, max_iterations, convergence_threshold):
        ss_model = self.ss_model
        centers = np.linspace(B1, B2, int((B2 - B1) / resol + 1), endpoint=True)
        center_size = len(centers)  # center_size = # of windows used

        U = torch.zeros(center_size, m, dtype=torch.float64).to('cuda')
        u = torch.zeros(center_size, dtype=torch.float64).to('cuda')  # len(centers) = # of windows
        middle_index = int(m / 2)
        delta = []  # delta records the difference of norm of q between each iteration it+1 and it

        # 1. Initial Guess
        q = 0.1 * torch.ones(center_size, 2, m, dtype=torch.cdouble).to('cuda')

        # 2. EM Algorithm
        iterations = max_iterations

        for it in range(max_iterations):
            q[:, 0, :] = q[:, -1, :]
            W_inv = torch.zeros(center_size, n, n, dtype=torch.cdouble).to('cuda')
            diag_squared_q = torch.diag_embed(torch.square(q[:, 0, :]), offset=0, dim1=-2, dim2=-1).to('cuda')
            W_inv = A @ diag_squared_q @ A_H

            id = torch.eye(ss_model.n, dtype=torch.cdouble).to('cuda')
            id = id.reshape((1, ss_model.n, ss_model.n))
            id_batch = id.repeat(center_size, 1, 1)

            W_inv = W_inv + (r / l) * id_batch

            W = torch.zeros(center_size, n, n, dtype=torch.cdouble).to('cuda')
            W = torch.linalg.inv(W_inv)

            mean = torch.zeros(center_size, m, dtype=torch.cdouble).to('cuda')
            variance = torch.zeros(center_size, m, dtype=torch.cdouble).to('cuda')

            mean = abs(diag_squared_q @ A_H @ W @ y)

            M1 = torch.diag_embed(torch.pow(q[:, 0, :], 4), offset=0, dim1=-2, dim2=-1).to('cuda')
            M2 = torch.diagonal(A_H @ W @ A, dim1=-2, dim2=-1).to('cuda')
            M2 = M2[:, :, None]

            variance = abs(torch.pow(q[:, 0, :], 2) - torch.squeeze(torch.bmm(M1, M2)))

            q[:, -1, :] = torch.sqrt(torch.square(mean) + variance)

            if max(torch.linalg.norm(q[:, 0, :] - q[:, -1, :], dim=1)) < convergence_threshold:
                q[:, -1, :] = q[:, -1, :]
                iterations = it
                break
            else:
                delta.append(max(torch.linalg.norm(q[:, 0, :] - q[:, -1, :], dim=1)))

        # 3. MAP Estimator
        U = abs(torch.diag_embed(torch.square(q[:, -1, :]), offset=0, dim1=-2, dim2=-1) @ A_H @ W @ y).to('cuda')
        u = U[:, middle_index]

        Spectrum = u.cpu()
        Spectrum_all = U.cpu()
        m2 = len(Spectrum)

        x = np.linspace(0, m2, m2, endpoint=False)
        DOA_pred, _ = find_peaks(Spectrum)
        DOA_pred = list(DOA_pred)
        DOA_pred.sort(key=lambda x: u[x], reverse=True)
        DOA = []
        spec_pred = []
        for i in DOA_pred[0:k]:
            DOA.append(centers[i])
            spec_pred.append(Spectrum_all[i, :])
        DOA.sort()

        return [u, DOA, spec_pred, iterations, delta]

    def PMSE(self, pred, DOA):
        prmse_list = []
        for p in list(permutations(pred, len(pred))):
            p = np.array(p)
            DOA = np.array(DOA)
            prmse_val = mean_squared_error(pred, DOA)
            prmse_list.append(prmse_val)
        return np.min(prmse_list)
