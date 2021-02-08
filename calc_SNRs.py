import numpy as np
import gw_detections_functions as gwfn
import quantum_matched_filter_functions as qmfn
import h5py, multiprocessing, os
from functools import partial

M = 2**19
snrfunc = qmfn.snr_given_index

Data = np.load('data/noise.npy')
psd = np.load('data/noise_psd.npy')

dt = 1./4096

bank, _, _ = qmfn.get_paras(M)

M = len(bank['mass1'])

paras = iter(np.array(list(bank.values())).T)

cores = 3

pool = multiprocessing.Pool(cores)
func = partial(snrfunc, Data, psd, dt)
SNRs = np.array(pool.map(func, paras))
pool.close()

np.save('data/SNRs_noise_spins', SNRs)
