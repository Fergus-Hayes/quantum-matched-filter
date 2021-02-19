import numpy as np
import gw_detections_functions as gwfn
import quantum_matched_filter_functions as qmfn
import h5py, multiprocessing, os, argparse
from functools import partial

parser = argparse.ArgumentParser(usage='', description="Calculate SNRs of template bank")
parser.add_argument('--datafile', help="", type=str, default='data/signal.npy')
parser.add_argument('--psdfile', help="", type=str, default='data/psd.npy')
parser.add_argument('--loc', help="", type=str, default='./')

opt = parser.parse_args()

M = 2**19
snrfunc = qmfn.snr_given_index

Data = np.load(opt.datafile)
psd = np.load(opt.psdfile)

dt = 1./4096

bank, _, _ = qmfn.get_paras(M)

M = len(bank['mass1'])

cores = 3

batch = 1000
i = 0

while i < M:
    if os.path.isfile(opt.loc+'SNRs_noise_spins.npy'):
        SNRs_ = np.load(opt.loc+'SNRs_noise_spins.npy')
        i = len(SNRs_)
        j = i + batch
    else:
        SNRs_ = []
        i = 0
        j = batch
    if j>M:
        j = M

    print(i,j)

    paras_ = np.array(list(bank.values())).T
    paras = iter(paras_[i:j])

    pool = multiprocessing.Pool(cores)
    func = partial(snrfunc, Data, psd, dt)
    SNRs = np.array(pool.map(func, paras))
    SNRs = np.concatenate((SNRs,SNRs_))
    pool.close()

    np.save(opt.loc+'SNRs_noise_spins', SNRs)
