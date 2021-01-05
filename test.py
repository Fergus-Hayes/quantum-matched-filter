import gw_detections_functions as gwfunc
from pycbc.waveform import get_td_waveform
import matplotlib.pyplot as plt
import numpy as np

m1, m2 = 33., 33.
fs = 4096
T = 28
inv_psd = 1./np.load('data/psd.npy')

_, temp = gwfunc.make_template(m1,m2,fs,T,inv_psd)

pycbc_temp_, _ = get_td_waveform(approximant="IMRPhenomPv3",
                                mass1=m1, mass2=m2,
                                spin1z=0., spin2z=0.,
                                delta_t=1./fs, f_lower=20.)

pycbc_temp = np.zeros(int(fs*T))
pycbc_temp[:len(pycbc_temp_)] = pycbc_temp_
idx = np.argmax(np.abs(pycbc_temp))
pycbc_temp = np.roll(pycbc_temp,-idx+(int(T*fs)//2))

plt.plot(temp)
plt.plot(pycbc_temp)
plt.show()
