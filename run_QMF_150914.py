import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from functools import partial

from pycbc.waveform import get_td_waveform
from pycbc.filter import matched_filter, make_frequency_series
from pycbc.types import TimeSeries, FrequencySeries

import time, os, h5py, multiprocessing

debug = False

# Set the random number generator seed 
np.random.seed(int(time.time()))

matplotlib.rc('text', usetex = True)
matplotlib.rc('font', **{'family' : "sans-serif"})
params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
figsize=(10,8)

# Set the threshold SNR
SNR = 12.
# Load the data (GW150914) and the corresponding psd
Data = np.load('signal.npy')
psd = np.load('psd.npy')

fig = plt.figure(figsize=figsize)
#plt.title('Figure 1')
Data_t = np.fft.irfft(Data)
plt.plot(np.arange(len(Data_t))/2048., Data_t/1.e-25, lw=.2, color='black')
plt.xlabel(r'$t$ (s)')
plt.ylabel(r'$x(t)$ ($10^{-25}$)')
plt.savefig('plots/strain.png')
plt.close()

# Qubits for templates
Mqubits = 18
# Qubits for Grover's
Pqubits = 7

# N is the number of data points
N = len(Data)
# Number of templates
M = int(2**(Mqubits))
# Number of states in ancillary register
P = int(2**(Pqubits))

# Generate indicies over time
n_inds = np.arange(N)
# Generate indicies over templates
m_inds = np.arange(M)
# Create indicies for each Grover's iteration
p_inds = np.arange(P)

def get_paras(M):
    '''
    Get mass/spins given index
    ''' 
    
    bank_file = 'H1L1-BANK2HDF-1134450017-1202400.hdf'#'H1L1-BANK2HDF_SPLITBANK_BANK0_INJECTIONS-1134450017-1202400.hdf'
    full_bank = h5py.File(bank_file,'r')

    bank_size = full_bank['mass1'].size

    if bank_size>M:
        indexes = np.arange(bank_size)[::int(bank_size//M)][:M]
    else:
        indexes = np.arange(bank_size)
    
    bank = {}
    bank['mass1'] = np.array(full_bank['mass1'])[indexes]   
    bank['mass2'] = np.array(full_bank['mass2'])[indexes]
    bank['spin1z'] = np.array(full_bank['spin1z'])[indexes]
    bank['spin2z'] = np.array(full_bank['spin2z'])[indexes]  
    #bank['f_lower'] = np.array(full_bank['f_lower'])[indexes]

    return bank

paras = get_paras(M)

def SNR_given_index(data, psd, dt, df, paras):
    '''
    Perform matched filtering on some data
    given the psd, delta f (df), delta t (dt)
    and parameters to make a waveform
    (mass1, mass2, spin1z, spin2z, f low).
    '''

    m1,m2,s1,s2 = paras
    N = len(data)
    waveform_ = get_td_waveform(approximant="IMRPhenomPv3",
                                mass1=m1, mass2=m2,
                                spin1z=s1, spin2z=s2,
                                delta_t=dt, f_lower=20.)
    waveform_[0].resize(int(2*(N - 1)))
    template_ = waveform_[0].cyclic_time_shift(waveform_[0].start_time)
    template = np.array(make_frequency_series(TimeSeries(template_,delta_t=dt)))
    SNR = matched_filter(FrequencySeries(template, delta_f=df),
                FrequencySeries(data, delta_f=df),
                FrequencySeries(psd, delta_f=df))
    SNR = np.array(SNR.crop(8, 4))
    return np.max(np.abs(SNR))

def k_12(M, Data, dt=1./2048, f_low=20., threshold=12):
    '''
    Make wavforms from index
    '''

    N = len(Data)
    bank = get_paras(M)

    psd = np.load('psd.npy')
    freqs = np.fft.fftfreq(2*(N-1))*1./dt
    df = np.abs(freqs[1]-freqs[0])

    m1s = bank['mass1']
    m2s = bank['mass2']
    s1s = bank['spin1z']
    s2s = bank['spin2z']
    #fls = bank['f_lower']

    #paras = iter(np.array([m1s,m2s,s1s,s2s,fls]).T)
    paras = iter(np.array([m1s,m2s,s1s,s2s]).T)

    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    func = partial(SNR_given_index, Data, psd, dt, df)
    SNRs = np.array(pool.map(func, paras))
    pool.close()
    w = np.where(SNRs>=threshold,-1.,1.)

    if len(w)<M:
        w = np.pad(w,(0,M-len(w)),constant_values=1.).astype(float)
    w*=1./np.sqrt(M)
    return w

# Define the psi ini state (currently not used)
psi_ini = None # Just used as a placeholder
# Apply k12 to the states
w = k_12(M, Data, threshold=SNR)

# Initialising joint state psi_0 which we apply Grover's algorithm over given w
psi_0 = np.ones((M,P))/np.sqrt(M*P)

def U_w(w):
    '''
    Defines the matrix operation applied to psi_0 given w.
    This operator flips the phase of a state that corresponds to w.
    '''
    M = len(w)
    return np.sqrt(M)*np.multiply(np.eye(M),w)

def U_s(w):
    '''
    Defines the matrix operation applied to psi_0 given w.
    This is "Grovers diffusion operator".
    '''
    M = len(w)
    s_i = np.array([np.ones(M)/np.sqrt(M)]).T
    return 2*np.outer(s_i,s_i)-np.eye(M)

def Grovers1(w, psi_0):
    '''
    This corresponds to the first part of Grovers algorithm.
    Here operators U_w and U_s are applied to psi_0 "p" times for each of the P states.
    '''
    M, P = psi_0.shape[0], psi_0.shape[1]
    psi_1 = np.zeros((M,P))
    for p in np.arange(P):
        psi_1[:,p] = np.dot(np.linalg.matrix_power(np.matmul(U_s(w),U_w(w)),p),psi_0[:,p])
    return psi_1

# Apply the first part of Grover's algorithm
psi_1 = Grovers1(w,psi_0)

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

X, Y = np.meshgrid(p_inds, m_inds)
Z=psi_1
fig = plt.figure(figsize=(14,10))
ax = fig.gca(projection='3d')
#ax.set_title('Figure 4')

surf = ax.plot_surface(X, Y, Z, cmap=cm.bwr,
                       linewidth=0, antialiased=False)

ax.set_xlabel(r'$P$')
#ax.set_ylabel(r'$M$')
ax.set_ylabel(r'$N$')
ax.set_zlabel(r'$|\psi_{1}\rangle$')

plt.tight_layout()
plt.savefig('plots/amplitudes.png')
plt.close()

def IQFT(P):
    '''
    Constructing the inverse fourier transform of size PxP.
    '''
    return np.round(
            np.array([[(np.exp(-2.*i*j*1j*np.pi/P))
                   for j in np.arange(P)]
                    for i in np.arange(P)]),1)*1./np.sqrt(P)

def Grovers2(psi_1):
    '''
    The seconds part of Grovers algorithm. 
    Here we apply the inverse quantum Fourier transform to psi_2.
    '''
    M, P = psi_1.shape[0], psi_1.shape[1]
    return np.dot(IQFT(P),psi_1.T)

# We apply the second part of Grover's algoirthm to state psi_1
psi_2 = Grovers2(psi_1)

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from matplotlib import cm

X, Y = np.meshgrid(p_inds, m_inds)
Z=np.absolute(psi_2).T**2
fig = plt.figure(figsize=(14,10))
ax = fig.gca(projection='3d')
#ax.set_title('Figure 5')

surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

ax.set_xlabel(r'$P$')
ax.set_ylabel(r'$M$')
ax.set_zlabel(r'$|\langle\psi_{2}\rangle|^{2}$')
plt.tight_layout()
plt.savefig('plots/IQFT.png')
plt.close()

# Measure each of the qubits corresponding to template states M and period states P
measurement = np.unravel_index(np.argmax(np.absolute(psi_2)**2), psi_2.shape)

# The resulting state measured from the P qubits tells us about the number of matching templates
P_measured = measurement[0]

# The result from measuring P corresponds to the frequency of the sinesoid shown in Figure 3
N_templates = int(np.round(M*np.sin(P_measured*np.pi/P)**2))

print('Number of matching templates:',N_templates)

# Calculate the optimum number of Grover's applications to find a
# matching template.
opt_p = int(np.round((np.pi/4.)* np.sqrt(M/N_templates)))

print('Optimum number of Grover applications:',opt_p)

psi_0_ = np.ones((M))/np.sqrt(M)
psi_1_opt = np.dot(np.linalg.matrix_power(np.matmul(U_s(w),U_w(w)),opt_p),psi_0_)

probs_out = []
for p in np.arange(4*opt_p+1):
    psi_1_ = np.dot(np.linalg.matrix_power(np.matmul(U_s(w),U_w(w)),p),psi_0_)
    probs_out.append(np.max(np.abs(psi_1_)**2))

probs_out = np.array(probs_out)

#plt.title('Figure 6')
fig = plt.figure(figsize=figsize)
plt.plot(probs_out, color='black')
plt.xlabel(r"\textrm{Grover's iterations}")
plt.ylabel(r'$max(|\langle\psi_{2}\rangle|^{2})$')#Probability of max probability state')
plt.axvline(opt_p, ls=':', color='black', label=r'$k_{\mathrm{t}}$')
#plt.axvline(opt_p, ls=':', color='black', label=r'$p_{opt}$')
plt.legend(fontsize=15)
plt.savefig('plots/popt.png')
plt.close()

fig = plt.figure(figsize=(10,7))
ax = plt.subplot(111)
#plt.title('Figure 8')
sc=ax.scatter(np.append(paras['mass1'], 0.), np.append(paras['mass2'], 40.), c=np.append(np.abs(psi_1_opt)**2, -.05), marker='o', lw=1, cmap=cm.Greys)
cb = plt.colorbar(sc, ax=[ax], label=r'$|\langle\psi_{2}\rangle|^{2}$')
plt.xlabel(r'$m_{1}$')
plt.ylabel(r'$m_{2}$')
plt.savefig('plots/mass_dists.png')
plt.close()

fig = plt.figure(figsize=(10,7))
ax = plt.subplot(111)
#plt.title('Figure 8')
sc=ax.scatter(np.append(paras['spin1z'], 1.), np.append(paras['spin2z'], 1.), c=np.append(np.abs(psi_1_opt)**2, -.05), marker='o', lw=1, cmap=cm.Greys)
cb = plt.colorbar(sc, ax=[ax], label=r'$|\langle\psi_{2}\rangle|^{2}$')
plt.xlabel(r'$s_{1}$')
plt.ylabel(r'$s_{2}$')
plt.savefig('plots/spin_dists.png')
plt.close()

Data = np.load('noise.npy')
# N is the number of data points
N = len(Data)

w = k_12(M, Data, threshold=SNR)

# Initialising joint state psi_0 which we apply Grover's algorithm over given w
psi_0_n = np.ones((M,P))/np.sqrt(M*P)
# Apply the first part of Grover's algorithm
psi_1_n = Grovers1(w,psi_0_n)
# We apply the second part of Grover's algoirthm to state psi_1
psi_2_n = Grovers2(psi_1_n)

# Measure each of the qubits corresponding to template states M and period states P
measurement_n = np.unravel_index(np.argmax(np.absolute(psi_2_n)**2), psi_2_n.shape)

# The resulting state measured from the P qubits tells us about the number of matching templates
P_measured_n = measurement_n[0]

# The result from measuring P corresponds to the frequency of the sinesoid shown in Figure 3
N_templates_n = int(np.round(M*np.sin(P_measured_n*np.pi/P)**2))

print('Number of matching templates given noise:',N_templates_n)
