import numpy as np
import gw_detections_functions as gwfn

from pycbc.waveform import get_td_waveform

import h5py, multiprocessing
from functools import partial

def get_paras(M, temp_file='data/template_bank.hdf'):
    '''
    Get mass/spins given index
    ''' 
    
    full_bank = h5py.File(temp_file,'r')

    bank_size = full_bank['mass1'].size
    
    if bank_size>M:
        indexes = np.arange(bank_size)[::int(bank_size//M)][:M]
    else:
        indexes = np.arange(bank_size)
        
    bank = {}
    bank['mass1'] = np.array(full_bank['mass1'],dtype=float)[indexes]   
    bank['mass2'] = np.array(full_bank['mass2'],dtype=float)[indexes]
    bank['spin1z'] = np.array(full_bank['spin1z'],dtype=float)[indexes]
    bank['spin2z'] = np.array(full_bank['spin2z'],dtype=float)[indexes] 
    bank['f_lower'] = np.array(full_bank['f_lower'],dtype=float)[indexes]

    return bank, None, None

def get_mass_grid(M, mmin=4., mmax=90., temp_file=None):
    '''
    Get mass/spins given index
    ''' 
    
    Nms = int(np.sqrt(2.*M))
    ms = np.linspace(mmin, mmax, Nms)
    M1s, M2s = np.meshgrid(ms,ms)
    m1s = M1s[M1s>=M2s]
    m2s = M2s[M1s>=M2s]
    bank = {}
    bank['mass1'] = m1s
    bank['mass2'] = m2s
    
    return bank, M1s, M2s

def snr_given_index(data, psd, dt, df, paras):
    '''
    Perform matched filtering on some data
    given the psd, delta f (df), delta t (dt)
    and parameters to make a waveform 
    (mass1, mass2, spin1z, spin2z, f low).
    '''
    
    fs = int(1./dt)
    
    m1,m2,s1,s2,flow = paras
    N = len(data)
    T = int(N/fs)
    _, temp = gwfn.make_template(m1, m2, fs, T, 1./psd, s1=s1, s2=s2, f_low=flow)
    if np.any(temp):
        SNR = gwfn.get_snr(data, temp, fs)
        SNR = SNR[8*fs:-4*fs]
        return np.max(np.abs(SNR))
    else:
        return 0

def k_12(index_states, Data, psd, dt=1./4096, f_low=20., threshold=12, spins=False, bankfunc=get_paras, temp_file=None):
    '''
    Make wavforms from index
    '''

    M = index_states.shape[0]
    N = len(Data)

    freqs = np.fft.fftfreq(2*(N-1))*1./dt
    df = np.abs(freqs[1]-freqs[0])

    bank, _, _ = get_paras(M, temp_file=temp_file)
    m1s = bank['mass1']
    m2s = bank['mass2']
    s1s = np.zeros(len(m1s))
    s2s = np.zeros(len(m1s))
    flows = bank['f_lower']
    if spins:
        s1s = bank['spin1z']
        s2s = bank['spin2z']

    paras = iter(np.array([m1s,m2s,s1s,s2s,flows]).T)

    #SNRs = []
    #for para in paras:
    #    print(para)
    #    SNRs.append(snr_given_index(Data, psd, dt, df, para))
    #    print(SNRs[-1])
    #SNRs = np.array(SNRs)

    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    func = partial(snr_given_index, Data, psd, dt, df)
    SNRs = np.array(pool.map(func, paras))
    pool.close()

    w = np.where(SNRs>=threshold, -1., 1.)

    M = len(w)
    w*=1./np.sqrt(M)
    return w, SNRs

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

def quantum_counting0(w, psi_0, dtype='float64'):
    M = len(w)
    U_psi_0 = psi_0.astype(dtype)*w*np.sqrt(M)
    s = np.ones(M).astype(dtype)/np.sqrt(M)
    return 2*np.dot(s,U_psi_0).astype(dtype)*s - U_psi_0

def iquantum_counting0(w, psi_0, itts, dtype='float64'):
    M = len(w)
    psi_ = psi_0.astype(dtype)
    for itt in np.arange(itts):
        psi_ = quantum_counting0(w, psi_, dtype)
    return psi_

def quantum_counting1(w, psi_0, dtype='float64'):
    '''
    This corresponds to the first part of Grovers algorithm.
    Here operators U_w and U_s are applied to psi_0 "p" times for each of the P states.
    '''
    M, P = psi_0.shape
    psi_1 = np.zeros((M,P)).astype(dtype)
    psi_1[:,0] = np.ones(M).astype(dtype)/np.sqrt(M*P)
    for p in np.arange(1,P):
        psi_1[:,p] = quantum_counting0(w, psi_1[:,p-1], dtype)
    return psi_1

def IQFT(P):
    '''
    Constructing the inverse fourier transform of size PxP.
    '''
    return np.round(
            np.array([[(np.exp(-2.*i*j*1j*np.pi/P))
                   for j in np.arange(P)]
                    for i in np.arange(P)]),1)*1./np.sqrt(P)

def quantum_counting2(psi_1):
    '''
    The seconds part of Grovers algorithm. 
    Here we apply the inverse quantum Fourier transform to psi_2.
    '''
    M, P = psi_1.shape[0], psi_1.shape[1]
    return np.dot(IQFT(P),psi_1.T)

def QMF(Data, psd, M, P, tag='out', path='./', SNR_threshold=12., bankfunc=get_paras, table=False, save_states=False, dtype='float64', temp_file='data/template_bank.hdf'):
    '''
    Runs the full algorithm given the input data, number of template qubits M,
    number of precision qubits P, output path and threshold SNR (default 12).
    Returns the template states after optimal applcations are applied.
    '''
    
    N = len(Data)

    # Create equal superposition across template states
    index_states = np.ones(M).astype(dtype)/np.sqrt(M)
    
    # Apply k1 and k2 to get w states
    w, _ = k_12(index_states, Data, psd, threshold=SNR_threshold, bankfunc=bankfunc, temp_file=temp_file)
    M = len(w)

    # Apply the first step to quantum counting
    psi_ = quantum_counting1(w,np.ones((M,P)).astype(dtype)/np.sqrt(M*P), dtype=dtype)
    
    # Apply the first step to quantum counting
    psi = quantum_counting2(psi_)
    
    # Save states
    if save_states:
        np.save(path+'/psi1_in_'+tag,psi_)
        np.save(path+'/psi2_in_'+tag,psi)

    # Measure the most probable state
    measurement = np.unravel_index(np.argmax(np.absolute(psi)**2), psi.shape)
    
    # From the measured state, we figure out the number of matches
    N_templates = int(np.round(M*np.sin(measurement[0]*np.pi/P)**2))
    
    if N_templates==0.:
        return measurement[0], np.ones((M))/np.sqrt(M), np.sum(w<0.)

    # The estimated optimal number of Grover's applications can then be determined
    opt_p = int(np.round(((np.pi/4)/np.arcsin(np.sqrt(N_templates/M))) - 1./2))
    
    # As we know the correct number of matches, we can calculate the true optimal applications
    opt_t = int(np.round(((np.pi/4)/np.arcsin(np.sqrt(int(np.sum(w<0.))/M))) - 1./2))
    
    # The state of the template states after the optimal applications is then determined
    psi_opt = iquantum_counting0(w, np.ones(M).astype(dtype)/np.sqrt(M), opt_p, dtype)
    
    # The number of matching template can then be determined
    p_correct = np.sum(np.array(np.abs(psi_opt)**2)[np.abs(psi_opt)**2>np.mean(np.abs(psi_opt)**2)])
    
    if table:
        print(SNR_threshold, '&', np.log2(P), '&', opt_p, '&', opt_t, '&', N_templates, '&',int(np.sum(w<0.)),'&',2**(np.log2(P)-1)+opt_p,'&',np.round(p_correct,2),r'$\\ \\$')
    
    return measurement[0], psi_opt, np.sum(w<0.)