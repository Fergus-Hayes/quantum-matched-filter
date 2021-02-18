import numpy as np
import gw_detections_functions as gwfn

import h5py, multiprocessing, os
from functools import partial

def get_paras(M, temp_file='data/template_bank.hdf', spins=True, flow=20.):
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
    if spins:
        bank['spin1z'] = np.array(full_bank['spin1z'],dtype=float)[indexes]
        bank['spin2z'] = np.array(full_bank['spin2z'],dtype=float)[indexes] 
    else:
        bank['spin1z'] = np.zeros(M)
        bank['spin2z'] = np.zeros(M)
        
    if not flow:
        bank['f_lower'] = np.array(full_bank['f_lower'],dtype=float)[indexes]
    else:
        if not isinstance(flow,float):
            flow=20.
        bank['f_lower'] = np.ones(M)*flow

    return bank, None, None

def get_mass_grid(M, mmin=4., mmax=200., temp_file=None, spins=False):
    '''
    Get mass/spins given index
    ''' 
    
    Nms = int(np.sqrt(2.*M))
    ms = np.linspace(mmin, mmax, Nms)
    M1s, M2s = np.meshgrid(ms,ms)
    m1s = M1s[M1s>=M2s]
    m2s = M2s[M1s>=M2s]
    bank = {}
    bank['mass1'] = m1s[:M]
    bank['mass2'] = m2s[:M]
    bank['spin1z'] = np.zeros(M)
    bank['spin2z'] = np.zeros(M)
    bank['f_lower'] = np.ones(M)*20.
    
    return bank, M1s[:M,:M], M2s[:M,:M]

def snr_given_index(data, psd, dt, paras):
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

def k_12(index_states, Data, psd, dt=1./4096, f_low=20., threshold=12, spins=True, bankfunc=get_paras, temp_file=None, cores=1, snrfunc=snr_given_index):
    '''
    Make wavforms from index
    '''

    M = index_states.shape[0]
    N = len(Data)

    freqs = np.fft.fftfreq(2*(N-1))*1./dt
    df = np.abs(freqs[1]-freqs[0])

    bank, _, _ = bankfunc(M, temp_file=temp_file, spins=spins)
    
    paras = iter(np.array(list(bank.values())).T)
    
    if cores == None:
        cores = multiprocessing.cpu_count()

    pool = multiprocessing.Pool(cores)
    func = partial(snrfunc, Data, psd, dt)
    SNRs = np.array(pool.map(func, paras))
    pool.close()

    w = np.where(SNRs>=threshold, -1., 1.)
    w = np.pad(w,(0,M-len(w)),'constant',constant_values=(0,1))
    
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

def QMF(Data, psd, M, P, tag='out', path='./', SNR_threshold=12., bankfunc=get_paras, table=False, save_states=False, load_states=False, dtype='float64', temp_file='data/template_bank.hdf', spins=True, cores=1):
    '''
    Runs the full algorithm given the input data, number of template qubits M,
    number of precision qubits P, output path and threshold SNR (default 12).
    Returns the template states after optimal applcations are applied.
    '''
    
    N = len(Data)

    # Create equal superposition across template states
    index_states = np.ones(M).astype(dtype)/np.sqrt(M)
    
    if load_states and os.path.isfile(path+'/snrs_'+tag+'.npy'):
        print('Loading SNR')
        snrs = np.load(path+'/snrs_'+tag+'.npy')
        w = np.where(snrs>=SNR_threshold,-1.,1.)/np.sqrt(M)
    elif load_states and os.path.isfile('data/SNRs.npy'):
        print('Loading SNR (from data)')
        snrs_ = np.load('data/SNRs.npy')
        if len(snrs_)>M:
            snrs = snrs_[::len(snrs_)//M]
            if len(snrs)>M:
                snrs = snrs[:M]
            elif len(snrs)<M:
                snrs__ = snrs_[1:][::len(snrs)//M][:M-len(snrs)]
                snrs = np.concatenate((snrs,snrs__))
        snrs = snrs_
        w = np.where(snrs>=SNR_threshold,-1.,1.)/np.sqrt(M)


    else:
        print('Calculating SNR')
        # Apply k1 and k2 to get w states
        w, snrs = k_12(index_states, Data, psd, threshold=SNR_threshold, bankfunc=bankfunc, temp_file=temp_file, spins=spins, cores=cores)
    
    M = len(w)

    print(w)

    # Apply the first step to quantum counting
    psi_ = quantum_counting1(w,np.ones((M,P)).astype(dtype)/np.sqrt(M*P), dtype=dtype)
    
    # Apply the second step to quantum counting
    psi = quantum_counting2(psi_)
    
    # Save states
    if save_states:
        np.save(path+'/snrs_'+tag,snrs)
        np.save(path+'/psi1_in_'+tag,psi_)
        np.save(path+'/psi2_in_'+tag,psi)

    # Measure the most probable state
    measurement = np.unravel_index(np.argmax(np.absolute(psi)**2), psi.shape)
    
    # From the measured state, we figure out the number of matches
    N_templates = int(np.round(M*np.sin(measurement[0]*np.pi/P)**2))
    
    if N_templates==0.:
        return measurement[0], np.ones((M))/np.sqrt(M), np.sum(w<0.)

    # The estimated optimal number of Grover's applications can then be determined
    b = measurement[0]
    if b!=0 or b!=P-1:
        opt_p = int(np.round((P/(4*b))-0.5))
    else:
        opt_p = 1.
    
    # As we know the correct number of matches, we can calculate the true optimal applications
    opt_t = int(np.round(((np.pi/4)/np.arcsin(np.sqrt(int(np.sum(w<0.))/M))) - 1./2))
    
    # The state of the template states after the optimal applications is then determined
    psi_opt = iquantum_counting0(w, np.ones(M).astype(dtype)/np.sqrt(M), opt_p, dtype)
    
    # The probability of returning a matching template can then be determined
    p_correct = np.sum(np.array(np.abs(psi_opt)**2)[np.abs(psi_opt)**2>np.mean(np.abs(psi_opt)**2)])
    
    if table:
        print(SNR_threshold, '&', np.log2(P), '&', opt_p, '&', opt_t, '&', N_templates, '&',int(np.sum(w<0.)),'&',2**(np.log2(P)-1)+opt_p,'&',np.round(p_correct,2),r'$\\ \\$')
    
    return b, psi_opt, np.sum(w<0.)
