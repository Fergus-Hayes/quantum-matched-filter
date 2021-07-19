import numpy as np
import matplotlib.pyplot as plt
import argparse, time, matplotlib, weakref
import quantum_matched_filter_functions as qmffn
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

np.random.seed(int(time.time()))

def IQFT(P):
    '''
    Constructing the inverse fourier transform of size PxP.
    '''
    return np.round(
            np.array([[(np.exp(-2.*i*j*1j*np.pi/P))
                   for j in np.arange(P)]
                    for i in np.arange(P)]),1)*1./np.sqrt(P)

def main(N, P, t, outpath, fontsize=28, ticksize=22, figsize=(12,8)):

    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'
    
    fig = plt.figure(figsize=figsize)
    ax = [fig.gca()]
    
    lw=1.5
    fontsize=16
    
    theta = np.arcsin(np.sqrt(float(t)/N))
    state = [np.array([np.sqrt(float(t)/N),np.sqrt((float(N)-t)/N)])]

    ax[0].plot(np.arange(P), np.sin(((2.*np.arange(P)) + 1.)*theta), color='blue', label=r'$\rho\geq\rho_{\regular{th}}$ for '+str(t)+' matches')
    ax[0].plot(np.arange(P), np.sin(((2.*np.arange(P)) + 1.)*np.arcsin(np.sqrt(float(6.*t)/N))), color='red', label=r'$\rho\geq\rho_{\regular{th}}$ for '+str(6*t)+' matches')

    ax[0].set_xlabel(r'$k$', fontsize=fontsize)
    ax[0].set_ylabel('Amplitude', fontsize=fontsize)
    ax[0].tick_params(axis='both', labelsize=ticksize)

    theta = np.arcsin(np.sqrt(t/N))
    k_t = ((np.pi/(2.*theta))-1.)/2.

    ax[0].annotate(text='', xy=(k_t,1.005), xytext=(k_t,0.), arrowprops=dict(arrowstyle='-', ls=':'))
    ax[0].annotate(text=str(np.round(k_t).astype(int)), xy=(k_t,-0.05), fontsize=fontsize)

    k_t2 = ((np.pi/(2.*np.arcsin(np.sqrt(float(6.*t)/N))))-1.)/2.

    ax[0].annotate(text='', xy=(k_t2,1.005), xytext=(k_t2,0.), arrowprops=dict(arrowstyle='-', ls=':'))
    ax[0].annotate(text=str(np.round(k_t2).astype(int)), xy=(k_t2,-0.05), fontsize=fontsize)

    #ax[0].set_xticks(np.arange(P).astype(int))

    leg = ax[0].legend(fontsize=fontsize, loc='upper right')
    fig.savefig(outpath+'diff_t.png', bbox_inches='tight')
    plt.close()

    b_amp1 = np.matmul(IQFT(P),np.sin(((2.*np.arange(P)) + 1.)*theta))
    b_amp2 = np.matmul(IQFT(P),np.sin(((2.*np.arange(P)) + 1.)*np.arcsin(np.sqrt(float(6.*t)/N))))

    b_prob1 = np.abs(b_amp1)**2
    b_prob2 = np.abs(b_amp2)**2
    bs = np.arange(P)

    bs = bs[1:]
    b_prob1 = b_prob1[1:]/np.sum(b_prob1)
    b_prob2 = b_prob2[1:]/np.sum(b_prob2)

    k_bs_ = (P/(4.*bs))-0.5
    k_bs__ = (P/(4.*(P-bs)))-0.5

    k_bs = np.where(k_bs_>0.,np.round(k_bs_),np.where(k_bs__,np.round(k_bs__),1.))
    k_probs1 = np.zeros(len(np.unique(k_bs)))
    k_probs2 = np.zeros(len(np.unique(k_bs)))

    for j,k_b in enumerate(np.unique(k_bs)):
        k_probs1[j] = np.sum(b_prob1[k_b==k_bs])
        k_probs2[j] = np.sum(b_prob2[k_b==k_bs])

    fig = plt.figure(figsize=figsize)
    ax = fig.gca()   

    #ax.set_xscale('symlog')
    #ax.set_xlim(-0.1)
    ax.set_ylim(0.,1.)
    ax.set_xlabel(r'$k$', fontsize=fontsize)
    ax.set_ylabel(r'$P(k_{obs}=k)$', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=ticksize)

    ax.plot(np.unique(k_bs), k_probs1, color='blue', label=r'$\rho\geq\rho_{\regular{th}}$ for '+str(t)+' matches', marker='o', lw=0.)
    ax.plot(np.unique(k_bs), k_probs2, color='red', label=r'$\rho\geq\rho_{\regular{th}}$ for '+str(int(6.*t))+' matches', marker='o', lw=0.)

    for k_b, k_prob in zip(np.unique(k_bs), k_probs2):
        ax.axvline(k_b, ymin=0., ymax=k_prob, color='red', lw=2)

    for k_b, k_prob in zip(np.unique(k_bs), k_probs1):
        ax.axvline(k_b, ymin=0., ymax=k_prob, color='blue', lw=2)

    ax.legend(fontsize=fontsize, loc='upper right')
    fig.savefig(outpath+'diff_k.png', bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Perform one QMF on GW150914 given a template bank")
    parser.add_argument('--N', help="", type=int, required=True)
    parser.add_argument('--P', help="", type=int, required=True)
    parser.add_argument('--t', help="", type=int, required=True)
    parser.add_argument('--outpath', help="", type=str, required=True)

    opt = parser.parse_args()
 
    main(opt.N, opt.P, opt.t, opt.outpath)
