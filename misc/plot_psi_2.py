import numpy as np
import matplotlib.pyplot as plt
import argparse, time, matplotlib

np.random.seed(int(time.time()))

def main(infile, outpath, fontsize=28, ticksize=22, figsize=(10,8)):

    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'

    psi_2 = np.load(infile)

    fig = plt.figure(figsize=figsize)
    ax = fig.gca()

    #ax.scatter(np.arange(psi_2.shape[0]),np.sum(np.absolute(psi_2).T**2,axis=0), color='black', marker='.', lw=2)
    M = psi_2.shape[1]
    P = psi_2.shape[0]
    bs = np.arange(psi_2.shape[0])
    probs = np.sum(np.absolute(psi_2).T**2,axis=0)
    ax.bar(np.arange(psi_2.shape[0]), np.sum(np.absolute(psi_2).T**2,axis=0), color='black')
    
    theta_t = np.where(bs<P/2,bs*np.pi/P,(bs*np.pi/P)+np.pi)
    matches = np.floor(M*np.sin(theta_t)**2)

    match_probs = np.zeros(len(np.unique(matches)))

    for i,match in enumerate(np.unique(matches)):
        match_probs[i] = np.sum(probs[match==matches])

    measurement = np.unravel_index(np.argmax(np.absolute(psi_2)**2), psi_2.shape)
    ax.set_xlabel(r'$b$', fontsize=fontsize)
    ax.set_ylabel(r'$p(b)$', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=ticksize)
    fig.savefig(outpath+'.'.join(infile.split('/')[-1].split('.')[:-1])+'_psi2.png')
    
    fig2 = plt.figure(figsize=figsize)
    ax2 = fig2.gca()
    
    #ax2.bar(np.unique(matches), match_probs, color='black', width=100)
    ax2.plot(np.unique(matches), match_probs, color='black', marker='o', lw=0., ms=5)

    ax2.set_xlabel(r'Number of matching templates from counting', fontsize=fontsize)
    ax2.set_ylabel(r'Probability', fontsize=fontsize)
    ax2.tick_params(axis='both', labelsize=ticksize)
    ax2.set_xscale('log')
    fig.tight_layout()
    fig2.savefig(outpath+'.'.join(infile.split('/')[-1].split('.')[:-1])+'_psi2_matches.png')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Perform one QMF on GW150914 given a template bank")
    parser.add_argument('--infile', help="", type=str, required=True)
    parser.add_argument('--outpath', help="", type=str, required=True)

    opt = parser.parse_args()
 
    main(opt.infile, opt.outpath)
