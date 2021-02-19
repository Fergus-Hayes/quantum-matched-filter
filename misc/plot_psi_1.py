import numpy as np
import matplotlib.pyplot as plt
import argparse, time, matplotlib

np.random.seed(int(time.time()))

def main(infiles, outpath, fontsize=28, ticksize=22, figsize=(10,8)):

    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'

    SNRs = []
    psi_1s = []
    for infile in infiles:
        SNRs.append(float(infile.split('/')[-1].split('_')[4]))
        psi_1s.append(np.load(infile))
    
    snr_inds = np.argsort(SNRs)
    SNRs = np.array(SNRs)[snr_inds]
    infiles = np.array(infiles)[snr_inds]

    label_match = r'$\rho \geq $'
    label_nmatch = r'$\rho < $'
    ylabel = r'Probability amplitude'
    xlabel = r'$i$'
    cmap = plt.cm.jet

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(figsize[0],figsize[1]))

    colors = iter(cmap(np.linspace(0,1,len(SNRs)+1)))
    col = next(colors)

    for i,psi_1 in enumerate(psi_1s):

        col = next(colors)
        match_ind = np.argmax(np.sum(np.abs(psi_1)**2,axis=1))
        nmatch_ind = np.argmin(np.sum(np.abs(psi_1)**2,axis=1))

        psi_match = psi_1[match_ind]
        psi_nmatch = psi_1[nmatch_ind]

        ax.plot(psi_1[match_ind], lw=2, label=label_match+str(int(SNRs[i])), color=col)
        #axes[1].plot(psi_1[nmatch_ind], lw=2, label=label_nmatch+'='+str(int(SNRs[i])), color=col, ls='--')

    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    #axes[1].set_ylabel(ylabel, fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=ticksize)
    #axes[1].tick_params(axis='both', labelsize=ticksize)
    leg = ax.legend(fontsize=3*fontsize//4)#, loc='upper right')
    leg.get_frame().set_linewidth(0.0)
    fig.savefig(outpath+'_'.join(infiles[-1].split('/')[-1].split('_')[:4])+'_snr_'+'_'.join(SNRs.astype(int).astype(str))+'_psi1_all_match_paper.png')

    ylabel = 'Probability'

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=figsize)#nrows=len(SNRs)//2, ncols=2, figsize=(2*figsize[0],len(SNRs)*figsize[1]//2))
    
    axes = [axes]

    colors = iter(cmap(np.linspace(0,1,len(SNRs)+1)))
    col = next(colors)

    for i,psi_1 in enumerate([psi_1s[-1]]):

        col = 'black'#next(colors)
        match_ind = np.argmax(np.sum(np.abs(psi_1)**2,axis=1))
        nmatch_ind = np.argmin(np.sum(np.abs(psi_1)**2,axis=1))

        psi_match = psi_1[match_ind]
        psi_nmatch = psi_1[nmatch_ind]

        matches_ind = np.abs(psi_1[:,np.argmax(psi_match)])**2>np.mean(np.abs(psi_1[:,np.argmax(psi_match)])**2)
        nmatches_ind = np.abs(psi_1[:,np.argmax(psi_match)])**2<=np.mean(np.abs(psi_1[:,np.argmax(psi_match)])**2)

        psi_matches = psi_1[matches_ind]
        psi_nmatches = psi_1[nmatches_ind]

        axes[i].plot(np.sum(np.abs(psi_matches)**2, axis=0), lw=2, label=label_match+str(int(SNRs[-1])), color=col)
        axes[i].plot(np.sum(np.abs(psi_nmatches)**2, axis=0), lw=2, label=label_nmatch+str(int(SNRs[-1])), color=col, ls='--')
        axes[i].set_ylabel(ylabel, fontsize=fontsize)
        axes[i].set_xlabel(xlabel, fontsize=fontsize)
        axes[i].tick_params(axis='both', labelsize=ticksize)
        #leg = axes[i].legend(fontsize=3*fontsize//4)
        #leg.get_frame().set_linewidth(0.0)


    fig.tight_layout()
    #axes[-1].set_xlabel(xlabel, fontsize=fontsize)
    fig.savefig(outpath+'_'.join(infiles[-1].split('/')[-1].split('_')[:4])+'_snr_'+'_'.join(SNRs.astype(int).astype(str))+'_psi1_paper_v2.png')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Perform one QMF on GW150914 given a template bank")
    parser.add_argument('--infile', help="", type=str, nargs='+', required=True)
    parser.add_argument('--outpath', help="", type=str, required=True)

    opt = parser.parse_args()
 
    main(opt.infile, opt.outpath)
