import numpy as np
import matplotlib.pyplot as plt
import argparse, time, matplotlib, os

np.random.seed(int(time.time()))

def main(infiles, outpath, noisefile=False, fontsize=28, ticksize=22, figsize=(12,8)):

    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'

    SNRs = []
    for infile in infiles:
        SNRs.append(float(infile.split('/')[-1].split('_')[4]))
   
    snr_inds = np.argsort(SNRs)
    SNRs = np.array(SNRs)[snr_inds]
    infiles = np.array(infiles)[snr_inds]

    ylabel = r'$p(k)$'
    xlabel = r'$k$'
    ylabel2 = r'$p(t)$'
    xlabel2 = r'$t$'

    cmap = plt.cm.jet
    lw=3
    ms=6
    
    fig, (ax, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(figsize[0],2*figsize[1]))

    colors = iter(cmap(np.linspace(0,1,len(SNRs)+1)))
    col = next(colors)

    if noisefile:
        infiles = np.concatenate(([noisefile],infiles))
        SNRs = np.concatenate((['noise'],SNRs))

    for i,infile in enumerate(infiles):

        probs = np.sum(np.abs(np.load(infile))**2,axis=1)[1:]
        
        P,M = np.load(infile).shape
        bs = np.arange(P)[1:]
        k_bs_ = (P/(4.*bs))-0.5
        k_bs__ = (P/(4.*(P-bs)))-0.5

        k_bs = np.where(k_bs_>0.,np.floor(k_bs_),np.where(k_bs__,np.floor(k_bs__),1.))
        k_probs = np.zeros(len(np.unique(k_bs)))

        for j,k_b in enumerate(np.unique(k_bs)):
            k_probs[j] = np.sum(probs[k_b==k_bs])
        
        if i==0 and noisefile:
            label = r'$\rho_{\regular{th}}$=8 without signal'
            col = 'black'
            ls = '--'
        else:
            label=r'$\rho_{\regular{th}}$='+str(int(SNRs[i]))
            col = next(colors)
            ls = '-'

        ax.plot(np.unique(k_bs), k_probs, color=col, label=label, marker='o', lw=0.)

        delta = 0.0
        for k_b, k_prob in zip(np.unique(k_bs), k_probs):
            ax.axvline(k_b, ymin=0.+delta, ymax=k_prob+delta, color=col, ls=ls, lw=lw)

        theta_t = np.where(bs<P/2,bs*np.pi/P,(bs*np.pi/P)+np.pi)
        matches = np.floor(M*np.sin(theta_t)**2)
        match_probs = np.zeros(len(np.unique(matches)))

        for j,match in enumerate(np.unique(matches)):
            match_probs[j] = np.sum(probs[match==matches])

        ax2.plot(np.unique(matches), match_probs, color=col, marker='o', lw=0., ms=ms)

        delta = 0.0
        for match, match_prob in zip(np.unique(matches), match_probs):
            ax2.axvline(match, ymin=0.+delta, ymax=match_prob+delta, color=col, ls=ls, lw=lw)

        filename = infile.split('/')[-1]
        snr_filename = 'snrs_'+'_'.join(filename.split('_')[2:])

        if os.path.isfile('/'.join(infile.split('/')[:-1])+'/'+snr_filename):
            snrs = np.load('/'.join(infile.split('/')[:-1])+'/'+snr_filename)
            t_t = np.sum(snrs>=SNRs[i])
            k_t = np.round(((np.pi/(2.*np.arcsin(np.sqrt(t_t/M))))-1.)/2.)
            ax.axvline(k_t, ymin=0., ymax=1., color=col, ls=':')
            ax2.axvline(t_t, ymin=0., ymax=1., color=col, ls=':')


    ax.set_xscale('symlog')
    ax.set_xlim(-0.1)
    ax.set_ylim(0.,1.)

    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=ticksize, top=False, right=False)

    ax2.set_xscale('symlog')
    ax2.set_xlim(-0.1)
    ax2.set_ylim(0.,1.)
    ax2.set_xlabel(xlabel2, fontsize=fontsize)
    ax2.set_ylabel(ylabel2, fontsize=fontsize)
    ax2.tick_params(axis='both', labelsize=ticksize, top=False, right=False)

    leg = fig.legend(fontsize=3*fontsize//4, bbox_to_anchor=[0.3, 0.75])
    leg.get_frame().set_linewidth(0.0)
    fig.savefig(outpath+'_'.join(infiles[-1].split('/')[-1].split('_')[:4])+'_snr_'+'_'.join(SNRs.astype(int).astype(str))+'_psi2_kt_paper.png',bbox_inches='tight')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Perform one QMF on GW150914 given a template bank")
    parser.add_argument('--infile', help="", type=str, nargs='+', required=True)
    parser.add_argument('--outpath', help="", type=str, required=True)
    parser.add_argument('--noisefile', help="", type=str, default=False)

    opt = parser.parse_args()
 
    main(opt.infile, opt.outpath, noisefile=opt.noisefile)
