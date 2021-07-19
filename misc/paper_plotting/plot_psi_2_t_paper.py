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

    ylabel = r'$P(t_{\regular{obs}}=t)$'
    xlabel = r'$t$'
    cmap = plt.cm.jet
    lw=3
    ms=6
    
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(figsize[0],figsize[1]))

    colors = ['#bae4bc','#7bccc4','#43a2ca','#0868ac']#['orange','violet','blue','green']#iter(cmap(np.linspace(0,1,len(SNRs)+1)))
    #col = next(colors)

    if noisefile:
        infiles = np.concatenate(([noisefile],infiles))
        SNRs = np.concatenate((['noise'],SNRs))

    for i,infile in enumerate(infiles):

        probs = np.sum(np.abs(np.load(infile))**2,axis=1)
        P,M = np.load(infile).shape
        bs = np.arange(P)[1:]
        matches = np.round(M*np.sin(bs*np.pi/P)**2).astype(int)
        matches = np.append(0,np.where(matches==0,1,matches))
        match_probs = np.zeros(len(np.unique(matches)))

        for j,match in enumerate(np.unique(matches)):
            match_probs[j] = np.sum(probs[match==matches])

        if i==0 and noisefile:
            label = r'$\rho_{\regular{thr}}$=8 without signal'
            col = 'black'
            ls = '--'
        else:
            label=r'$\rho_{\regular{thr}}$='+str(int(SNRs[i]))
            col = colors[i]#next(colors)
            ls = '-'

        ax.plot(np.unique(matches), match_probs, color=col, marker='o', lw=0., ms=ms, label=label)
        
        delta = 0.0
        for match, match_prob in zip(np.unique(matches), match_probs):
            ax.axvline(match, ymin=0.+delta, ymax=match_prob+delta, color=col, ls=ls, lw=lw)

        filename = infile.split('/')[-1]
        snr_filename = 'snrs_'+'_'.join(filename.split('_')[2:])

        if os.path.isfile('/'.join(infile.split('/')[:-1])+'/'+snr_filename):
            snrs = np.load('/'.join(infile.split('/')[:-1])+'/'+snr_filename)
            t_t = np.sum(snrs>=SNRs[i])
            ax.axvline(t_t, ymin=0., ymax=1., color=col, ls=':')

    ax.set_xscale('symlog')
    ax.set_xlim(-0.1)
    ax.set_ylim(0.,1.)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=ticksize, top=False, right=False)

    leg = fig.legend(fontsize=3*fontsize//4, bbox_to_anchor=[0.3, 0.85])
    leg.get_frame().set_linewidth(0.0)
    fig.savefig(outpath+'_'.join(infiles[-1].split('/')[-1].split('_')[:4])+'_snr_'+'_'.join(SNRs.astype(int).astype(str))+'_psi2_matches_paper.png',bbox_inches='tight')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Perform one QMF on GW150914 given a template bank")
    parser.add_argument('--infile', help="", type=str, nargs='+', required=True)
    parser.add_argument('--outpath', help="", type=str, required=True)
    parser.add_argument('--noisefile', help="", type=str, default=False)

    opt = parser.parse_args()
 
    main(opt.infile, opt.outpath, noisefile=opt.noisefile)
