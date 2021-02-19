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

    ylabel = r'$P(b_{obs}=b)$'
    xlabel = r'$b$'
    cmap = plt.cm.jet
    
    fig, (ax, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(figsize[0],figsize[1]))

    colors = iter(cmap(np.linspace(0,1,len(SNRs)+1)))
    col = next(colors)

    if noisefile:
        infiles = np.concatenate(([noisefile],infiles))
        SNRs = np.concatenate(([0],SNRs))

    for i,infile in enumerate(infiles):

        prob = np.sum(np.abs(np.load(infile))**2,axis=1)
        P,M = np.load(infile).shape

        if i==0 and noisefile:
            label = r'$\rho_{\regular{th}}$=8 without signal'
            col = 'black'
            ls = '--'

        else:
            label=r'$\rho_{\regular{th}}$='+str(int(SNRs[i]))
            col = next(colors)
            ls = '-'

        ax.bar(np.arange(len(prob)),prob, color=col)
        ax2.bar(np.arange(len(prob)),prob, label=label, color=col)

        filename = infile.split('/')[-1]
        snr_filename = 'snrs_'+'_'.join(filename.split('_')[2:])

        if os.path.isfile('/'.join(infile.split('/')[:-1])+'/'+snr_filename):
            snrs = np.load('/'.join(infile.split('/')[:-1])+'/'+snr_filename)
            t_t = np.sum(snrs>=SNRs[i])
            theta = np.arcsin(np.sqrt(t_t/M))
            k_t = ((np.pi/(2.*theta))-1.)/2.
            b_t1 = np.round(P*theta/np.pi)
            b_t2 = np.round(P*(np.pi-theta)/np.pi)
            ax.axvline(b_t1, ymin=0., ymax=1., color=col, ls=':')
            ax2.axvline(b_t2, ymin=0., ymax=1., color=col, ls=':')


    stretch = 60 #len(prob)//9

    ax.set_xlim(0., stretch)  # outliers only
    ax2.set_xlim(len(prob)-stretch, len(prob))  # most of the data
    
    # hide the spines between ax and ax2
    ax.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax.yaxis.tick_left()
    ax.tick_params(labelright=False)
    ax2.yaxis.tick_right()
    
    d = .015 # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((1-d,1+d), (-d,+d), **kwargs)
    ax.plot((1-d,1+d),(1-d,1+d), **kwargs)

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d,+d), (1-d,1+d), **kwargs)
    ax2.plot((-d,+d), (-d,+d), **kwargs)

    #ax.set_xlabel(xlabel, fontsize=fontsize)
    #fig.suptitle(xlabel, y=-0.1, fontsize=fontsize)
    fig.text(.5, .02, xlabel, ha='center', fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    #axes[1].set_ylabel(ylabel, fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=ticksize)
    ax2.tick_params(axis='both', labelsize=ticksize, right=False, labelright=False)

    #axes[1].tick_params(axis='both', labelsize=ticksize)
    leg = fig.legend(fontsize=3*fontsize//4, bbox_to_anchor=[0.6, 0.8])#, loc='upper right')
    leg.get_frame().set_linewidth(0.0)
    fig.savefig(outpath+'_'.join(infiles[-1].split('/')[-1].split('_')[:4])+'_snr_'+'_'.join(SNRs.astype(int).astype(str))+'_psi2_b_paper.png',bbox_inches='tight')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Perform one QMF on GW150914 given a template bank")
    parser.add_argument('--infile', help="", type=str, nargs='+', required=True)
    parser.add_argument('--outpath', help="", type=str, required=True)
    parser.add_argument('--noisefile', help="", type=str, default=False)

    opt = parser.parse_args()
 
    main(opt.infile, opt.outpath, noisefile=opt.noisefile)
