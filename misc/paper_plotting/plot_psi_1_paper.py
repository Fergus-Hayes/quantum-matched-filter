import numpy as np
import matplotlib.pyplot as plt
import argparse, time, matplotlib

np.random.seed(int(time.time()))

def main(infiles, outpath, noisefile=False, fontsize=28, ticksize=22, figsize=(12,8)):

    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'

    SNRs = []
    for infile in infiles:
        SNRs.append(float(infile.split('/')[-1].split('_')[4]))
    
    snr_inds = np.argsort(SNRs)
    SNRs = np.array(SNRs)[snr_inds]

    label_match = r'$\rho_{\regular{thr}}=$' #r'$\rho \geq $'
    label_nmatch = r'$\rho < $'
    ylabel = r'Probability amplitude'
    xlabel = r'Counting register state index $j$'
    cmap = plt.cm.jet

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(figsize[0],figsize[1]))

    colors = iter(cmap(np.linspace(0,1,len(SNRs)+1)))
    col = next(colors)

    for i,infile in enumerate(infiles):

        psi_1 = np.load(infile)
        P = psi_1.shape[1]
        psi_1 = psi_1[np.argmax(np.sum(np.abs(psi_1)**2,axis=1))]

        col = next(colors)

        ax.plot(psi_1, lw=2, label=label_match+str(int(SNRs[i])), color=col)
        del psi_1

    ax.set_xlim(0,P)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    #axes[1].set_ylabel(ylabel, fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=ticksize)
    #axes[1].tick_params(axis='both', labelsize=ticksize)
    leg = ax.legend(fontsize=2*fontsize//3, framealpha=1., loc='upper right')
    #leg.get_frame().set_linewidth(0.0)
    fig.savefig(outpath+'_'.join(infiles[-1].split('/')[-1].split('_')[:4])+'_snr_'+'_'.join(SNRs.astype(int).astype(str))+'_psi1_all_match_paper.png', bbox_inches='tight')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Perform one QMF on GW150914 given a template bank")
    parser.add_argument('--infile', help="", type=str, nargs='+', required=True)
    parser.add_argument('--outpath', help="", type=str, required=True)
    parser.add_argument('--noisefile', help="", type=str, default=False)

    opt = parser.parse_args()
 
    main(opt.infile, opt.outpath, noisefile=opt.noisefile)
