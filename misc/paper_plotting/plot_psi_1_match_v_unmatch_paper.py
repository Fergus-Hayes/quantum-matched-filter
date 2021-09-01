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
    infiles = np.array(infiles)[snr_inds]

    label_match = r'$P(\regular{Match})$'#r'$\rho \geq \rho_{\regular{th}} (= 18)$'
    label_nmatch = r'1-$P(\regular{Match})$'#r'$\rho < \rho_{\regular{th}} (= 18)$'
    ylabel = r'Probability'
    xlabel = r'Number of Grover operations'
    cmap = plt.cm.jet

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=figsize)#nrows=len(SNRs)//2, ncols=2, figsize=(2*figsize[0],len(SNRs)*figsize[1]//2))
    
    axes = [axes]

    colors = iter(cmap(np.linspace(0,1,len(SNRs)+1)))
    col = next(colors)

    for i,infile in enumerate([infiles[-1]]):

        psi_1 = np.load(infile)

        col = 'black'#next(colors)

        psi_match = psi_1[np.argmax(np.sum(np.abs(psi_1)**2,axis=1))]
        psi_nmatch = psi_1[np.argmin(np.sum(np.abs(psi_1)**2,axis=1))]

        psi_matches = psi_1[np.abs(psi_1[:,np.argmax(psi_match)])**2>np.mean(np.abs(psi_1[:,np.argmax(psi_match)])**2)]
        psi_nmatches = psi_1[np.abs(psi_1[:,np.argmax(psi_match)])**2<=np.mean(np.abs(psi_1[:,np.argmax(psi_match)])**2)]

        norm = 1./(np.sum(np.abs(psi_matches)**2, axis=0)+np.sum(np.abs(psi_nmatches)**2, axis=0))

        #print(np.argmax(np.sum(np.abs(psi_matches)**2, axis=0)))
        
        t = psi_matches.shape[0]
        N = psi_1.shape[0]
        theta = np.arcsin(np.sqrt(t/N))
        k_t = ((np.pi/(2.*theta))-1.)/2.

        axes[i].annotate(text='', xy=(k_t,1.005), xytext=(1,1.005), arrowprops=dict(arrowstyle='<->', lw=1.5))
        axes[i].annotate(text=r'$k$', xy=((k_t/2)-5,.95), fontsize=fontsize)

        axes[i].plot(norm*np.sum(np.abs(psi_matches)**2, axis=0), lw=2, label=label_match, color=col)
        axes[i].plot(norm*np.sum(np.abs(psi_nmatches)**2, axis=0), lw=2, label=label_nmatch, color=col, ls='--')
        axes[i].set_ylabel(ylabel, fontsize=fontsize)
        axes[i].set_xlabel(xlabel, fontsize=fontsize)
        axes[i].tick_params(axis='both', labelsize=ticksize)
        axes[i].set_xlim(0,275)
        leg = axes[i].legend(fontsize=2*fontsize//3, framealpha=0.)
        #leg.get_frame().set_linewidth(0.0)
        leg.set_bbox_to_anchor((.74, .84))

    #axes[-1].set_xlabel(xlabel, fontsize=fontsize)
    fig.savefig(outpath+'_'.join(infiles[-1].split('/')[-1].split('_')[:4])+'_snr_'+'_'.join(SNRs.astype(int).astype(str))+'_psi1_paper_v2.png',bbox_inches='tight')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Perform one QMF on GW150914 given a template bank")
    parser.add_argument('--infile', help="", type=str, nargs='+', required=True)
    parser.add_argument('--outpath', help="", type=str, required=True)
    parser.add_argument('--noisefile', help="", type=str, default=False)

    opt = parser.parse_args()
 
    main(opt.infile, opt.outpath, noisefile=opt.noisefile)
