import numpy as np
import matplotlib.pyplot as plt
import argparse, time, matplotlib
import quantum_matched_filter_functions as qmffn
import matplotlib.animation as animation
import scipy.stats as ss
import seaborn as sns
import pandas as pd

np.random.seed(int(time.time()))

def main(infiles, outpath, bank='bank', fontsize=28, ticksize=22, figsize=(18,11), tempfile=None):

    if bank=='bank':
        bankfunc = qmffn.get_paras
    elif bank=='grid':
        bankfunc = qmffn.get_mass_grid
        spins=False
    else:
        raise ValueError(bank+' is not an option. Try either "bank" or "grid".')
        exit()

    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'

    SNRs = []
    psi_opts = []
    for infile in infiles:
        SNRs.append(float(infile.split('/')[-1].split('_')[4]))
        psi_opts.append(np.load(infile))
    snr_inds = np.argsort(SNRs)
    SNRs = np.array(SNRs)[snr_inds]
    infiles = np.array(infiles)[snr_inds]
    psi_opts = np.array(psi_opts)[snr_inds]

    M = psi_opts.shape[1]
    temp_bank, M1, M2 = bankfunc(M, temp_file=tempfile)

    effspin = (temp_bank['mass1']*temp_bank['spin1z'] + temp_bank['mass2']*temp_bank['spin2z'])/(temp_bank['mass1']+temp_bank['mass2'])

    effstr = r'$\chi_{\regular{eff}}$'

    points = pd.DataFrame()
    points[r'$m_{1}$']=temp_bank['mass1']
    points[r'$m_{2}$']=temp_bank['mass2']
    points[effstr]=effspin

    cond_str ='> '+str(int(SNRs[0]))

    SNR_label = np.tile(cond_str,M)
    SNR_list = np.zeros(M)

    for i in np.arange(psi_opts.shape[0]):
        probs = np.abs(psi_opts[i])**2
        matches = probs>np.mean(probs)
        snr_str = int(SNRs[i])
        SNR_label[matches] = snr_str
        SNR_list[matches] = SNRs[i]

    points[r'$\rho_{th}$'] = SNR_label

    cond = SNR_list!=0.

    points_m = pd.DataFrame()

    points_m[r'$m_{1}$']=points[r'$m_{1}$'][cond]
    points_m[r'$m_{2}$']=points[r'$m_{2}$'][cond]
    points_m[effstr]=points[effstr][cond]
    points_m[r'$\rho_{th}$']=points[r'$\rho_{th}$'][cond]

    points_m = points_m.sort_values(r'$\rho_{th}$')

    size_max = 200.
    size_min = 10.

    size = size_max * ((points_m[r'$\rho_{th}$'].astype(float) -np.min(SNRs)) / np.max(SNRs)) + size_min

    para_strs = [r'$m_{1}$',r'$m_{2}$',effstr]

    ndim=len(para_strs)

    fig, axes = plt.subplots(nrows=ndim-1, ncols=ndim-1, figsize=(12,12))
  
    count=0
    coords = np.array(np.array(np.tril_indices(ndim-1)).T)
    paras_inds = np.array(np.array(np.tril_indices(ndim,k=-1)).T)
    
    ticksize=16
    fontsize=22

    alpha=0.7
    cmap = plt.cm.jet

    for coord,para_inds in zip(coords,paras_inds):
        colors = iter(cmap(np.linspace(0,1,len(SNRs)+1)))
        col = next(colors)
        count+=1
        scs = []
        plt_labels=[]
        for SNR in SNRs:
            if count==len(coords):
                plt_labels.append(r'$\rho_{th}=$'+str(int(SNR)))
                label=r'$\rho_{th}=$'+str(int(SNR))
            else:
                label=None
            col = next(colors)
            snr_cond = np.array(points_m[r'$\rho_{th}$'])==str(int(SNR))
            scs.append(axes[coord[0],coord[1]].scatter(points_m[para_strs[para_inds[1]]][snr_cond], points_m[para_strs[para_inds[0]]][snr_cond], s=size[snr_cond], color=col, marker='.', alpha=alpha, label=label))
        axes[coord[0],coord[1]].spines['right'].set_visible(False)
        axes[coord[0],coord[1]].spines['top'].set_visible(False)
        axes[coord[0],coord[1]].tick_params(axis='both', labelsize=ticksize)
        if coord[0]==ndim-2:
            axes[coord[0],coord[1]].set_xlabel(para_strs[para_inds[1]], fontsize=fontsize)
        else:
            axes[coord[0],coord[1]].set_xticks([])
        if coord[1]==0:
            axes[coord[0],coord[1]].set_ylabel(para_strs[para_inds[0]], fontsize=fontsize)
        else:
            axes[coord[0],coord[1]].set_yticks([])

    for coord in np.array(np.array(np.triu_indices(ndim-1, k=1)).T):
        axes[coord[0],coord[1]].axis('off')
        
    leg = fig.legend(scs, plt_labels, loc='upper right', fontsize=fontsize)
    leg.get_frame().set_linewidth(0.0)
    leg.set_bbox_to_anchor((.85, .85))

    for handle in leg.legendHandles:
        handle.set_sizes([size_max])

    fig.tight_layout()

    fig.savefig(outpath+'_'.join(infiles[-1].split('/')[-1].split('_')[:4])+'_snr_'+'_'.join(SNRs.astype(int).astype(str))+'_thrs_paper.png',bbox_inches='tight') 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Perform one QMF on GW150914 given a template bank")
    parser.add_argument('--infile', help="", type=str, nargs='+', required=True)
    parser.add_argument('--outpath', help="", type=str, required=True)
    parser.add_argument('--bank', help="", type=str, default='bank')
    parser.add_argument('--tempfile', help="", type=str, default=None)

    opt = parser.parse_args()
 
    main(opt.infile, opt.outpath, bank=opt.bank, tempfile=opt.tempfile)
