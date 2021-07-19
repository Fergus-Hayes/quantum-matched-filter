import numpy as np
import matplotlib.pyplot as plt
import argparse, time, matplotlib, os
import quantum_matched_filter_functions as qmffn
import matplotlib.animation as animation
import scipy.stats as ss
import seaborn as sns
import pandas as pd

np.random.seed(int(time.time()))

def main(infiles, outpath, noisefile=False, bank='bank', fontsize=28, ticksize=22, figsize=(18,11), tempfile=None):

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
    
    snr_max = False
    
    for infile in infiles:
        filename = infile.split('/')[-1]
        snr_filename = 'snrs_'+'_'.join(filename.split('_')[2:])
        if os.path.isfile('/'.join(infile.split('/')[:-1])+'/'+snr_filename):
            snr_maxv = np.max(np.load('/'.join(infile.split('/')[:-1])+'/'+snr_filename))
            snr_maxind = np.argmax(np.load('/'.join(infile.split('/')[:-1])+'/'+snr_filename))
            snr_max = True
        SNRs.append(float(infile.split('/')[-1].split('_')[4]))
        psi_opts.append(np.load(infile))
    
    snr_inds = np.argsort(SNRs)
    SNRs = np.array(SNRs)[snr_inds]
    infiles = np.array(infiles)[snr_inds]
    psi_opts = np.array(psi_opts)[snr_inds]

    M = psi_opts.shape[1]
    temp_bank, M1, M2 = bankfunc(M, temp_file=tempfile)

    effspin = (temp_bank['mass1']*temp_bank['spin1z'] + temp_bank['mass2']*temp_bank['spin2z'])/(temp_bank['mass1']+temp_bank['mass2'])
    totalM = temp_bank['mass1']+temp_bank['mass2']
    chirpM = np.power(temp_bank['mass1']*temp_bank['mass2'],3./5)/np.power(temp_bank['mass1']+temp_bank['mass2'],1./5)

    snr_label = r'$\rho_{\regular{th}}$'
    effstr = r'$\chi_{\regular{eff}}$'
    m1_label = r'$m_{1}$ ($M_{\odot})$'
    m2_label = r'$m_{2}$ ($M_{\odot})$'
    Mt_label = r'$M_{\regular{tot}}$ ($M_{\odot})$'
    Mc_label = r'$\mathcal{M}$ ($M_{\odot})$'

    if snr_max:
        snr_dict = {}
        snr_dict[m1_label] = temp_bank['mass1'][snr_maxind]
        snr_dict[m2_label] = temp_bank['mass2'][snr_maxind]
        snr_dict[effstr] = effspin[snr_maxind]
        snr_dict[Mt_label] = totalM[snr_maxind]
        snr_dict[Mc_label] = chirpM[snr_maxind]

        print('Best SNR fit- m1:',snr_dict[m1_label], 'm2:', snr_dict[m2_label], 's1:', temp_bank['spin1z'][snr_maxind], 's2:', temp_bank['spin2z'][snr_maxind])

    points = pd.DataFrame()
    points[m1_label]=temp_bank['mass1']
    points[m2_label]=temp_bank['mass2']
    points[effstr]=effspin
    points[Mt_label]=totalM
    points[Mc_label]=chirpM

    cond_str ='> '+str(int(SNRs[0]))

    SNR_label = np.tile(cond_str,M)
    SNR_list = np.zeros(M)
    prob_label = np.tile(cond_str,M)

    for i in np.arange(psi_opts.shape[0]):
        probs = np.abs(psi_opts[i])**2
        matches = probs>np.mean(probs)
        snr_str = int(SNRs[i])
        SNR_label[matches] = snr_str
        prob_label[matches] = str(np.log(np.unique(probs[matches])[0]))
        SNR_list[matches] = SNRs[i]

    points[snr_label] = SNR_label
    points['lnp'] = prob_label

    cond = SNR_list!=0.

    points_m = pd.DataFrame()

    points_m[m1_label]=points[m1_label][cond]
    points_m[m2_label]=points[m2_label][cond]
    points_m[effstr]=points[effstr][cond]
    points_m[snr_label]=points[snr_label][cond]
    points_m[Mt_label]=points[Mt_label][cond]
    points_m[Mc_label]=points[Mc_label][cond]
    points_m['lnp']=points['lnp'][cond]

    points_m = points_m.sort_values(snr_label)

    size_max = 125.
    size_min = 15.

    size = size_max * ((points_m['lnp'].astype(float) - np.min(points_m['lnp'].astype(float))) / (np.max(points_m['lnp'].astype(float)) - np.min(points_m['lnp'].astype(float)))) + size_min

    para_strs = [m1_label, m2_label, effstr]#[Mc_label, Mt_label, effstr]#[m1_label,m2_label,Mt_label,effstr,Mc_label]

    ndim=len(para_strs)

    fig, axes = plt.subplots(nrows=ndim-1, ncols=ndim-1, figsize=(12,12))

    if ndim==2:
        axes = np.array([[axes]])

    count=0
    coords = np.array(np.array(np.tril_indices(ndim-1)).T)
    paras_inds = np.array(np.array(np.tril_indices(ndim,k=-1)).T)
    
    ticksize=16
    fontsize=22

    alpha2 = 0.3
    alpha = 1.
    cmap = plt.cm.jet

    for coord,para_inds in zip(coords,paras_inds):
        colors = iter(['#0173B2','#029E73','#DE8F05','#D55E00'])#iter(cmap(np.linspace(0,1,len(SNRs)+1)))
        #col = next(colors)
        count+=1
        scs = []
        plt_labels=[]
        if snr_max:
            label=None
            if count==len(coords):
                label=r'$\max(\{\rho^{\regular{max}}_{0},...,\rho^{\regular{max}}_{N}\})=$'+str(np.round(snr_maxv,2))
            sc_snr = axes[coord[0],coord[1]].scatter(snr_dict[para_strs[para_inds[1]]], snr_dict[para_strs[para_inds[0]]], color='black', marker='o', lw=8., alpha=alpha2, label=label)
        
        for SNR in SNRs:
            if count==len(coords):
                plt_labels.append(snr_label+'='+str(int(SNR)))
                label=snr_label+'='+str(int(SNR))
            else:
                label=None
            col = next(colors)
            snr_cond = np.array(points_m[snr_label])==str(int(SNR))
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
        if snr_max:
            if count==len(coords):
                label=r'$\max(\{\rho^{\regular{max}}_{0},...,\rho^{\regular{max}}_{N}\})=$'+str(np.round(snr_maxv,2))
                plt_labels.append(label)
                scs.append(sc_snr)

    for coord in np.array(np.array(np.triu_indices(ndim-1, k=1)).T):
        axes[coord[0],coord[1]].axis('off')
        
    leg = fig.legend(scs, plt_labels, loc='upper right', fontsize=fontsize)
    leg.get_frame().set_linewidth(0.0)
    leg.set_bbox_to_anchor((1, .9))
    #leg.set_bbox_to_anchor((.95, .35))#leg.set_bbox_to_anchor((.85, .85))
    

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
    parser.add_argument('--noisefile', help="", type=str, default=False)

    opt = parser.parse_args()
 
    main(opt.infile, opt.outpath, noisefile=opt.noisefile, bank=opt.bank, tempfile=opt.tempfile)
