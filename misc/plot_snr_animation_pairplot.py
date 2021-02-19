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

    #sns.pairplot(points, hue=r'$\rho_{th}$', markers='.')
    #plt.savefig(outpath+'_'.join(infiles[-1].split('/')[-1].split('_')[:4])+'_snr_'+'_'.join(SNRs.astype(int).astype(str))+'_thrs_pairplot.png',bbox_inches='tight')

    points_m = pd.DataFrame()


    points_m[r'$m_{1}$']=points[r'$m_{1}$'][cond]
    points_m[r'$m_{2}$']=points[r'$m_{2}$'][cond]
    points_m[effstr]=points[effstr][cond]
    points_m[r'$\rho_{th}$']=points[r'$\rho_{th}$'][cond]

    points_m = points_m.sort_values(r'$\rho_{th}$')

    #points_m[r'$m_{1}$']=points_m[r'$m_{1}$'].iloc[sort_inds]
    #points_m[r'$m_{2}$']=points_m[r'$m_{2}$'].iloc[sort_inds]
    #points_m[r'$\chi_{eff}$']=points_m[r'$\chi_{eff}$'].iloc[sort_inds]
    #points_m[r'$\rho_{th}$']=points_m[r'$\rho_{th}$'].iloc[sort_inds]

    print(points_m)

    sns.set_context("paper", rc={"font.size":20,"axes.titlesize":20,"axes.labelsize":20})

    size_max = 75.
    size_min = 10.

    size = size_max * ((points_m[r'$\rho_{th}$'].astype(float) -np.min(SNRs)) / np.max(SNRs)) + size_min

    g = sns.pairplot(points_m, hue=r'$\rho_{th}$', markers='.', corner=True, hue_order=SNRs.astype(int).astype(str), plot_kws={'s': size, 'alpha': 0.7}, diag_kind='hist')
    
    plt.setp(g._legend.get_texts(), fontsize='19') # for legend text
    plt.setp(g._legend.get_title(), fontsize='19')
    g._legend.set_bbox_to_anchor((.8, .8))

    #g.map_diag(plt.hist, density=True)

    g.fig.savefig(outpath+'_'.join(infiles[-1].split('/')[-1].split('_')[:4])+'_snr_'+'_'.join(SNRs.astype(int).astype(str))+'_thrs_pairplot_matches_only.png',bbox_inches='tight') 

    exit()

    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=2, metadata=dict(artist='Me'), bitrate=1800)
    im_ani = animation.ArtistAnimation(fig, ims, interval=50, repeat_delay=300, blit=True)
    im_ani.save(outpath+'_'.join(infiles[-1].split('/')[-1].split('_')[:4])+'_snr_'+'_'.join(SNRs.astype(int).astype(str))+'_thrs.mp4', writer=writer)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Perform one QMF on GW150914 given a template bank")
    parser.add_argument('--infile', help="", type=str, nargs='+', required=True)
    parser.add_argument('--outpath', help="", type=str, required=True)
    parser.add_argument('--bank', help="", type=str, default='bank')
    parser.add_argument('--tempfile', help="", type=str, default=None)

    opt = parser.parse_args()
 
    main(opt.infile, opt.outpath, bank=opt.bank, tempfile=opt.tempfile)
