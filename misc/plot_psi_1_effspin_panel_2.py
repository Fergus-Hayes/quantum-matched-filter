import numpy as np
import matplotlib.pyplot as plt
import argparse, time, matplotlib
import quantum_matched_filter_functions as qmffn
import matplotlib.animation as animation
import scipy.stats as ss

np.random.seed(int(time.time()))

def main(infiles, outpath, bank='bank', fontsize=28, ticksize=22, figsize=(24,8), tempfile=None):

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

    print(SNRs)
    print(M)

    effspin = (temp_bank['mass1']*temp_bank['spin1z'] + temp_bank['mass2']*temp_bank['spin2z'])/(temp_bank['mass1']+temp_bank['mass2'])

    fig = plt.figure(figsize=figsize)
    gs = matplotlib.gridspec.GridSpec(1, 3, figure=fig)

    axs = []
    axs.append(fig.add_subplot(gs[0,0]))
    axs.append(fig.add_subplot(gs[0,1]))
    axs.append(fig.add_subplot(gs[0,2]))

    axs[0].set_xlim(np.min(temp_bank['mass1']),np.max(temp_bank['mass1']))
    axs[0].set_ylim(np.min(temp_bank['mass2']),np.max(temp_bank['mass2']))
    axs[1].set_xlim(np.min(temp_bank['mass1']),np.max(temp_bank['mass1']))
    axs[1].set_ylim(np.min(temp_bank['mass2']),np.max(temp_bank['mass2']))
    axs[2].set_xlim(np.min(temp_bank['mass1']),np.max(temp_bank['mass1']))
    axs[2].set_ylim(np.min(temp_bank['mass2']),np.max(temp_bank['mass2']))

    axs[0].set_xlabel(r'$m_{1}$ $(M_{\odot})$', fontsize=fontsize, labelpad=10)
    axs[0].set_ylabel(r'$m_{2}$ $(M_{\odot})$', fontsize=fontsize, labelpad=10)
    axs[1].set_xlabel(r'$m_{1}$ $(M_{\odot})$', fontsize=fontsize, labelpad=10)
    #axs[1].set_ylabel(r'$m_{2}$ $(M_{\odot})$', fontsize=fontsize, labelpad=10)
    axs[2].set_xlabel(r'$m_{1}$ $(M_{\odot})$', fontsize=fontsize, labelpad=10)
    #axs[2].set_ylabel(r'$m_{2}$ $(M_{\odot})$', fontsize=fontsize, labelpad=10)

    axs[0].tick_params(axis='both', labelsize=ticksize)
    axs[1].tick_params(axis='both', labelsize=ticksize)
    axs[2].tick_params(axis='both', labelsize=ticksize)


    titles=[r'$\chi_{\regular{eff}}\leq 0$', r'$0\leq\chi_{\regular{eff}}<0.75$', r'$0.75\leq\chi_{\regular{eff}}$']

    axs[0].set_title(titles[0], fontsize=fontsize)
    axs[1].set_title(titles[1], fontsize=fontsize)
    axs[2].set_title(titles[2], fontsize=fontsize)

    cmap = plt.cm.OrRd
    alpha = 1
    lw = 1
    marker = 'o'
    colour = iter(cmap(np.linspace(0.,1.,psi_opts.shape[0]+1)))
    ims = []

    col = next(colour)

    conditions = [effspin<=0.,(0.<effspin)*(effspin<0.75),0.75<=effspin]

    for count,psi_opt in enumerate(psi_opts):
        col = next(colour)
        probs = np.abs(psi_opt)**2
        inds = probs>np.mean(probs)
        for i,condition in enumerate(conditions):
            sc = axs[i].scatter(temp_bank['mass1'][inds*condition], temp_bank['mass2'][inds*condition], color=col, lw=lw, marker=marker, cmap=cmap, vmin=-1, vmax=psi_opts.shape[0]-1, alpha=alpha, label=r'$\rho_{\mathregular{th}}=$'+str(int(SNRs[count])))

    leg = axs[-1].legend(loc='upper right', fontsize=fontsize)
    #leg.get_frame().set_linewidth(0.0)


    fig.savefig(outpath+'_'.join(infiles[-1].split('/')[-1].split('_')[:4])+'_snr_'+'_'.join(SNRs.astype(int).astype(str))+'_thrs_panel.png',bbox_inches='tight')

    exit()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Perform one QMF on GW150914 given a template bank")
    parser.add_argument('--infile', help="", type=str, nargs='+', required=True)
    parser.add_argument('--outpath', help="", type=str, required=True)
    parser.add_argument('--bank', help="", type=str, default='bank')
    parser.add_argument('--tempfile', help="", type=str, default=None)

    opt = parser.parse_args()
 
    main(opt.infile, opt.outpath, bank=opt.bank, tempfile=opt.tempfile)
