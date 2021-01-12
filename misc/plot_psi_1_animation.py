import numpy as np
import matplotlib.pyplot as plt
import argparse, time, matplotlib
import quantum_matched_filter_functions as qmffn
import matplotlib.animation as animation

np.random.seed(int(time.time()))

def main(infile, outpath, bank='bank', fontsize=28, ticksize=22, figsize=(14,11), tempfile=None):

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

    psi_1 = np.load(infile)

    match_ind = np.argmax(np.sum(np.abs(psi_1)**2,axis=1))
    nmatch_ind = np.argmin(np.sum(np.abs(psi_1)**2,axis=1))

    psi_match = psi_1[match_ind]
    psi_nmatch = psi_1[nmatch_ind]

    label_match = r'$\rho \geq \rho_{\mathregular{th}}$'
    label_nmatch = r'$\rho < \rho_{\mathregular{th}}$'
    ylabel = r'Probability amplitude'
    xlabel = r'$i$'

    M = psi_1.shape[0]

    temp_bank, M1, M2 = bankfunc(M, temp_file=tempfile)

    if True:
        fig2, axs2 = plt.subplots(2,1, figsize=figsize)
        axs2[0].set_xlabel(r'$m_{1}$', fontsize=fontsize)
        axs2[0].set_ylabel(r'$m_{2}$', fontsize=fontsize)

        axs2[0].tick_params(axis='both', labelsize=ticksize)
        axs2[1].tick_params(axis='both', labelsize=ticksize)

        axs2[1].plot(psi_1[match_ind], color='black', lw=2, label=label_match)
        axs2[1].plot(psi_1[nmatch_ind], ls='--', color='black', lw=2, label=label_nmatch)
        axs2[1].set_xlabel(r'$i$', fontsize=fontsize)
        axs2[1].set_ylabel(r'Probability amplitude', fontsize=fontsize)

        fig2.tight_layout()
       
        opt = np.argmax(psi_match)

        sc2 = axs2[0].scatter(temp_bank['mass1'], temp_bank['mass2'], c=psi_1[:,opt], marker='.', lw=3*(np.abs(psi_1[:,opt])**2)/np.max(np.abs(psi_1)**2), alpha=0.5, vmin=np.min(psi_1), vmax=np.max(psi_1), label=r'$p=$'+str(opt), cmap=matplotlib.cm.twilight_shifted)
        psi_minmax = np.sort([psi_match[opt],psi_nmatch[opt]])
        axs2[1].vlines(opt,ymin=psi_minmax[0],ymax=psi_minmax[1], color='black')

        leg2 = axs2[1].legend(loc='lower left', fontsize=fontsize)
        leg2.get_frame().set_linewidth(0.0)
        cblabel='Probability amplitude'
        cb2 = plt.colorbar(sc2, ax=[axs2])#[axs[0]])
        cb2.set_label(label=cblabel, fontsize=fontsize)
        cb2.ax.tick_params(labelsize=ticksize)
        fig2.savefig(outpath+'.'.join(infile.split('/')[-1].split('.')[:-1])+'_mass_ani_still.png')

    exit()

    fig, axs = plt.subplots(2,1, figsize=figsize)

    axs[0].set_xlabel(r'$m_{1}$', fontsize=fontsize)
    axs[0].set_ylabel(r'$m_{2}$', fontsize=fontsize)
    #axs[0].set_xlim(0.,150)

    axs[0].tick_params(axis='both', labelsize=ticksize)
    axs[1].tick_params(axis='both', labelsize=ticksize)

    axs[1].plot(psi_1[match_ind], color='black', lw=2, label=label_match)
    axs[1].plot(psi_1[nmatch_ind], ls='--', color='black', lw=2, label=label_nmatch)
    axs[1].set_xlabel(r'$i$', fontsize=fontsize)
    axs[1].set_ylabel(r'Probability amplitude', fontsize=fontsize)

    fig.tight_layout()

    ims = []
    for p in np.arange(psi_1.shape[1]):
        sc = axs[0].scatter(temp_bank['mass1'], temp_bank['mass2'], c=psi_1[:,p], marker='.', lw=3*(np.abs(psi_1[:,p])**2)/np.max(np.abs(psi_1)**2), alpha=0.5, vmin=np.min(psi_1), vmax=np.max(psi_1), label=r'$p=$'+str(p), cmap=matplotlib.cm.twilight_shifted)
        psi_minmax = np.sort([psi_match[p],psi_nmatch[p]])
        ims.append((sc,axs[1].vlines(p,ymin=psi_minmax[0],ymax=psi_minmax[1], color='black'),))

    leg = axs[1].legend(loc='lower left', fontsize=fontsize)
    leg.get_frame().set_linewidth(0.0)
    cblabel='Probability amplitude'
    cb = plt.colorbar(sc, ax=[axs])#[axs[0]])
    cb.set_label(label=cblabel, fontsize=fontsize)
    cb.ax.tick_params(labelsize=ticksize)
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=1800)
    im_ani = animation.ArtistAnimation(fig, ims, interval=50, repeat_delay=3000, blit=True)
    im_ani.save(outpath+'.'.join(infile.split('/')[-1].split('.')[:-1])+'_mass_ani.mp4', writer=writer)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Perform one QMF on GW150914 given a template bank")
    parser.add_argument('--infile', help="", type=str, required=True)
    parser.add_argument('--outpath', help="", type=str, required=True)
    parser.add_argument('--bank', help="", type=str, default='bank')
    parser.add_argument('--tempfile', help="", type=str, default=None)

    opt = parser.parse_args()
 
    main(opt.infile, opt.outpath, bank=opt.bank, tempfile=opt.tempfile)
