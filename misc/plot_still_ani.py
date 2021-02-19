import numpy as np
import matplotlib.pyplot as plt
import argparse, time, matplotlib, weakref
import quantum_matched_filter_functions as qmffn
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

np.random.seed(int(time.time()))

def main(infile, outpath, bank='bank', fontsize=28, ticksize=22, figsize=(20,25), tempfile=None):

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

    effspin = (temp_bank['mass1']*temp_bank['spin1z'] + temp_bank['mass2']*temp_bank['spin2z'])/(temp_bank['mass1']+temp_bank['mass2'])

    if True:
        fig2 = plt.figure(figsize=figsize)
        gs = matplotlib.gridspec.GridSpec(3, 2, figure=fig2) 
        axs=[]
        axs.append(fig2.add_subplot(gs[:2,:], projection='3d'))
        axs.append(fig2.add_subplot(gs[2,:]))
        axs[0].set_ylabel(r'$m_{1}$ $(M_{\odot})$', fontsize=fontsize, labelpad=10)
        axs[0].set_xlabel(r'$\chi_{\regular{eff}}$', fontsize=fontsize, labelpad=10)
        axs[0].set_zlabel(r'$m_{2}$ $(M_{\odot})$', fontsize=fontsize, labelpad=10)

        axs[0].view_init(30, 45)
        axs[0].invert_xaxis()
        axs[0].tick_params(axis='both', labelsize=ticksize)
        
        opt = np.argmax(psi_match)
        lws = 3.+10.*(np.abs(psi_1[:,opt])**2)/np.max(np.abs(psi_1)**2)
        sc = axs[0].scatter(effspin, temp_bank['mass1'], temp_bank['mass2'], c=psi_1[:,opt], marker='.', s=lws, alpha=0.5, vmin=np.min(psi_1), vmax=np.max(psi_1), label=r'$p=$'+str(opt), cmap=matplotlib.cm.twilight_shifted)
        psi_minmax = np.sort([psi_match[opt],psi_nmatch[opt]])
        cblabel='Probability amplitude'
        cb = plt.colorbar(sc, ax=[axs])#[axs[0]])
        cb.set_label(label=cblabel, fontsize=fontsize)
        cb.ax.tick_params(labelsize=ticksize)

        axs[1].tick_params(axis='both', labelsize=ticksize)
        axs[1].plot(psi_1[match_ind], color='black', lw=2, label=label_match)
        axs[1].plot(psi_1[nmatch_ind], ls='--', color='black', lw=2, label=label_nmatch)
        axs[1].set_xlabel(r'Grover iterations', fontsize=fontsize)
        axs[1].set_ylabel(r'Probability amplitude', fontsize=fontsize)
       
        axs[1].vlines(opt,ymin=psi_minmax[0],ymax=psi_minmax[1], color='black')

        leg = axs[1].legend(loc='lower left', fontsize=fontsize)
        leg.get_frame().set_linewidth(0.0)
         
        fig2.savefig(outpath+'.'.join(infile.split('/')[-1].split('.')[:-1])+'_mass_effspin_ani_still.png',bbox_inches='tight')
        plt.close()

    fig = plt.figure(figsize=figsize)
    fig3 = plt.figure(figsize=figsize)

    #gs = matplotlib.gridspec.GridSpec(3, 2, figure=fig)
    axs=[]
    axs.append(fig.add_subplot(111, projection='3d'))
    axs.append(fig3.add_subplot(111))

    #axs.append(fig.add_subplot(gs[:2,:], projection='3d'))
    #axs.append(fig3.add_subplot(gs[2,:]))
    axs[0].set_ylabel(r'$m_{1}$ $(M_{\odot})$', fontsize=fontsize, labelpad=10)
    axs[0].set_xlabel(r'$\chi_{\regular{eff}}$', fontsize=fontsize, labelpad=10)
    axs[0].set_zlabel(r'$m_{2}$ $(M_{\odot})$', fontsize=fontsize, labelpad=10)

    #axs[0].view_init(30, 45)
    axs[0].invert_xaxis()
    axs[0].tick_params(axis='both', labelsize=ticksize)

    axs[1].tick_params(axis='both', labelsize=ticksize)
    axs[1].plot(psi_1[match_ind], color='black', lw=2, label=label_match)
    axs[1].plot(psi_1[nmatch_ind], ls='--', color='black', lw=2, label=label_nmatch)
    axs[1].set_xlabel(r'Grover iterations', fontsize=fontsize)
    axs[1].set_ylabel(r'Probability amplitude', fontsize=fontsize)

    leg = axs[1].legend(loc='lower left', fontsize=fontsize)
    leg.get_frame().set_linewidth(0.0)

    p=0
    axs[0].view_init(30, p*(360)/psi_1.shape[1])
    sc = axs[0].scatter(effspin, temp_bank['mass1'], temp_bank['mass2'], c=psi_1[:,p], marker='.', s=lws, alpha=0.5, vmin=np.min(psi_1), vmax=np.max(psi_1), label=r'$p=$'+str(p), cmap=matplotlib.cm.twilight_shifted)
    cblabel='Probability amplitude'
    #cb = plt.colorbar(sc, ax=[axs])#[axs[0]])
    #cb.set_label(label=cblabel, fontsize=fontsize)
    #cb.ax.tick_params(labelsize=ticksize)

    line1, line2 = [], []

    def init():
        p=0
        lws = 2+(8.*(np.abs(psi_1[:,p])**2)/np.max(np.abs(psi_1)**2))
        psi_minmax = np.sort([psi_match[p],psi_nmatch[p]])
        axs[0].view_init(30, 45+p*(360)/psi_1.shape[1])
        line1.append(axs[1].scatter(p, psi_match[p], c=psi_match[p], vmin=np.min(psi_1), vmax=np.max(psi_1), cmap=matplotlib.cm.twilight_shifted, lw=5, marker='o'))
        line2.append(axs[1].scatter(p, psi_nmatch[p], c=psi_nmatch[p], vmin=np.min(psi_1), vmax=np.max(psi_1), cmap=matplotlib.cm.twilight_shifted, lw=5, marker='o'))
        axs[0].scatter(effspin, temp_bank['mass1'], temp_bank['mass2'], c=psi_1[:,p], marker='.', s=lws, alpha=0.5, vmin=np.min(psi_1), vmax=np.max(psi_1), label=r'$p=$'+str(p), cmap=matplotlib.cm.twilight_shifted)
        return fig,

    def animate(p):
        lws = (2+(4.*(np.abs(psi_1[:,p])**2)/np.max(np.abs(psi_1)**2)))**2
        psi_minmax = np.sort([psi_match[p],psi_nmatch[p]])
        axs[0].view_init(30, 45+p*(360)/psi_1.shape[1])
        line1[-1].remove()
        line2[-1].remove()
        line1.append(axs[1].scatter(p, psi_match[p], c=psi_match[p], vmin=np.min(psi_1), vmax=np.max(psi_1), cmap=matplotlib.cm.twilight_shifted, lw=5, marker='o'))
        line2.append(axs[1].scatter(p, psi_nmatch[p], c=psi_nmatch[p], vmin=np.min(psi_1), vmax=np.max(psi_1), cmap=matplotlib.cm.twilight_shifted, lw=5, marker='o'))
        p=0
        axs[0].scatter(effspin, temp_bank['mass1'], temp_bank['mass2'], c=psi_1[:,p], marker='.', s=lws, alpha=0.5, vmin=np.min(psi_1), vmax=np.max(psi_1), label=r'$p=$'+str(p), cmap=matplotlib.cm.twilight_shifted)
        return fig,

    frames = psi_1.shape[1]
    im_ani = matplotlib.animation.FuncAnimation(fig, animate, init_func=init, frames=psi_1.shape[1], interval=100, blit=True) 

    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=1800)
    im_ani.save(outpath+'.'.join(infile.split('/')[-1].split('.')[:-1])+'_mass_ani_still.mp4', writer=writer)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Perform one QMF on GW150914 given a template bank")
    parser.add_argument('--infile', help="", type=str, required=True)
    parser.add_argument('--outpath', help="", type=str, required=True)
    parser.add_argument('--bank', help="", type=str, default='bank')
    parser.add_argument('--tempfile', help="", type=str, default=None)

    opt = parser.parse_args()
 
    main(opt.infile, opt.outpath, bank=opt.bank, tempfile=opt.tempfile)
