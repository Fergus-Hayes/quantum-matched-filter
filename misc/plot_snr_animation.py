import numpy as np
import matplotlib.pyplot as plt
import argparse, time, matplotlib
import quantum_matched_filter_functions as qmffn
import matplotlib.animation as animation
import scipy.stats as ss

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

    fig = plt.figure(figsize=figsize)
    ax = fig.gca()

    ax.set_xlabel(r'$m_{1}$ $(M_{\odot})$', fontsize=fontsize)
    ax.set_ylabel(r'$m_{2}$ $(M_{\odot})$', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=ticksize)

    if True:
        samples = np.load('data/GW150914_posterior_samples.npy')
        m1_ins, m2_ins = samples[0], samples[1] # m1, m2
        x = m1_ins[m1_ins>=m2_ins]
        y = m2_ins[m1_ins>=m2_ins]

        xmin, xmax = 4., 150#np.min(x), np.max(x)
        ymin, ymax = 4., 150#np.min(y), np.max(y)

        xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
        positions = np.vstack([xx.ravel(), yy.ravel()])
        values = np.vstack([x,y])
        kernel = ss.gaussian_kde(values)
        f = np.reshape(kernel(positions).T, xx.shape)
        cset = ax.contour(xx, yy, f, colors='k', alpha = 0.3, linewidths=3)
        cset.collections[0].set_label('GW150914 posterior')
        delta = .5
        ax.fill([xmin+delta,xmin+delta,xmax-delta],[ymin+delta+1,ymax-delta+1,ymax-delta+1], color='white', zorder=2)

    cmap = plt.cm.OrRd
    alpha = 1
    lw = 1
    marker = 'o'
    colour = iter(cmap(np.linspace(0.,1.,psi_opts.shape[0]+1)[::-1]))
    ims = []
    for i in np.arange(psi_opts.shape[0])[::-1]:
        col = next(colour)
        sc = ax.scatter(200, 4, color=col, lw=lw*10, marker=marker, label=r'$\rho_{\mathregular{th}}=$'+str(int(SNRs[i])), alpha=alpha)

    for count,psi_opt in enumerate(psi_opts):
        colours = np.zeros(M)
        for i in np.arange(count+1):
            probs = np.abs(psi_opts[i])**2
            colours = np.where(probs>np.mean(probs),i,colours)
        sc = ax.scatter(temp_bank['mass1'], temp_bank['mass2'], c=colours, lw=lw, marker=marker, cmap=cmap, vmin=-1, vmax=psi_opts.shape[0]-1, alpha=alpha)
        ims.append((sc,))
        ims.append((sc,))

    leg = fig.legend(loc='upper left', bbox_to_anchor=(0.15, 0.85), fontsize=fontsize)
    leg.get_frame().set_linewidth(0.0)

    fig.savefig(outpath+'_'.join(infiles[-1].split('/')[-1].split('_')[:4])+'_snr_'+'_'.join(SNRs.astype(int).astype(str))+'_thrs.png',bbox_inches='tight')

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
