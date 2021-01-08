import numpy as np
import matplotlib.pyplot as plt
import argparse, time, matplotlib

np.random.seed(int(time.time()))

def main(infile, outpath, fontsize=28, ticksize=22, figsize=(15,10)):

    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'

    psi_1 = np.load(infile)

    fig = plt.figure(figsize=figsize)
    ax = fig.gca()

    match_ind = np.argmax(np.sum(np.abs(psi_1)**2,axis=1))
    nmatch_ind = np.argmin(np.sum(np.abs(psi_1)**2,axis=1))

    psi_match = psi_1[match_ind]
    psi_nmatch = psi_1[nmatch_ind]

    label_match = r'$\rho \geq \rho_{\mathregular{th}}$'
    label_nmatch = r'$\rho < \rho_{\mathregular{th}}$'
    ylabel = r'Probability amplitude'
    xlabel = r'$i$'

    ax.plot(psi_1[match_ind], color='black', lw=2, label=label_match)
    ax.plot(psi_1[nmatch_ind], ls='--', color='black', lw=2, label=label_nmatch)

    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=ticksize)
    leg = ax.legend(fontsize=fontsize)#, loc='upper right')
    leg.get_frame().set_linewidth(0.0)
    fig.savefig(outpath+'.'.join(infile.split('/')[-1].split('.')[:-1])+'_psi1.png')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Perform one QMF on GW150914 given a template bank")
    parser.add_argument('--infile', help="", type=str, required=True)
    parser.add_argument('--outpath', help="", type=str, required=True)

    opt = parser.parse_args()
 
    main(opt.infile, opt.outpath)
