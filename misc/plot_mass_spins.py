import numpy as np
import matplotlib.pyplot as plt
import argparse, time, matplotlib
import quantum_matched_filter_functions as qmffn

np.random.seed(int(time.time()))

def main(infile, outpath, bank='bank', fontsize=28, ticksize=22, figsize=(15,10), tempfile=None):

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

    psi_2 = np.load(infile)

    M = psi_2.shape[0]

    temp_bank, M1, M2 = bankfunc(M,tempfile)

    fig = plt.figure(figsize=figsize)
    ax = fig.gca()

    sc=ax.scatter(temp_bank['mass1'], temp_bank['mass2'], c=np.log(np.abs(psi_2)**2).T, marker='x', lw=1, cmap=matplotlib.cm.Reds, alpha=0.5)
    cb = plt.colorbar(sc, ax=[ax], label=r'ln$|\langle\psi_{2}\rangle|^{2}$')
    ax.set_xlabel(r'$m_{1}$', fontsize=fontsize)
    ax.set_ylabel(r'$m_{2}$', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=ticksize)
    fig.savefig(outpath+'.'.join(infile.split('/')[-1].split('.')[:-1])+'_mass.png')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Perform one QMF on GW150914 given a template bank")
    parser.add_argument('--infile', help="", type=str, required=True)
    parser.add_argument('--outpath', help="", type=str, required=True)
    parser.add_argument('--bank', help="", type=str, default='bank')
    parser.add_argument('--tempfile', help="", type=str, default=None)

    opt = parser.parse_args()
 
    main(opt.infile, opt.outpath, bank=opt.bank, tempfile=opt.tempfile)
