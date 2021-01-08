import numpy as np
import matplotlib.pyplot as plt
import argparse, time, matplotlib

np.random.seed(int(time.time()))

def main(infile, outpath, fontsize=28, ticksize=22, figsize=(15,10)):

    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'

    psi_2 = np.load(infile)

    fig = plt.figure(figsize=figsize)
    ax = fig.gca()

    ax.scatter(np.arange(psi_2.shape[0]),np.sum(np.absolute(psi_2).T**2,axis=0), color='black', marker='.', lw=2)
    ax.set_xlabel(r'$P$', fontsize=fontsize)
    ax.set_ylabel(r'Probability of recovering correct template', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=ticksize)
    fig.savefig(outpath+'.'.join(infile.split('/')[-1].split('.')[:-1])+'_psi2.png')
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Perform one QMF on GW150914 given a template bank")
    parser.add_argument('--infile', help="", type=str, required=True)
    parser.add_argument('--outpath', help="", type=str, required=True)

    opt = parser.parse_args()
 
    main(opt.infile, opt.outpath)
