import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
import argparse

def main(psi_file, outfile, fontsize=24, ticksize=19):

    matplotlib.rc('text', usetex = True)
    matplotlib.rc('font', **{'family' : "sans-serif", 'size' : fontsize})

    psi_1 = np.load(psi_file)

    psi_1 = psi_1[:399952]

    print(psi_1.shape)

    fig = plt.figure(figsize=(15,10))
    ax = fig.gca()

    #matches=np.sum(np.sum(np.abs(psi_1)**2,axis=1)>np.mean(np.sum(np.abs(psi_1)**2,axis=1)))
    #total=len(np.sum(np.abs(psi_1)**2,axis=1))
    #print(np.log2((2*5*np.pi/7)*np.sqrt(total/matches)))
    #exit()

    match_ind = np.argmax(np.sum(np.abs(psi_1)**2,axis=1))
    nmatch_ind = np.argmin(np.sum(np.abs(psi_1)**2,axis=1))

    psi_match = psi_1[match_ind]
    psi_nmatch = psi_1[nmatch_ind]

    ax.plot(psi_1[match_ind], label=r'$\rho \ge \rho_{\normalsize \textrm{th}}$', color='black', lw=2)
    ax.plot(psi_1[nmatch_ind], label=r'$\rho < \rho_{\normalsize \textrm{th}}$', ls='--', color='black', lw=2)

    ax.set_xlabel(r'$i$', fontsize=fontsize)
    ax.set_ylabel(r'\textrm{Probability amplitude}', fontsize=fontsize)
    leg = ax.legend(fontsize=fontsize)#, loc='upper right')
    leg.get_frame().set_linewidth(0.0)
    fig.savefig(outfile)

    fig = plt.figure(figsize=(15,10))
    ax = fig.gca()
    surf = ax.imshow(psi_1, cmap=cm.seismic, aspect='auto')
    cbar = plt.colorbar(surf)
    cbar.set_label(r'\textrm{Probability amplitude}')
    ax.set_xlabel(r'$a$', fontsize=fontsize ,labelpad=20)
    ax.set_ylabel(r'$i$', fontsize=fontsize, labelpad=20)
    #ax.set_zlabel(r'\textrm{Probability amplitude}', fontsize=fontsize ,labelpad=25)
    ax.tick_params(labelsize=ticksize, pad=12)
    fig.savefig('interp.png')

    fig3 = plt.figure(figsize=(15,10))
    ax3 = fig3.gca()
    surf = ax3.imshow(psi_1, cmap=cm.seismic, interpolation='none', aspect='auto')
    cbar = plt.colorbar(surf)
    cbar.set_label(r'\textrm{Probability amplitude}')
    ax3.set_xlabel(r'$a$', fontsize=fontsize ,labelpad=20)
    ax3.set_ylabel(r'$i$', fontsize=fontsize, labelpad=20)
    #ax.set_zlabel(r'\textrm{Probability amplitude}', fontsize=fontsize ,labelpad=25)
    ax3.tick_params(labelsize=ticksize, pad=12)
    fig3.savefig('ninterp.png')

    fig2 = plt.figure(figsize=(15,10))
    ax2 = fig2.gca()

    ax2.plot(np.sum(np.abs(psi_1)**2,axis=1))
    ax2.set_xlabel(r'$a$', fontsize=fontsize ,labelpad=20)
    ax2.set_ylabel(r'$\textrm{Probability}$', fontsize=fontsize ,labelpad=20)
    fig2.savefig('Nstates.png')

main('../psi1_in_None.npy', './psi1.png')
