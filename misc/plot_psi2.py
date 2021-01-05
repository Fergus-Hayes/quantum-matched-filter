import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
import argparse
import h5py

def main(psi_file, outfile, bank_file, fontsize=24, ticksize=19):

    matplotlib.rc('text', usetex = True)
    matplotlib.rc('font', **{'family' : "sans-serif", 'size' : fontsize})

    psi_2 = np.load(psi_file)
    psi_2 = psi_2[:399952]

    full_bank = h5py.File(bank_file, 'r')

    bank = {}
    bank['mass1'] = np.array(full_bank['mass1'])   
    bank['mass2'] = np.array(full_bank['mass2'])
    bank['spin1z'] = np.array(full_bank['spin1z'])
    bank['spin2z'] = np.array(full_bank['spin2z'])
    bank['f_lower'] = np.array(full_bank['f_lower'])

    fig = plt.figure(figsize=(10,7))
    ax = plt.subplot(111)
    #plt.title('Figure 8')
    #sc=ax.scatter(np.append(bank['mass1'], 0.), np.append(bank['mass2'], 40.), c=np.append(np.abs(psi_2)**2, -.01), marker='.', lw=1, cmap=cm.Greys)
    #sc=ax.scatter(bank['mass1'], bank['mass2'], marker='x', cmap=cm.Reds)
    sc=ax.scatter(bank['mass1'], bank['mass2'], c=np.log(np.abs(psi_2)**2), marker='x', lw=1, cmap=cm.Reds, alpha=0.5)
    
    print(np.sum(np.mean(np.abs(psi_2)**2)>=np.abs(psi_2)**2))
    print(np.sum(np.mean(np.abs(psi_2)**2)<np.abs(psi_2)**2))
    print(np.sum(np.mean(np.abs(psi_2)**2)<np.abs(psi_2)**2)/len(psi_2))

    #ax.scatter(paras['mass1'][np.abs(psi_1_opt)**2>np.mean(np.abs(psi_1_opt)**2)], paras['mass2'][np.abs(psi_1_opt)**2>np.mean(np.abs(psi_1_opt)**2)], c='black', marker='.', lw=1, zorder=1)
    cb = plt.colorbar(sc, ax=[ax], label=r'$\textrm{ln}|\langle\psi_{2}\rangle|^{2}$')
    plt.xlabel(r'$m_{1}$', fontsize=fontsize)
    plt.ylabel(r'$m_{2}$', fontsize=fontsize)
    plt.xticks(fontsize=ticksize)
    plt.yticks(fontsize=ticksize)
    plt.savefig('mass_dists.png')
    plt.show()
    
    exit()

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

main('../psi_opt_None.npy', './psi2.png', './H1L1-BANK2HDF-1134450017-1202400.hdf')
