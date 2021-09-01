import numpy as np
from tqdm import tqdm
import gw_detections_functions as gwfn
import quantum_matched_filter_functions as qmffn
import argparse, time
import matplotlib.pyplot as plt
import argparse, time, matplotlib, os, sys

np.random.seed(int(time.time()))

outpath = '/home/fergus.hayes/public_html/LVC/Misc/QMF/'

P = 2**11
N = 28*4096
M = 2**17 

SD_cost = 1.#(P-1)#*((N*np.log(N))+np.log(M))

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
fontsize=28
ticksize=22
figsize=(12,8)

colors = ['blue', '#e41a1c', 'blue']#['#377eb8', '#e41a1c', '#999999']#'#ff7f00']

#trials1 = np.load('output/simulation_out_cost_1.npy')/SD_cost
#trials0 = np.load('output/simulation_out_cost_0.npy')/SD_cost
trials1 = np.load('output/simulation_out_1.npy')#/(P-1)
trials0 = np.load('output/simulation_out_0.npy')#/(P-1)
#trials3 = np.load('output/simulation_-1.npy') - ((2**11) - 1)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(figsize[0],figsize[1]))
labels = [r'trials$=1$', r'trials$=$'+'\u221e', r'trials$=(2^{p}-1)/k_{\regular{obs}}$']

#bins = np.linspace(np.log10(1),np.log10(100),41)#np.log10(np.logspace(8, 10, 41)/SD_cost)

bins = np.linspace(np.log10(1*(P-1)),np.log10(100*(P-1)),41)
bins = np.power(10.,bins)
ax.hist(trials1, color=colors[0], alpha=0.5, bins=bins, label=labels[0])
ax.hist(trials0, color=colors[1], alpha=0.5, bins=bins, label=labels[1])
ax.axvline(np.mean(trials1), color=colors[0], ls='--', lw=3)
ax.axvline(np.mean(trials0), color=colors[1], ls='--', lw=3)
ax.axvline(M/SD_cost, color='black', lw=2, ls=':')

#ax.axvline((N*M*np.log(N))/SD_cost, color='black', lw=2, ls=':')

print(np.mean(trials1), np.mean(trials0))
print(M/SD_cost)

ax.set_xlabel(r'Evaluations of $f$', fontsize=fontsize)
ax.set_ylabel(r'Frequency', fontsize=fontsize)
#ax.set_yscale('log')
ax.set_xscale('log')
ax.tick_params(axis='both', labelsize=ticksize, top=False, right=False)

fig.savefig(outpath+'simulation_'+sys.argv[1]+'_scenarios.png', bbox_inches='tight')
