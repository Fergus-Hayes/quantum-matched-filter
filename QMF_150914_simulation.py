import numpy as np
from tqdm import tqdm
import gw_detections_functions as gwfn
import quantum_matched_filter_functions as qmffn
import argparse, time
import matplotlib.pyplot as plt
import argparse, time, matplotlib, os

np.random.seed(int(time.time()))

def main(infile, outpath, tag, fontsize=28, ticksize=22, figsize=(12,8), dtype='float64', all_trials=[1], runs=1000):

    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'

    SNRs = [float(infile.split('/')[-1].split('_')[4])]

    snr_inds = np.argsort(SNRs)
    SNRs = np.array(SNRs)[snr_inds]

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(figsize[0],figsize[1]))

    ax.set_xlabel(r'Number of Grover operations', fontsize=fontsize)
    ax.set_ylabel(r'Frequency', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=ticksize, top=False, right=False)

    lw=3
    ms=6

    #cmap = plt.cm.
    colors = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3','#999999', '#e41a1c', '#dede00']
    #colors = ['#377eb8', '#e41a1c', '#ff7f00']
    #iter(cmap(np.linspace(0,1,len(infiles)+1)[::-1]))
    #col = next(colors)

    labels = []

    for i in np.arange(len(all_trials)):
        
        trials = all_trials[i]
        if trials>0:
            labels.append(r'$B=$'+str(trials))
        if trials==0:
            labels.append(r'$B=$'+'\u221e')
        if trials<0.:
            labels.append(r'$B=2^{p}/k_{\regular{obs}}$')

        filename = infile.split('/')[-1]
        snr_filename = 'snrs_'+'_'.join(filename.split('_')[2:])

        if os.path.isfile('/'.join(infile.split('/')[:-1])+'/'+snr_filename):
            snrs = np.load('/'.join(infile.split('/')[:-1])+'/'+snr_filename)

        SNR_threshold = float(infile.split('/')[-1].split('_')[4])

        prob = np.load(infile)

        P = prob.shape[0]
        M = prob.shape[1]
        prob = np.sum(np.abs(prob)**2,axis=1)
        prob = prob/np.sum(prob)

        w = np.where(snrs>=SNR_threshold,-1.,1.)/np.sqrt(M)

        run_ops = []
        for run in tqdm(np.arange(runs)):
            switch = True
            count = 0
            ops = 0
            while switch:
                count+=1
                ops+=(P-1)

                b_obs = 0
                while b_obs==0:
                    b_obs = np.random.choice(P, 1, p=prob)
                k_obs = np.round(((np.pi/4.)*np.sqrt(M/b_obs))-0.5).astype(int)#np.round((P/(4*b_obs))-0.5).astype(int)
                t_obs = np.round(M*np.sin(b_obs*np.pi/P)**2).astype(int)
                if t_obs == 0:
                    t_obs=np.array([1])
                psi_k = qmffn.iquantum_counting0(w, np.ones(M).astype(dtype)/np.sqrt(M), k_obs, dtype)
 
                prob_k = np.abs(psi_k)**2
                prob_k = prob_k/np.sum(prob_k)

                if trials < 0:
                    trials_ = int(P/k_obs)
                    for trial in np.arange(1,trials_+1):
                        ops+=k_obs
                        temp_ind = np.random.choice(M, 1, p=prob_k)
                        if snrs[temp_ind]>SNR_threshold:
                            switch=False
                            break

                elif trials>0:
                    trials_ = trials
                    for trial in np.arange(1,trials_+1):
                        ops+=k_obs
                        temp_ind = np.random.choice(M, 1, p=prob_k)
                        if snrs[temp_ind]>SNR_threshold:
                            switch=False
                            break
                
                elif trials==0:
                    switch2 = True
                    while switch2:
                        ops+=k_obs
                        temp_ind = np.random.choice(M, 1, p=prob_k)
                        if snrs[temp_ind]>SNR_threshold:
                            switch=False
                            switch2=False
                            break
            
            run_ops.append(ops)
        
        run_ops = np.array(run_ops)
        
        np.save(outpath+'simulation_out_'+str(i)+'_'+tag,run_ops)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Perform one QMF on GW150914 given a template bank")
    parser.add_argument('--infile', help="", type=str, required=True)
    parser.add_argument('--outpath', help="", type=str, required=True)
    parser.add_argument('--tag', help="", type=str, required=True)
    parser.add_argument('--trials', help="", type=int, nargs='+', required=True)

    opt = parser.parse_args()

    main(opt.infile, opt.outpath, opt.tag, all_trials=opt.trials)
