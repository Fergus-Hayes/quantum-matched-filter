import numpy as np
import gw_detections_functions as gwfn
import quantum_matched_filter_functions as qmffn
import argparse, time

np.random.seed(int(time.time()))

def main(Mqubits, Pqubits, tag='out', path='./output/', SNR_threshold=12., data_file='data/signal.npy', psd_file='data/psd.npy', template_file='data/template_bank.hdf'):
    
    Data = np.load(data_file)
    psd = np.load(psd_file)
    
    M = 2**Mqubits
    P = 2**Pqubits

    tag = str(M)+'_'+str(P)+'_'+str(SNR_threshold).replace('.','_')+'_'+tag

    measurement, psi_opt, Nmatches = qmffn.QMF(Data, psd, M, P, tag=tag, path=path, SNR_threshold=SNR_threshold, bankfunc=qmffn.get_paras, table=False, save_states=True, dtype='float64', temp_file=template_file)

    print(measurement)
    print(np.sum(np.abs(psi_opt)**2>np.mean(np.abs(psi_opt)**2)), Nmatches)

    np.save('./output/psi_opt_'+tag,psi_opt)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Perform one QMF on GW150914 given a template bank")
    parser.add_argument('--tag', help="Name tag", type=str, required=True)
    parser.add_argument('--Mq', help="Template qubits", type=int, default=8)
    parser.add_argument('--Pq', help="Precision qubits", type=int, default=7)
    parser.add_argument('--SNR-thr', help="SNR threshold", type=float, default=12.)
    parser.add_argument('--data-file', help="", type=str, default='data/signal.npy')
    parser.add_argument('--psd-file', help="", type=str, default='data/psd.npy')
    parser.add_argument('--temp-file', help="", type=str, default='data/template_bank.hdf')
    parser.add_argument('--path', help="", type=str, default='./output/')


    opt = parser.parse_args()
 
    main(opt.Mq, opt.Pq, tag=opt.tag, path=opt.path, SNR_threshold=opt.SNR_thr, data_file=opt.data_file, psd_file=opt.psd_file, template_file=opt.temp_file)
