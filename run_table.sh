#!/bin/bash

M=19
m=$((2**${M}))

#declare -a SNRS=("8" "12" "16" "18")
#declare -a PS=("9" "10")

declare -a SNRS=("18")
declare -a PS=("10")

LOC="${PWD}/runs"

for P in "${PS[@]}" ; do
for SNR in "${SNRS[@]}" ; do

TAG="loading_snr"
RUN="M_17_P_${P}_${SNR}_${TAG}_multi"

mkdir -p ${LOC}/${RUN}/

p=$((2**${P}))

cp data/SNRs_signal_spins.npy ${LOC}/${RUN}/snrs_${m}_${p}_${SNR}_0_${TAG}.npy

echo "P: ${P} SNR: ${SNR}"

bash wrapper.sh ~/src/anaconda3/bin/activate QMF QMF_150914.py --Mq ${M} --Pq ${P} --bank bank --tag ${TAG} --SNR-thr ${SNR} --cores 3 --data-file data/signal.npy --psd-file data/psd.npy --path ${LOC}/${RUN}/

done
done
