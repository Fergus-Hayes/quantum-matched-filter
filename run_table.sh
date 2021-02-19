#!/bin/bash

M=18
m=$((2**${M}))

declare -a SNRS=("8" "12" "16" "18")
declare -a PS=("12" "13" "14")

ENVLOC="/home/fergus.hayes/src/conda/anaconda3/bin/activate"
ENV="QMF"
CORES=1

#declare -a SNRS=("18")
#declare -a PS=("12")

LOC="${PWD}/runs"

for P in "${PS[@]}" ; do
for SNR in "${SNRS[@]}" ; do

TAG="loading_snr"
RUN="M_${M}_P_${P}_${SNR}_${TAG}_multi"

mkdir -p ${LOC}/${RUN}/

p=$((2**${P}))

cp data/SNRs_signal_spins.npy ${LOC}/${RUN}/snrs_${m}_${p}_${SNR}_0_${TAG}.npy

if [ "$(echo $HOSTNAME)" == 'cl8' ]; then

echo "P: ${P} SNR: ${SNR}"
echo "universe = vanilla">"${LOC}/${RUN}/sample.sub"
echo "log = ${LOC}/${RUN}/sample.log">>"${LOC}/${RUN}/sample.sub"
echo "error = ${LOC}/${RUN}/sample.err">>"${LOC}/${RUN}/sample.sub"
echo "output = ${LOC}/${RUN}/sample.out">>"${LOC}/${RUN}/sample.sub"
echo "executable = ${PWD}/wrapper.sh">>"${LOC}/${RUN}/sample.sub"
echo "arguments = ${ENVLOC} ${ENV} ${PWD}/QMF_150914.py --Mq ${M} --Pq ${P} --bank bank --tag ${TAG} --SNR-thr ${SNR} --cores ${CORES} --data-file data/signal.npy --psd-file data/psd.npy --path ${LOC}/${RUN}/">>"${LOC}/${RUN}/sample.sub"
echo "accounting_group = aluk.sim.o3.cbc.pe.lalinference">>"${LOC}/${RUN}/sample.sub"
echo "transfer_input_files = ">>"${LOC}/${RUN}/sample.sub"
echo "request_cpus = ${CORES}">>"${LOC}/${RUN}/sample.sub"
echo "request_memory = 10GB">>"${LOC}/${RUN}/sample.sub"
echo "request_disk = 2GB">>"${LOC}/${RUN}/sample.sub"
echo "queue 1">>"${LOC}/${RUN}/sample.sub"
condor_submit "${LOC}/${RUN}/sample.sub"

else
bash ${PWD}/wrapper.sh ${ENVLOC} ${ENV} ${PWD}/QMF_150914.py --Mq ${M} --Pq ${P} --bank bank --tag ${TAG} --SNR-thr ${SNR} --cores ${CORES} --data-file data/signal.npy --psd-file data/psd.npy --path ${LOC}/${RUN}/
fi

done
done
