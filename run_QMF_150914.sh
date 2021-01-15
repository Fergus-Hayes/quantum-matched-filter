#!/bin/bash

LOC="${PWD}/runs"
ENVLOC="/home/fergus.hayes/src/conda/anaconda3/bin/activate"
ENV="QMF"

M=17
P=6

SNR=8

CORES=3

echo 'Enter a word to describe the run and then hit ENTER:'

read TAG

RUN="M_${M}_P_${P}_${TAG}_bank"

mkdir -p ${LOC}/${RUN}

echo "universe = vanilla">"${LOC}/${RUN}/sample.sub"
echo "log = ${LOC}/${RUN}/sample.log">>"${LOC}/${RUN}/sample.sub"
echo "error = ${LOC}/${RUN}/sample.err">>"${LOC}/${RUN}/sample.sub"
echo "output = ${LOC}/${RUN}/sample.out">>"${LOC}/${RUN}/sample.sub"
echo "executable = ${PWD}/wrapper.sh">>"${LOC}/${RUN}/sample.sub"
echo "arguments = ${ENVLOC} ${ENV} ${PWD}/QMF_150914.py --Mq ${M} --Pq ${P} --tag ${TAG} --data-file data/signal.npy --SNR-thr ${SNR} --psd-file data/psd.npy --cores ${CORES} --path ${LOC}/${RUN}/">>"${LOC}/${RUN}/sample.sub"
echo "accounting_group = aluk.sim.o3.cbc.pe.lalinference">>"${LOC}/${RUN}/sample.sub"
echo "transfer_input_files = ">>"${LOC}/${RUN}/sample.sub"
echo "request_cpus = ${CORES}">>"${LOC}/${RUN}/sample.sub"
echo "request_memory = 2GB">>"${LOC}/${RUN}/sample.sub"
echo "request_disk = 2GB">>"${LOC}/${RUN}/sample.sub"
echo "queue 1">>"${LOC}/${RUN}/sample.sub"
condor_submit "${LOC}/${RUN}/sample.sub"

