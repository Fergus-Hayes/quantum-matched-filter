#!/bin/bash

TARGET=${1}
FTARGET="formated_${TARGET}"

cp ${TARGET} ${FTARGET}

sed -i 's/{max}/{\\text{max}}/g' ${FTARGET}
sed -i 's/_{th}/_{\\text{th}}/g' ${FTARGET}
sed -i 's/_{obs}/_{\\text{obs}}/g' ${FTARGET}
sed -i 's/\rho_{\\text{th}}/\\rho_{\\text{thr}}/g' ${FTARGET}
