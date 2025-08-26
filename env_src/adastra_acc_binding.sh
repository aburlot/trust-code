#!/bin/bash
set -eu

LOCAL_RANK_INDEX="${SLURM_LOCALID}"

function Adastra_MI300_4TasksWith24ThreadsAnd1GPU() {
    AFFINITY_NUMACTL=('0-23' '24-47' '48-71' '72-95')
    AFFINITY_GPU=('0' '1' '2' '3')
}

function Adastra_MI250_8TasksWith8ThreadsAnd1GPU() {
    AFFINITY_NUMACTL=('48-55' '56-63' '16-23' '24-31' '0-7' '8-15' '32-39' '40-47')
    AFFINITY_GPU=('0' '1' '2' '3' '4' '5' '6' '7')
}

function Adastra_GENOA_24TasksWith8Threads() {
    AFFINITY_NUMACTL=('0-7' '8-15' '16-23' '24-31' '32-39' '40-47' '48-55' '56-63' '64-71' '72-79' '80-87' '88-95' '96-103' '104-111' '112-119' '120-127' '128-135' '136-143' '144-151' '152-159' '160-167' '168-175' '176-183' '184-191')
}

function Adastra_GENOA_48TasksWith4Threads() {
    AFFINITY_NUMACTL=('0-3' '4-7' '8-11' '12-15' '16-19' '20-23' '24-27' '28-31' '32-35' '36-39' '40-43' '44-47' '48-51' '52-55' '56-59' '60-63' '64-67' '68-71' '72-75' '76-79' '80-83' '84-87' '88-91' '92-95' '96-99' '100-103' '104-107' '108-111' '112-115' '116-119' '120-123' '124-127' '128-131' '132-135' '136-139' '140-143' '144-147' '148-151' '152-155' '156-159' '160-163' '164-167' '168-171' '172-175' '176-179' '180-183' '184-187' '188-191')
}

export MPICH_OFI_NIC_POLICY="NUMA"
export OMP_PROC_BIND="TRUE"

Adastra_MI250_8TasksWith8ThreadsAnd1GPU
# Adastra_GENOA_24TasksWith8Threads
# Adastra_GENOA_48TasksWith4Threads

CPU_SET="${AFFINITY_NUMACTL[$((${LOCAL_RANK_INDEX} % ${#AFFINITY_NUMACTL[@]}))]}"
if [ ! -z ${AFFINITY_GPU+x} ]; then
    GPU_SET="${AFFINITY_GPU[$((${LOCAL_RANK_INDEX} % ${#AFFINITY_GPU[@]}))]}"
    export ROCR_VISIBLE_DEVICES="${GPU_SET}"
fi
echo physcpubind="${CPU_SET}"
exec numactl --localalloc --physcpubind="${CPU_SET}" -- "${@}"
