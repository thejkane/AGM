#!/bin/bash
#
#PBS -m abe
#PBS -q cpu 
#PBS -N agm_weak_2_64
#PBS -l walltime=24:00:00
#PBS -l nodes=64:ppn=32
#PBS -j oe
#PBS -S /bin/bash

cd $PBS_O_WORKDIR


export MPICH_NEMESIS_ASYNC_PROGRESS=SC
export MPICH_MAX_THREAD_SAFETY=multiple

export ALGO=AGMWeak
export Q_NAME=cpu
export NUM_NODES=64
export JOB_NAME=agm_weak_2_64
export WALL_TIME="24:00:00"
export MAX_THREADS=32
export NUM_SOURCES=24
export POLL_TASKS=4
export FLUSH=22
export RECV_DEPTH=4
export COALESCING_SIZES=43000,86000
export DELTAS=1,2,3,5
export KLEVELS=1,2,3,5
export DATA_STRUCTURES=ms_nodeq,ms_numaq,threadq,buffer,nodeq,numaq
export REDUCTION_CACHE_SIZES=18
export PERF_BINARY=./performance_test

let NUM_THREADS=$MAX_THREADS-1

declare -a node_array

node_array[30]=64
node_array[29]=32
node_array[28]=16
node_array[27]=8
node_array[26]=4
node_array[25]=2

for scale in 30 29 28 27 26 25
do
    echo "######### Running Scale : $scale, Nodes : ${node_array[${scale}]}, Threads : $NUM_THREADS, Delta : $DELTAS, Coalescing : $COALESCING_SIZES, Poll : $POLL_TASKS, Flush : $FLUSH, Depth : $RECV_DEPTH #######"
    echo "######### Binary : $PERF_BINARY, Data Structures : $DATA_STRUCTURES, Reduction Cache Sizes : $REDUCTION_CACHE_SIZES #######"
    echo "#######################[DELTA STEPPING]###########################"
#    aprun -b -n ${node_array[${scale}]}  -N 1 -d $NUM_THREADS -r 1 $PERF_BINARY --poll-task $POLL_TASKS --threads $NUM_THREADS --scale $scale --degree 16 --num-sources $NUM_SOURCES --coalescing-size $COALESCING_SIZES --flush $FLUSH --delta $DELTAS --ds $DATA_STRUCTURES --eager-limit 10 --max-weight 100 --receive-depth $RECV_DEPTH --priority_coalescing_size 43000 --with-no-reductions --without-per-thread-reductions --run_ds
    
    echo "\n"
    echo "#######################[KLA]###########################"
#    aprun -b -n ${node_array[${scale}]}  -N 1 -d $NUM_THREADS -r 1 $PERF_BINARY --poll-task $POLL_TASKS --threads $NUM_THREADS --scale $scale --degree 16 --num-sources $NUM_SOURCES --coalescing-size $COALESCING_SIZES --flush $FLUSH --ds $DATA_STRUCTURES --klevel $KLEVELS --eager-limit 10 --max-weight 100 --receive-depth $RECV_DEPTH --priority_coalescing_size 43000 --with-no-reductions --without-per-thread-reductions --run_kla

    echo "\n"
    echo "#######################[Chaotic]###########################"
#    aprun -b -n ${node_array[${scale}]}  -N 1 -d $NUM_THREADS -r 1 $PERF_BINARY --poll-task $POLL_TASKS --threads $NUM_THREADS --scale $scale --degree 16 --num-sources $NUM_SOURCES --coalescing-size $COALESCING_SIZES --flush $FLUSH --ds ms_nodeq,nodeq,ms_numaq,numaq,threadq --eager-limit 10 --max-weight 100 --receive-depth $RECV_DEPTH --priority_coalescing_size 43000 --with-no-reductions --without-per-thread-reductions --run_dc 
done
