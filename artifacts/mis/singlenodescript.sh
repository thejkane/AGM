#!/bin/bash -l
#PBS -l nodes=1:ppn=32
#PBS -l walltime=03:00:00
#PBS -N mis-rmat1-no-pq-1
#PBS -V
#PBS -j oe
#PBS -m abe
#PBS -M thejkane@indiana.edu
#PBS -S /bin/bash


PT=./performance_test_wopq_rmat1
set -x
declare -a scalemap
scalemap[32]=24
scalemap[16]=23
scalemap[8]=22
scalemap[4]=21
scalemap[2]=20
scalemap[1]=19

set -x

#export MPICH_NEMESIS_ASYNC_PROGRESS=SC
export MPICH_MAX_THREAD_SAFETY=multiple
ulimit -c unlimited

cd $PBS_O_WORKDIR
echo 'Running script\n'
for i in 16 8 4 2 1; do
    aprun -b -n 1 -N 1 -d $i $PT --poll-task 1 --threads $i --scale ${scalemap[$i]} --degree 16 --num-sources 3 --coalescing-size 40000 --flush 20 --ds threadq --eager-limit 10 --max-weight 100 --receive-depth 4 --priority_coalescing_size 43000 --with-no-reductions --without-per-thrad-reductions --run_ss_mis
done
