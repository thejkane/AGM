#!/bin/bash
#
#PBS -m abe
#PBS -N initial
#PBS -l walltime=01:00:00
#PBS -l nodes=2:ppn=32
#PBS -j oe
#PBS -S /bin/bash

cd bin

export MPICH_NEMESIS_ASYNC_PROGRESS=SC
export MPICH_MAX_THREAD_SAFETY=multiple

aprun -b -n 1  -N 1 -d 31 -r 1 ./performance_test --poll-task 1 --threads 31 --scale 10 --degree 16 --num-sources 3 --coalescing-size 42000 --flush 18 --eager-limit 10 --max-weight 100 --receive-depth 4 --priority_coalescing_size 43000 --with-no-reductions --without-per-thrad-reductions --run_cc --run_dd_cc --run_sv_cc --run_sv_ps_cc
aprun -b -n 1  -N 1 -d 31 -r 1 ./performance_test --poll-task 1 --threads 31 --scale 11 --degree 16 --num-sources 3 --coalescing-size 42000 --flush 18 --eager-limit 10 --max-weight 100 --receive-depth 4 --priority_coalescing_size 43000 --with-no-reductions --without-per-thrad-reductions --run_cc --run_dd_cc --run_sv_cc --run_sv_ps_cc
aprun -b -n 1  -N 1 -d 31 -r 1 ./performance_test --poll-task 1 --threads 31 --scale 12 --degree 16 --num-sources 3 --coalescing-size 42000 --flush 18 --eager-limit 10 --max-weight 100 --receive-depth 4 --priority_coalescing_size 43000 --with-no-reductions --without-per-thrad-reductions --run_cc --run_dd_cc --run_sv_cc --run_sv_ps_cc
aprun -b -n 1  -N 1 -d 31 -r 1 ./performance_test --poll-task 1 --threads 31 --scale 13 --degree 16 --num-sources 3 --coalescing-size 42000 --flush 18 --eager-limit 10 --max-weight 100 --receive-depth 4 --priority_coalescing_size 43000 --with-no-reductions --without-per-thrad-reductions --run_cc --run_dd_cc --run_sv_cc --run_sv_ps_cc
aprun -b -n 1  -N 1 -d 31 -r 1 ./performance_test --poll-task 1 --threads 31 --scale 14 --degree 16 --num-sources 3 --coalescing-size 42000 --flush 18 --eager-limit 10 --max-weight 100 --receive-depth 4 --priority_coalescing_size 43000 --with-no-reductions --without-per-thrad-reductions --run_cc --run_dd_cc --run_sv_cc --run_sv_ps_cc
aprun -b -n 1  -N 1 -d 31 -r 1 ./performance_test --poll-task 1 --threads 31 --scale 15 --degree 16 --num-sources 3 --coalescing-size 42000 --flush 18 --eager-limit 10 --max-weight 100 --receive-depth 4 --priority_coalescing_size 43000 --with-no-reductions --without-per-thrad-reductions --run_cc --run_dd_cc --run_sv_cc --run_sv_ps_cc
aprun -b -n 1  -N 1 -d 31 -r 1 ./performance_test --poll-task 1 --threads 31 --scale 16 --degree 16 --num-sources 3 --coalescing-size 42000 --flush 18 --eager-limit 10 --max-weight 100 --receive-depth 4 --priority_coalescing_size 43000 --with-no-reductions --without-per-thrad-reductions --run_cc --run_dd_cc --run_sv_cc --run_sv_ps_cc
aprun -b -n 1  -N 1 -d 31 -r 1 ./performance_test --poll-task 1 --threads 31 --scale 17 --degree 16 --num-sources 3 --coalescing-size 42000 --flush 18 --eager-limit 10 --max-weight 100 --receive-depth 4 --priority_coalescing_size 43000 --with-no-reductions --without-per-thrad-reductions --run_cc --run_dd_cc --run_sv_cc --run_sv_ps_cc
aprun -b -n 1  -N 1 -d 31 -r 1 ./performance_test --poll-task 1 --threads 31 --scale 18 --degree 16 --num-sources 3 --coalescing-size 42000 --flush 18 --eager-limit 10 --max-weight 100 --receive-depth 4 --priority_coalescing_size 43000 --with-no-reductions --without-per-thrad-reductions --run_cc --run_dd_cc --run_sv_cc --run_sv_ps_cc
aprun -b -n 1  -N 1 -d 31 -r 1 ./performance_test --poll-task 1 --threads 31 --scale 19 --degree 16 --num-sources 3 --coalescing-size 42000 --flush 18 --eager-limit 10 --max-weight 100 --receive-depth 4 --priority_coalescing_size 43000 --with-no-reductions --without-per-thrad-reductions --run_cc --run_dd_cc --run_sv_cc --run_sv_ps_cc
aprun -b -n 1  -N 1 -d 31 -r 1 ./performance_test --poll-task 1 --threads 31 --scale 20 --degree 16 --num-sources 3 --coalescing-size 42000 --flush 18 --eager-limit 10 --max-weight 100 --receive-depth 4 --priority_coalescing_size 43000 --with-no-reductions --without-per-thrad-reductions --run_cc --run_dd_cc --run_sv_cc --run_sv_ps_cc


