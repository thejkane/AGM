#!/bin/bash

PT=./performance_test_wopq_rmat1
set -x
declare -a scalemap
scalemap[128]=31
scalemap[64]=30
scalemap[32]=29
scalemap[16]=28
scalemap[8]=27
scalemap[4]=26
scalemap[2]=25
scalemap[1]=24

for NODES in 1 4 8 16 32 64 128; do
for THREADS in 15; do
qsub - <<EOF
#!/bin/bash -l
#PBS -l nodes=$NODES:ppn=32
#PBS -l walltime=02:00:00
#PBS -N mis-ss-rmat1-wopq-$NODES-${scalemap[$NODES]}
#PBS -V
#PBS -j oe
#PBS -m abe
#PBS -M thejkane@indiana.edu
#PBS -S /bin/bash

set -x

export MPICH_NEMESIS_ASYNC_PROGRESS=SC
export MPICH_MAX_THREAD_SAFETY=multiple
ulimit -c unlimited

cd \$PBS_O_WORKDIR
echo 'Running script\n'
for i in 1; do
    aprun -b -n $(( 2 * $NODES )) -N 2 -d 15 -r 1 $PT --poll-task 1 --threads $THREADS --scale ${scalemap[$NODES]} --degree 16 --num-sources 3 --coalescing-size 160000,180000 --flush 20 --ds threadq --eager-limit 10 --max-weight 100 --receive-depth 4 --priority_coalescing_size 43000 --with-no-reductions --without-per-thrad-reductions --run_ss_mis
done
EOF

done
done
