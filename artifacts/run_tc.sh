cd bin

#export MPICH_NEMESIS_ASYNC_PROGRESS=SC
export MPICH_NEMESIS_ASYNC_PROGRESS=MC
export MPICH_NEMESIS_ON_NODE_ASYNC_OPT=1
export MPICH_MAX_THREAD_SAFETY=multiple
export MPICH_CRAY_OPT_THREAD_SYNC=1
export MPICH_STATS_VERBOSITY=3
export MPICH_OPTIMIZED_MEMCPY=1

#export MPICH_GNI_USE_UNASSIGNED_CPUS=enabled
#export OMP_NUM_THREADS=16
#export PAT_RT_HWPC=2
export PAT_RT_PERFCTR=2
export LD_LIBRARY_PATH=/home/users/p02119/development/cds-2.1.0/bin/cc-amd64-linux-64/:$LD_LIBRARY_PATH


export CORES=4
export SCALE=23

aprun -b -n 128 -N 2 -d 22 ./tc_family_ordered --no-weight --threads 16 --file /lus/scratch/p02119/graphs/twitter/out.twitter --read-graph --id-distribution horizontal --warm-up-iterations 0 --preprocess --vertices 41652230 --edges 1468365182 --graph-type twitter --iterations 1 --run_tc_sucsuc --pred-block-size 100000,200000,400000 --coalescing-size 20000000,25000000,30000000 --suc-block-size 100000,200000,400000 --expected-triangles 34824916864 --verify

#aprun -n 2 -N 1 -d 16 ./tc_family --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 0 --preprocess --iterations 1 --scale 25 --run_tc_sucsuc --pred-block-size 1000 --coalescing-size 4000000 --suc-block-size 1000
#aprun -n 1 -N 1 -d 16 ./tc_family --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 0 --preprocess --iterations 1 --scale 25 --run_tc_sucsuc --pred-block-size 1000 --coalescing-size 4000000 --suc-block-size 1000
#aprun -n 4 -N 1 -d 22 ./tc_family --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 0 --preprocess --iterations 1 --scale 23 --run_tc_sucsuc --pred-block-size 1000 --coalescing-size 4000000 --suc-block-size 1000
#aprun -n 4 -N 1 -d 22 ./tc_family_ordered --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 0 --preprocess --iterations 1 --scale 23 --run_tc_sucsuc --pred-block-size 1000 --coalescing-size 4000000 --suc-block-size 1000
#aprun -n 4 -N 1 -d 22 ./tc_family_ordered --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 1 --preprocess --iterations 3 --scale 20 --run_tc_sucsuc --pred-block-size 1000 --coalescing-size 4000000 --suc-block-size 1000
#aprun -n 4 -N 1 -d 22 ./tc_family_prev --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 1 --preprocess --iterations 3 --scale 20 --run_tc_sucsuc --pred-block-size 1000 --coalescing-size 4000000 --suc-block-size 1000
#aprun -n 4 -N 1 -d 22 ./tc_family --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 1 --preprocess --iterations 3 --scale 20 --run_tc_blocked --pred-block-size 1000 --coalescing-size 4000000 --suc-block-size 10

#aprun -n 1 -N 1 -d 22 ./tc_family --file /lus/scratch/p02119/shared/tc/mm/ca-HepTh_adj.mmio --block-size 100,52,455 --read-graph --read-header --preprocess --no-weight --graph-type "mmio" --expected-triangles 28339 --threads 16 --id-distribution horizontal --warm-up-iterations 0 --iterations 3 --run_tc_striped --coalescing-size 40000 --verify
#aprun -n 1 -N 1 -d 22 ./tc_family --file /lus/scratch/p02119/shared/tc/mm/ca-HepTh_adj.mmio --read-graph --read-header --preprocess --no-weight --graph-type "mmio" --expected-triangles 28339 --threads 1 --id-distribution horizontal --warm-up-iterations 0 --iterations 3 --run_tc_level --coalescing-size 40000 --verify
#aprun -n 1 -N 1 -d 22 ./tc_family --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 1 --preprocess --iterations 3 --scale 18 --run_tc_blocked --coalescing-size 40000 --verify
#aprun -n 1 -N 1 -d 22 ./tc_family+pat --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 0 --preprocess --iterations 1 --scale 18 --run_tc_blocked --coalescing-size 40000
#aprun -n 1 -N 1 -d 22 ./tc_family --no-weight --threads 2 --id-distribution horizontal --warm-up-iterations 1 --preprocess --iterations 3 --scale 18 --run_tc_blocked --coalescing-size 40000
#aprun -n 1 -N 1 -d 22 ./tc_family --no-weight --threads 4 --id-distribution horizontal --warm-up-iterations 1 --preprocess --iterations 3 --scale 18 --run_tc_blocked --coalescing-size 40000
#aprun -n 1 -N 1 -d 22 ./tc_family --no-weight --threads 8 --id-distribution horizontal --warm-up-iterations 1 --preprocess --iterations 3 --scale 18 --run_tc_blocked --coalescing-size 40000
#aprun -n 1 -N 1 -d 22 ./tc_family --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 1 --preprocess --iterations 3 --scale 2 --run_tc_blocked --coalescing-size 40000
#aprun -n 4 -N 2 -d 16 ./tc_family --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 0 --preprocess --iterations 1 --scale 23 --run_tc_blocked --coalescing-size 40000 --suc-block-size 100 #--pred-block-size 100
#aprun -n 1 -N 1 -d 16 ./tc_family+pat --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 0 --preprocess --iterations 1 --scale 23 --run_tc_sucsuc --coalescing-size 3000000 --suc-block-size 250000 --pred-block-size 500000
#aprun -n 4 -N 2 -d 16 ./tc_family --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 0 --preprocess --iterations 1 --scale 23 --run_tc_blocked --coalescing-size 3000000 --suc-block-size 250634 --pred-block-size 3000
#aprun -n 2 -N 2 -d 22 ./tc_family --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 1 --preprocess --iterations 3 --scale 20 --run_tc_blocked --block-size 1000 --coalescing-size 440000 --suc-block-size 100
#aprun -n 1 -N 1 -d 22 ./tc_family+pat --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 0 --preprocess --iterations 1 --scale 23 --run_tc_sucsuc --pred-block-size 1000 --coalescing-size 440000 --suc-block-size 1000
#aprun -n 2 -N 2 -d 22 ./tc_family --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 0 --preprocess --iterations 1 --scale 23 --run_tc_sps --pred-block-size 10000 --coalescing-size 940000 --suc-block-size 100

#aprun -n 2 -N 2 -d 22 ./tc_family --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 1 --preprocess --iterations 3 --scale 23 --run_tc_blocked --pred-block-size 1000 --coalescing-size 440000 --suc-block-size 1000
#aprun -n 32 -N 2 -d 22 ./tc_family --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 0 --preprocess --iterations 1 --scale 29 --run_tc_blocked --block-size 10000 --coalescing-size 19200000 --suc-block-size 1000
#best configuration for scale 25
#aprun -n 32 -N 2 -d 22 ./tc_family --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 0 --preprocess --iterations 1 --scale 29 --run_tc_blocked --pred-block-size 10000 --coalescing-size 1920000 --suc-block-size 100
#aprun -n 2 -N 2 -d 22 ./tc_family --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 0 --preprocess --iterations 1 --scale 20 --run_tc_blocked --block-size 100 --coalescing-size 1100000 --suc-block-size 6500
#aprun -n 1 -N 1 -d 22 ./tc_family --no-weight --threads 1 --id-distribution horizontal --warm-up-iterations 0 --preprocess --iterations 1 --scale 21 --run_tc_blocked --pred-block-size 1000 --coalescing-size 1100000 --suc-block-size 100,1000,2000,3000
#aprun -n 8 -N 2 -d 22 ./tc_family --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 0 --preprocess --iterations 1 --scale 22 --seed 3432640012 --run_tc_blocked --block-size 100 --coalescing-size 1100000
#aprun -n 4 -N 2 -d 22 ./tc_family --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 0 --preprocess --iterations 1 --scale 22 --run_tc_blocked --block-size 10000 --coalescing-size 1100000 --flow-control 10 --suc-block-size 100
#aprun -n 4 -N 2 -d 22 ./tc_family --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 0 --preprocess --iterations 1 --scale 21 --run_tc_striped --block-size 10000 --coalescing-size 1100000 --flow-control 10

#aprun -n 1 -N 1 -d 44 ./tc_family --file /lus/scratch/p02119/graphs/twitter/out.twitter --read-graph --vertices 41652230 --edges 1468365182 --no-weight --graph-type "twitter" --threads 16 --id-distribution horizontal --warm-up-iterations 0 --iterations 3 --run_tc_blocked --coalescing-size 40000 --verify
#aprun -n 1 -N 1 -d 22 ./tc_family --file /lus/scratch/p02119/graphs/twitter/test.out --read-graph --preprocess --vertices 7  --edges 8 --no-weight --expected-triangles 5 --preprocess --graph-type "twitter" --threads 16 --id-distribution horizontal --warm-up-iterations 0 --iterations 1 --run_tc_level --coalescing-size 40000 --verify
#aprun -n 4 -N 2 -d 22 ./tc_family --file /lus/scratch/p02119/graphs/twitter/test.out --read-graph --preprocess --vertices 7  --edges 8 --no-weight --expected-triangles 5 --preprocess --graph-type "twitter" --threads 2 --id-distribution horizontal --warm-up-iterations 0 --iterations 3 --run_tc_blocked --coalescing-size 40000 --pred-block-size 1 --verify
#aprun -n 4 -N 2 -d 22 ./tc_family_no_ordering --file /lus/scratch/p02119/graphs/twitter/test.out --read-graph --preprocess --vertices 7  --edges 8 --no-weight --expected-triangles 5 --preprocess --graph-type "twitter" --threads 2 --id-distribution horizontal --warm-up-iterations 0 --iterations 1 --run_tc_blocked --coalescing-size 40000 --pred-block-size 1 --verify
#aprun -n 4 -N 2 -d 22 ./tc_family --file /lus/scratch/p02119/graphs/twitter/test.out --read-graph --preprocess --vertices 7  --edges 8 --no-weight --expected-triangles 5 --preprocess --graph-type "twitter" --threads 1 --id-distribution horizontal --warm-up-iterations 0 --iterations 1 --run_tc_blocked --coalescing-size 40000 --pred-block-size 1 --suc-block-size 1 --verify
#aprun -n 2 -N 2 -d 22 ./tc_family --file /lus/scratch/p02119/graphs/moreno_bison/out.moreno_bison_bison --read-graph --vertices 26 --edges 314 --preprocess --no-weight --graph-type "bison" --expected-triangles 1018 --threads 22 --id-distribution horizontal --warm-up-iterations 0 --iterations 1 --run_tc_blocked --block-size 100 --coalescing-size 40000 --suc-block-size 5 --verify
#aprun -n 1 -N 1 -d 1 ./tc_family --file /lus/scratch/p02119/graphs/moreno_bison/out.moreno_bison_bison --read-graph --vertices 26 --edges 314 --preprocess --no-weight --graph-type "bison" --expected-triangles 1018 --threads 44 --id-distribution horizontal --warm-up-iterations 0 --iterations 1 --run_tc_level --block-size 100000 --coalescing-size 40000 --verify
#aprun -n 1 -N 1 -d 22 ./tc_family --file /lus/scratch/p02119/graphs/citeseer/out.citeseer --read-graph --vertices 384413  --edges 1751463 --preprocess --no-weight --graph-type "citeseer" --expected-triangles 1351820 --threads 16 --id-distribution horizontal --warm-up-iterations 0 --iterations 3 --run_tc_blocked --coalescing-size 40000 --verify
#aprun -n 1 -N 1 -d 22 ./tc_family --file /lus/scratch/p02119/graphs/ca-HepTh.tsv --read-graph --vertices 9877  --edges 51971 --preprocess --no-weight --graph-type "ca-Hepth" --format zero --expected-triangles 28339 --threads 16 --id-distribution horizontal --warm-up-iterations 0 --iterations 3 --run_tc_striped --coalescing-size 40000 --verify
#aprun -n 2 -N 2 -d 22 ./tc_family --file /lus/scratch/p02119/graphs/ca-HepTh.tsv --read-graph --vertices 9877  --edges 51971 --preprocess --no-weight --graph-type "ca-Hepth" --format zero --expected-triangles 28339 --threads 16 --id-distribution horizontal --warm-up-iterations 0 --iterations 3 --run_tc_blocked --coalescing-size 40000 --verify
#aprun -n 1 -N 1 -d 22 ./tc_family+pat --file /lus/scratch/p02119/graphs/soc-LiveJournal1/out.soc-LiveJournal1 --read-graph --preprocess --vertices 4847571  --edges 68475391 --no-weight --graph-type "soclive" --expected-triangles 285730264 --threads 16 --id-distribution horizontal --warm-up-iterations 0 --iterations 3 --run_tc_blocked --coalescing-size 240000 --block-size 10000 --verify
#aprun -n 2 -N 2 -d 22 ./tc_family --file /lus/scratch/p02119/graphs/soc-LiveJournal1/out.soc-LiveJournal1 --read-graph --preprocess --vertices 4847571  --edges 68475391 --no-weight --graph-type "soclive" --expected-triangles 285730264 --threads 16 --id-distribution horizontal --warm-up-iterations 0 --iterations 1 --run_tc_sucsuc --coalescing-size 240000 --pred-block-size 100000 --suc-block-size 300 --verify
#aprun -n 1 -N 1 -d 22 ./tc_family --file /lus/scratch/p02119/graphs/soc-LiveJournal1/out.soc-LiveJournal1 --read-graph --vertices 4847571  --edges 68475391 --no-weight --graph-type "soclive" --expected-triangles 285730264 --threads 16 --id-distribution horizontal --warm-up-iterations 0 --iterations 1 --run_tc_sps --coalescing-size 240000 --pred-block-size 100000 --suc-block-size 300

#aprun -n 4 -N 2 -d 22 ./tc_family --file /lus/scratch/p02119/graphs/soc-LiveJournal1/out.soc-LiveJournal1 --read-graph --preprocess --vertices 4847571  --edges 68475391 --no-weight --graph-type "soclive" --expected-triangles 285730264 --threads 16 --id-distribution horizontal --warm-up-iterations 0 --iterations 3 --run_tc_sucsuc --coalescing-size 240000 --pred-block-size 10000 --suc-block-size 1000 --verify


#aprun -n 4 -N 2 -d 22 ./tc_family --file /lus/scratch/p02119/graphs/friendstr/friendster/out.friendster --read-graph --preprocess --vertices 4847571  --edges 68475391 --no-weight --graph-type "friendstr" --expected-triangles 285730264 --threads 16 --id-distribution horizontal --warm-up-iterations 0 --iterations 1 --run_tc_blocked --coalescing-size 2400000 --pred-block-size 1000 --suc-block-size 1000 --verify

#aprun -n 4 -N 2 -d 22 ./tc_family --file /lus/scratch/p02119/graphs/soc-LiveJournal1/out.soc-LiveJournal1 --read-graph --preprocess --vertices 4847571  --edges 68475391 --no-weight --graph-type "soclive" --expected-triangles 285730264 --threads 16 --id-distribution horizontal --warm-up-iterations 0 --iterations 1 --run_tc_blocked --coalescing-size 2400000 --pred-block-size 10,100,1000,10000,100000 --suc-block-size 80 --verify
#aprun -n 4 -N 2 -d 22 ./tc_family --file /lus/scratch/p02119/graphs/soc-LiveJournal1/out.soc-LiveJournal1 --read-graph --preprocess --vertices 4847571  --edges 68475391 --no-weight --graph-type "soclive" --expected-triangles 285730264 --threads 16 --id-distribution horizontal --warm-up-iterations 0 --iterations 1 --run_tc_sucsuc --coalescing-size 2400000 --pred-block-size 10,100,1000,10000,100000 --suc-block-size 1000 --verify
#aprun -n 4 -N 2 -d 22 ./tc_family_prev --file /lus/scratch/p02119/graphs/soc-LiveJournal1/out.soc-LiveJournal1 --read-graph --preprocess --vertices 4847571  --edges 68475391 --no-weight --graph-type "soclive" --expected-triangles 285730264 --threads 16 --id-distribution horizontal --warm-up-iterations 0 --iterations 1 --run_tc_sucsuc --coalescing-size 2400000 --pred-block-size 10,100,1000,10000,100000 --suc-block-size 80 --verify



#aprun -n 1 -N 1 -d 22 ./tc_family --file /lus/scratch/p02119/graphs/soc-LiveJournal1/out.soc-LiveJournal1 --read-graph --vertices 4847571 --edges 68475391 --preprocess --no-weight --graph-type "soclive" --expected-triangles 285730264 --threads 16 --id-distribution horizontal --warm-up-iterations 0 --iterations 1 --run_tc_level --block-size 100000 --coalescing-size 40000 --verify
#aprun -n 1 -N 1 -d 22 ./tc_family --file /lus/scratch/p02119/graphs/soc-LiveJournal1/out.soc-LiveJournal1 --read-graph --vertices 4847571 --edges 68475391 --preprocess --no-weight --graph-type "soclive" --expected-triangles 285730264 --threads 16 --id-distribution horizontal --warm-up-iterations 0 --iterations 1 --run_tc_blocked --coalescing-size 40000 --verify
#aprun -n 1 -N 1 -d 22 ./mis_family --threads 16 --warm-up-iterations 1 --iterations 3 --scale 20 --luby_algorithms AV2 --coalescing-size 20000
#aprun -n 1 -N 1 -d 22 ./mis_family --threads 16 --id-distribution vertical --warm-up-iterations 0 --iterations 1 --scale 20 --luby_algorithms A --coalescing-size 40000 --verify
#aprun -n 1 -N 1 -d 22 ./mis_family --threads 16 --id-distribution vertical --warm-up-iterations 0 --iterations 1 --scale 20 --luby_algorithms AV1 --coalescing-size 40000 --verify
#aprun -n 1 -N 1 -d 22 ./mis_family --threads 16 --id-distribution vertical --warm-up-iterations 0 --iterations 1 --scale 20 --luby_algorithms AV2 --coalescing-size 40000 --verify
#aprun -n 1 -N 1 -d 22 ./mis_family --threads 16 --id-distribution vertical --warm-up-iterations 0 --iterations 1 --scale 20 --luby_algorithms B --coalescing-size 40000 --verify
#aprun -n 1 -N 1 -d 22 ./mis_family --threads 2 --id-distribution horizontal --warm-up-iterations 0 --iterations 1 --scale 20 --luby_algorithms AV1 --coalescing-size 40000 
#aprun -n 2 -N 1 -d 22 ./mis_family --threads 16 --id-distribution horizontal --warm-up-iterations 0 --iterations 1 --scale 25 --luby_algorithms AV1 --coalescing-size 40000 
#aprun -n 2 -N 1 -d 22 ./mis_family --threads 16 --id-distribution horizontal --warm-up-iterations 0 --iterations 1 --scale 25 --luby_algorithms AV2 --coalescing-size 40000 
#aprun -n 1 -N 1 -d 22 ./mis_family --threads 2 --id-distribution horizontal --warm-up-iterations 0 --iterations 1 --scale 20 --luby_algorithms B --coalescing-size 40000 
#aprun -n 1 -N 1 -d 22 ./mis_family --threads 16 --warm-up-iterations 1 --iterations 3 --scale 20 --luby_algorithms A --coalescing-size 20000
#aprun -n 2 -N 1 -d 22 ./mis_family --threads 16 --warm-up-iterations 1 --iterations 3 --scale 20 --luby_algorithms B --coalescing-size 20000
