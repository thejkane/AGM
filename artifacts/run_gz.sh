cd bin

export LD_LIBRARY_PATH=/home/users/p02119/development/cds-2.1.0/bin/cc-amd64-linux-64/:$LD_LIBRARY_PATH


export CORES=4
export SCALE=23


aprun -n 1 -N 1 -d 16 ./gizmo --no-weight --threads 16 --id-distribution horizontal --warm-up-iterations 0 --rmata 2500 --rmatbc 2500 --iterations 1 --scale 10 --run_tc_sucsuc --pred-block-size 1000 --coalescing-size 4000000 --suc-block-size 1000 --gizmo-file /home/users/p02119/development/pbgl2/libs/graph_parallel/drivers/gizmo/config/gizmo.cfg

