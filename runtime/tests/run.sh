
#export MPICH_NEMESIS_ASYNC_PROGRESS=SC
export MPICH_MAX_THREAD_SAFETY=multiple
#export MPICH_GNI_USE_UNASSIGNED_CPUS=enabled
#export OMP_NUM_THREADS=16


aprun -n 2 -N 2 -d 16 -e MPICH_MAX_THREAD_SAFETY=multiple ./varbuffer 16
