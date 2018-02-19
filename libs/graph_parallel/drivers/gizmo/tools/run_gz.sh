export CORES=4
export SCALE=23


#aprun -n 32 -N 2 -d 1 ./calc --iterations 10000 --format
aprun -n 1 -N 1 -d 1 ./calc --iterations 10000 --format
#aprun -n 16 -N 2 -d 1 ./calc_old --iterations 10000

