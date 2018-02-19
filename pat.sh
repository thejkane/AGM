#!/bin/sh

export BINARY=mis_family

module load perftools-base
module load perftools

cd ../pat
rm *.txt

export PAT_RT_PERFCTR=PAPI_CA_INV,PAPI_CA_CLN,PAPI_CA_INV,PAPI_CA_ITV

cd ../build/bin
rm $BINARY*

cd ../
./b_rel.sh

cd bin
pat_build -g pthreads -g heap -w $BINARY
#pat_build a.out

#run
aprun -n 1 -d 22 ./a.out+pat

#pat_report -s ap2=no,th=ALL -O thread_times -o output.txt *.xf 
#pat_report -s ap2=no,th=ALL -o output.txt *.xf -- good 1
# -O profile = Only shows the function group table (1st table in default pat_report)
#TODO : Experiment with MPI
pat_report -s ap2=no,pe=ALL,th=ALL,pe.th=ALL,aggr_th=sum,regions=show -b gr,fu=ALL,th -d counters -O calltree,ct+src,ca+src,load_balance_function,load_imbalance_thread,program_time -o output.txt *.xf 
mv output.txt ../pat/1.txt
