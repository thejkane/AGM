#!/bin/bash

export LD_LIBRARY_PATH=/N/u/thejkane/BigRed2/crest/ampp/installation/ampp/lib:$LD_LIBRARY_PATH

aprun -b -n 10 -N 1 -d 32 ./europar_tests --delta-stepping --threads 32 --scale 15 --num-sources 10 --coalescing-size 3000 --max-weight 100 --delta 2 --with-per-thread-reductions --verify
