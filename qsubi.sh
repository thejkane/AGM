#!/bin/bash
#qsub interactive
qsub -I -l walltime=01:00:00 -l nodes=10:ppn=32 -q debug_cpu -k o -m abe -N pbgldebug -j oe
