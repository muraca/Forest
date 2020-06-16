#!/bin/sh
mpirun -np $1 --hostfile ./hostfile ./testForest $2
