#!/bin/sh
MPIRUN=$HOME/opt/usr/local/bin/mpirun
$MPIRUN -np $1 --hostfile ./hostfile ./testForest $2
