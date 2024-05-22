#!/bin/bash
#$ -V
#$ -j yes
#$ -cwd
    #
    # This script runs D-Emissions on Linux
    #
scriptdirname=`readlink \-f \$0`
scriptdir=`dirname $scriptdirname`
D3D_HOME=$scriptdir/..
bindir=$D3D_HOME/bin
$bindir/run_delwaq.sh $1 $2 $3 $4 $5 $6 $7 $8 $9 -dem
