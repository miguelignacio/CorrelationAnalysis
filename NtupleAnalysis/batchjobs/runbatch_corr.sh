#!/bin/bash

if [[ $# -eq 0 ]] ; then
    echo 'please give Track Skimming GeV as argument'
    exit 0
fi

date

#FIXME: call Ivan's script here

for p in {0,4,6}
do
    for i in {0..280..20} #Mix 300 events
    do
	mix_min=$i
	mix_max="$((i + 19))"
	sbatch -p shared-chos -t 12:00:00 runCorr.sh $1 $mix_min $mix_max $p
    #echo "$mix_min $mix_max $1"
    done
done