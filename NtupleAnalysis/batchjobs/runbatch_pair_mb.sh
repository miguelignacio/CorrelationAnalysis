#!/bin/bash

if [[ $# -eq 0 ]] ; then
    echo 'Please Give Datafile Argument [13d] [13e] [13f]'
    exit 0
fi

date


for i in {0..280..20} #Mix 300 events
do
    mix_min=$i
    mix_max="$((i + 19))"
    sbatch -p shared-chos -t 30:00:00 run_mb_pairing.sh $1 $mix_min $mix_max 4
    #./run_mb_pairing.sh $1 $mix_min $mix_max 4
    echo $1 $mix_min $mix_max 4
done
