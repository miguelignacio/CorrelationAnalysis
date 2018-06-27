#!/bin/bash

if [[ $# -eq 0 ]] ; then
    echo 'Please Give Aguments [13d, 13e, 13f]'
    exit 0
fi

date

echo "USing hardcoded track GeV of 0"
for i in {0..280..20} #Mix 300 events
do
    mix_min=$i
    mix_max="$((i + 19))"
    sbatch -p shared-chos -t 30:00:00 run_mb_pairing.sh $1 $mix_min $mix_max 0
    #./run_mb_pairing.sh $1 $mix_min $mix_max 0
    echo $1 $mix_min $mix_max 4
done
