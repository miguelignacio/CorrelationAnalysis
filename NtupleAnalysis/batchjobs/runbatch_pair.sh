#!/bin/bash

if [[ $# -eq 0 ]] ; then
    echo 'please give Track Skimming GeV as argument'
    exit 0
fi

date
for i in {0..280..20} #Mix 300 events
do
    mix_min=$i
    mix_max="$((i + 19))"
    sbatch -p shared-chos -t 20:00:00 runpairing.sh $mix_min $mix_max $1
done
