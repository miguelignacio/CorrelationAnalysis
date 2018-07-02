#!/bin/bash

if [[ $# -eq 0 ]] ; then
    echo 'Please Give Dataset Argument [13d] [13e] [13f]'
    exit 0
fi

date

#for p in {0,4,6}
for p in 0
do
    for i in {0..280..20} #Mix 300 events
    do
	mix_min=$i
	mix_max="$((i + 19))"
	sbatch -p shared-chos -t 20:00:00 runpairing.sh $1 $mix_min $mix_max $p
	#./runpairing.sh $1 $mix_min $mix_max $p
    done
done