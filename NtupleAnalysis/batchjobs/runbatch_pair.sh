#!/bin/bash

if [[ $# -eq 1 ]] ; then
    echo 'Please Give Dataset Argument [triggered root] [MB root] [full/test]'
    exit 0
fi

date

#for p in {0,4,6}
for p in 4
do
    for i in {0..280..20} #Mix 300 events
    do
	mix_min=$i
	mix_max="$((i + 19))"

        if [[ $3 == full ]]; then
	    sbatch -p shared-chos -t 24:00:00 runpairing.sh $1 $2 $mix_min $mix_max $p
        else
	./runpairing.sh $1 $2 $mix_min $mix_max $p
        fi

    done
done