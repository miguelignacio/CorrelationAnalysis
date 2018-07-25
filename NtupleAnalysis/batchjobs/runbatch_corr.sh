#!/bin/bash

if [[ $# -eq 1 ]] ; then
    echo 'please give [13d, 13e, 17q,...] and [MinBias.hdf5] [test/full] as arguments'
    exit 0
fi

date

module add ROOT/6.08.00

#for p in {0,4,6}
for p in 0
do
    name=$(echo "$filename" | cut -f 1 -d '.')
    if [ ! -f $name_${p}GeVTrack_paired.root ]; then
	echo "calling paired Injector with file $1.root and track Gev $p"
	./../pair_gale_shapley/paired_injector $1.root $p
    fi
    
    for i in {0..280..20} #Mix 300 events
    do
	mix_min=$i
	mix_max="$((i + 19))"
	if [[ $3 == full ]]; then
	    sbatch -p shared-chos -t 16:00:00 runCorr.sh $1_${p}GeVTrack_paired.root $2 $mix_min $mix_max $p	  
	else
	./runCorr.sh $1_${p}GeVTrack_paired.root $2 $mix_min $mix_max $p
	fi
    echo "$mix_min $mix_max $1"
    done
done