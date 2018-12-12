#!/bin/bash

if [[ $# -eq 1 ]] ; then
    echo 'please give [13d, 13e, 17q,...] and [MinBias.hdf5] [test/full] as arguments'
    exit 0
fi

date

module add ROOT/6.08.00

#for p in {0,4,6}
for p in 4
do
    echo $1
    name=${1%.*}
    #name=$(basename ../$1 .root)
    if [ ! -f ${name}_${p}GeVTrack_paired.root ]; then
	echo "calling paired Injector with file $name.root and track Gev $p"
	./../pair_gale_shapley/paired_injector $name.root $p
    fi
    
    for i in {0..280..20} #Mix 300 events
    do
	mix_min=$i
	mix_max="$((i + 19))"
	if [[ $3 == full ]]; then
	    sbatch -p shared-chos -t 16:00:00 runCorr.sh ${name}_${p}GeVTrack_paired.root $2 $mix_min $mix_max $p	  
	else
	./runCorr.sh ${name}_${p}GeVTrack_paired.root $2 $mix_min $mix_max $p
	fi
    echo "$mix_min $mix_max $name"
    done
done