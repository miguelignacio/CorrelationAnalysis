#!/bin/bash

if [[ $# -eq 0 ]] ; then
    echo 'please give Track Skimming GeV as argument'
    exit 0
fi

#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -q regular
#SBATCH --mail-user=fernando_tta@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH -t 10:30:00

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread


#run the application:
date
for i in {0..40..20} #Mix 300 events
do
    mix_min=$i
    mix_max="$((i + 19))"
    #SBATCH -J 6GeV_$mix_min_$mix_max_Correlations
    srun -n 1 -c 64 --cpu_bind=cores runCorr.sh $mix_min $mix_max $1
done
