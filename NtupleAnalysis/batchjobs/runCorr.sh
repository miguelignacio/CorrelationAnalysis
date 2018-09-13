#!/bin/bash

source ~/.bashrc.ext
module add ROOT/6.08.00
module load hdf5
which gcc
which root
module list
date
echo "Command runCorr with file = $1, hdf5 = $2 mix_start = $3, mix_end = $4, and TrackSkim GeV = $5"

#export cwd=$(pwd)
#cd $cwd/../
./Batch_Mix_Correlations $1 $2 $3 $4 $5
#./MC_Batch_Mix_Correlations $1 $2 $3 $4 $5
date
