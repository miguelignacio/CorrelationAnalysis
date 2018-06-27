#!/bin/bash

source ~/.bashrc.ext
module add ROOT/6.08.00
module load hdf5
which gcc
which root
module list
date
echo "Command runCorr with file = $1, mix_start = $2, mix_end = $3, and TrackSkim GeV = $4"

export cwd=$(pwd)
cd $cwd/../
#./Batch_Mix_Correlations InputData/$1_$4GeVTrack_paired.root InputData/$1_$4GeVTrack.hdf5 $1 $2 $3
./Batch_Mix_Correlations InputData/$1_$4GeVTrack_paired.root InputData/$1.hdf5 $2 $3 $4
date

#First file name depends on output of Ivan's script. Tentatively will call it 13$1_$4GeV_MinBias_paired.root
#13c_pass4_v1_4GevtrackSkim_paired.root
#"%s_%iGeVTrack_paired.root
./Batch_Mix_Correlations InputData/13c_pass4_v1_4GevtrackSkim_4GeVTrack_paired.root InputData/13c_pass4_v1_4GevtrackSkim.hdf5 0 19 4
