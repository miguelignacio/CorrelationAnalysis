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
./Batch_Mix_Correlations InputData/$1_$4GeVTrack_paired.root InputData/$1_minbias_$4GeVTracks.hdf5 $2 $3 $4
date
