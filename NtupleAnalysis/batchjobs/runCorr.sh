#!/bin/bash

source ~/.bashrc.ext
module add ROOT/6.08.00
module load hdf5
which gcc
which root
module list
date
echo "Command runCorr with mix_start = $1, mix_end = $2, and TrackSkim GeV = $3"

export cwd=$(pwd)
cd $cwd/../
./Batch_Mix_Correlations InputData/13defv1_c_$1_$2_$3GeV_TrackSkim_mixed.root InputData/13c_pass4_v1_$3GeV_Skim.hdf5 $1 $2 $3
date
