#!/bin/bash

source ~/.bashrc.ext

module add ROOT/6.08.00
module load hdf5
which gcc
which root
module list
date
echo "Command runpairing with mix_start = $1, mix_end = $2, and TrackSkim GeV = $3"

export cwd=$(pwd)
cd $cwd/../pair_gale_shapley/
./mix_gale_shapley ../InputData/13def.root ../InputData/13c_pass4_v1_$3GevtrackSkim.root $1 $2 $3
date 