#!/bin/bash

source ~/.bashrc.ext

module add ROOT/6.08.00
module load hdf5
which gcc
which root
module list
date
echo "Command runpairing with [$1] [$2] [mix_start = $3] [mix_end = $4] [TrackSkim = $5]"

export cwd=$(pwd)
cd $cwd/../pair_gale_shapley/
#./mix_gale_shapley ../InputData/$1.root ../InputData/$1_minbias_$4GeVTracks.root $2 $3 $4
./mix_gale_shapley $1 $2 $3 $4 $5
date 


#The first argument can be replaced with full file name. Written assuming use on 13d,13e, and 13f (separately)