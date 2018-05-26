#!/bin/bash

source ~/.bashrc.ext
#module load slurm
module list
module add ROOT/6.08.00
which gcc
date
cd /project/projectdirs/alice/ftorales/CorrelationAnalysis/NtupleAnalysis/pair_gale_shapley/
./mix_gale_shapley ../InputData/13def_v1_small.root ../InputData/13c_pass4_v1_$3GevtrackSkim.root $1 $2 $3
date 
#echo $1 $2
