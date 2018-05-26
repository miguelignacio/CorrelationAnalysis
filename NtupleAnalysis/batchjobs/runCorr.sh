#!/bin/bash

source ~/.bashrc.ext
#module load slurm
module list
module add ROOT/6.08.00
module load hdf5
which gcc
date
cd /project/projectdirs/alice/ftorales/CorrelationAnalysis/NtupleAnalysis/
./Batch_Mix_Correlations InputData/13defv1_c_$1_$2_$3GeV_TrackSkim_mixed.root InputData/13c_pass4_v1_$3GeV_Skim.hdf5 $1 $2 $3
#echo "$1 $2 Track GeV Skimmed: $3GeV"
date
#echo $1 $2
