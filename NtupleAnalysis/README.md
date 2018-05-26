# Gamma Hadron Correlation Code

This code is a mix of c++ and python code to extract a gamma hadron, and ultimatley fragmentation function involving isolated, single photons and associated charged hadrons from a hard scattering. This uses NTuples created in [Yue Shi Lai's NTupilizer](https://github.com/yslai/ntuple-gj). Gale Shapeley and hdf5 pairing also built from examples. This is meant to run and PDSF.

## Getting Started

### Files

Root files and hdf5 files should be placed in the `InputData` folder. There is a set of hdf5 and 13c root files for each track skimmed pT threshold (0, 4, and 6 GeV). These files are:

```
/project/projectdirs/alice/NTuples/pPb/13def_v1_small.root
/project/projectdirs/alice/NTuples/pPb/13def_v2_small.root

/project/projectdirs/alice/NTuples/pPb/13c/13c_pass4_v1_0GevtrackSkim.root
/project/projectdirs/alice/NTuples/pPb/13c/13c_pass4_v1_4GevtrackSkim.root
/project/projectdirs/alice/NTuples/pPb/13c/13c_pass4_v1_6GevtrackSkim.root

/project/projectdirs/alice/NTuples/hdf5/13c_pass4_v1_0GeV_Skim.hdf5
/project/projectdirs/alice/NTuples/hdf5/13c_pass4_v1_4GeV_Skim.hdf5
/project/projectdirs/alice/NTuples/hdf5/13c_pass4_v1_6GeV_Skim.hdf5
```

### Prerequisites

#### PDSF Modules
module add ROOT/6.08.00
module load hdf5

#### python requirements for notebook
```
jupyter notebook
rootpy
pandas
matplotlib
AstroML
```

## Running

### Batch Jobs
Paring must be done first, then correlation. 

At the moment, the pairing replicates the small 13def Ntuple for each mixing range into the InputData folder with an appropriate name. An exit script to merge this was not trivial, as jobs do not finish sequentially. This current implementaiton at least does not require any additional commands between running the pairing and the mixed correlation function.

#### Batch Paring
Currently hardcoded to run for the v1 clusterizer 13def gamma triggered data set.

Make mix_gale_shapley.cc in the `pair_gale_shapley/` folder

All batch scripts take a track skimming pT as an argument (0, 4 or 6 GeV)
In the `batchjobs/` folder do `./runbatch_pair.sh [Track GeV]`

#### Batch Correlation

Make Batch_Mix_Correlations.cc in the `NtupleAnalysis/` folder
./runbatch_corr.sh

#### Testing before submitting
can run individual test on the pairing or the correlation with 

runCorr.sh and runpairing.sh

however, the mixing range must be given. A very quick way to test would be:
`./runCorr.sh 0 1 6`

This runs for one mixed event, and uses 6GeV track skimmed 13c Ntuple.

### Combining Correlation function
simply `hadd` the various mixed correlation function outputs.
The final output *needs* to be called `Mix_Correlation_[0,4,6]GeVTracks_13defv1_ALL.root` for the last step to work

### Same_Mix_Ratio
Divides the same event correlation by the mixe event correlation in the appropriate zT bin, taking from all skimmed correlation functions.

Make Same_Mix_Ratio
./Same_Mix_Ratio

## Disclaimer
These have been tested extensiveley, however files were recently moved around in order to consolidate repositories, and make my and other users' life much easier when running. Any errors are likeley related to files being in the wrong spot. Apologies in advance.

#### Notebook
Simlpy start the notebook (jupyter-notebook at a terminal) and you should be able to run over different NTuples by removing the comment prefix for the appropriate file. Save images is a pics/ folder.

#### Correlation Code
Make sure the gamma triggered NTuple has had the mix_GS run on it from the [NTupilizer](https://github.com/yslai/ntuple-gj). This populates a branch with a list of mixed events to run over.

Command:
/Mix_Correlations [gamma triggered NTuple] [Min Bias NTuple]
/Frixione_Correlations [gamma triggered NTuple]

#### Divide Correlations
The naming of root files here is very strict. The same event correlation function will be called something along the lines of `Same_Event_Correlation.root`, however the last part. The Mixed event correlation functios will be multiple root files with the name `Mix_Correlation_[i]GeVTracks_13def.root`
where `i = {0,4,6,30}` at the moment. These contain the mixed event correlation functions that have been split up to address the lack of statistics in high zT bins.

The names of these files is hardcoded in the Correlation code. A check that the names of the root files created in the Correlation code should match what the divisios macro wants should be done regularly. But the name can alway be changed after the Correlation root file is updated.

The correlation code then picks which bin each file will contribute to, and normalizes the same event correlation function in the appropriate zT bin.

### More to come. Would like to create a doxygen for this as well. 
