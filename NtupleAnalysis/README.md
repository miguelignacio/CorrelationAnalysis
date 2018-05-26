# Gamma Hadron Correlation Code

This code is a mix of c++ and python code to extract a gamma hadron, and ultimatley fragmentation function involving isolated, single photons and associated charged hadrons from a hard scattering. This uses NTuples created in [Yue Shi Lai's NTupilizer](https://github.com/yslai/ntuple-gj). hdf5 and Gale Shapeley pairing also built from examples

## Getting Started

I will try to keep the repository updated with root files to at least be able to run the notebook. However, because the root files needed for running the c++ code (Mix_Correlatons.cc and Frixione_Correlations.cc) are larger than githubs single file limit, the .root extension has been added to the .git ignore. Essentially, simply cloning this repository should be sufficient to run the notebook, while running the c code will require additional NTuples.

### Prerequisites

For creating the `.root` files
```
root6
hdf5
```

For runnning the jupyter notebook
```
jupyter notebook
rootpy
pandas
matplotlib
AstroML
```

### Running

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
