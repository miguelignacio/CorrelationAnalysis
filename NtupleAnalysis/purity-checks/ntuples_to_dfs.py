import numpy as np
import pandas as pd

try:
    import root_pandas as rpd
except ImportError:
    import sys
    sys.path.insert(0, '/usr/common/software/rootpy/')
    import root_pandas as rpd

import os
import time


def getWeights(filename, dataframe):
    if '17g6a3' in filename:
        if 'pthat1' in filename:
            dataframe.eval('weights=4.47e-11', inplace=True)
        elif 'pthat2' in filename:
            dataframe.eval('weights=9.83e-11', inplace=True)
        elif 'pthat3' in filename:
            dataframe.eval('weights=1.04e-10', inplace=True)
        elif 'pthat4' in filename:
            dataframe.eval('weights=1.01e-10', inplace=True)
        elif 'pthat5' in filename:
            dataframe.eval('weights=6.93e-11', inplace=True)
    elif '17g6a1' in filename:
        if('pthat1' in filename):
            dataframe.eval('weights=1.60e-11', inplace=True)
        elif('pthat2' in filename):
            dataframe.eval('weights=2.72e-12', inplace=True)
        elif('pthat3' in filename):
            dataframe.eval('weights=3.96e-13', inplace=True)
        elif('pthat4' in filename):
            dataframe.eval('weights=6.14e-14', inplace=True)
        elif('pthat5' in filename):
            dataframe.eval('weights=1.27e-14', inplace=True)
    elif 'pthat' in filename:
        avg_eg_ntrial = str(np.mean(dataframe['eg_ntrial']))
        dataframe.eval('weights=eg_cross_section/'+avg_eg_ntrial, inplace=True)
    else:
        dataframe.eval('weights=1', inplace=True)


def getData(inputFiles, basedir='ntuples', maxEvents=None):
    scalarColumns = []
    scalarColumns.append('eg_cross_section')
    scalarColumns.append('eg_ntrial')
    scalarColumns.append('ue_estimate_its_const')

    arrayColumns = []
    arrayColumns.append('cluster_pt')
    arrayColumns.append('cluster_eta')
    arrayColumns.append('cluster_ncell')
    arrayColumns.append('cluster_e_cross')
    arrayColumns.append('cluster_e')
    arrayColumns.append('cluster_e_max')
    arrayColumns.append('cluster_tof')
    arrayColumns.append('cluster_nlocal_maxima')
    arrayColumns.append('cluster_distance_to_bad_channel')
    arrayColumns.append('cluster_iso_its_04')
    arrayColumns.append('cluster_iso_its_04_ue')
    arrayColumns.append('cluster_NN1')
    arrayColumns.append('cluster_Lambda')
    arrayColumns.append('cluster_iso_04_truth')

    columns = scalarColumns + arrayColumns

    dfs = []
    for inputFile in inputFiles:
        start = time.time()
        filename = os.path.basename(inputFile)
        dataframe = rpd.read_root(os.path.join(basedir, inputFile), columns=columns, flatten=arrayColumns, stop=maxEvents)
        getWeights(filename, dataframe)
        end = time.time()
        print 'Processed {0} in {1} seconds'.format(filename, end - start)
        dfs.append(dataframe)

    return pd.concat(dfs).drop_duplicates().reset_index(drop=True)


def applyCut(inputDataframe, cut, text=None):
    cutDataframe = inputDataframe.query(cut)
    if text:
        print '{0}: {1}'.format(text, cutDataframe.shape[0])
    return cutDataframe


def applyCuts(fullDataframe, ptRange=(15.0, 20.0)):
    fullDataframe.eval('cluster_ecross_over_e = cluster_e_cross/cluster_e', inplace=True)
    fullDataframe.eval('cluster_emax_over_e = cluster_e_max/cluster_e', inplace=True)
    # fullDataframe.eval('cluster_iso_its_04_raw = cluster_iso_its_04 + cluster_iso_its_04_ue', inplace=True)
    fullDataframe.eval('cluster_iso_its_04_sub = cluster_iso_its_04 + cluster_iso_its_04_ue - ue_estimate_its_const*0.4*0.4*3.1416', inplace=True)
    fullDataframe.eval('cluster_iso_04_truth_sub = cluster_iso_04_truth - cluster_pt', inplace=True)

    dataframe = fullDataframe
    dataframe = applyCut(dataframe, 'cluster_pt>{0} and cluster_pt<{1}'.format(*ptRange), '{0} < pt < {1}'.format(*ptRange))
    dataframe = applyCut(dataframe, 'abs(cluster_eta)<0.67', '|eta| < 0.67')
    dataframe = applyCut(dataframe, 'cluster_ncell>2', 'ncell > 2')
    dataframe = applyCut(dataframe, 'cluster_ecross_over_e>0.05', 'ecross/e > 0.05')
    dataframe = applyCut(dataframe, 'abs(cluster_tof)<20', '|tof| < 20')
    dataframe = applyCut(dataframe, 'cluster_nlocal_maxima<3', 'N local maxima < 3')
    dataframe = applyCut(dataframe, 'cluster_distance_to_bad_channel>=2.0', 'distance to bad channel >= 2.0')
    dataframe = applyCut(dataframe, 'cluster_NN1 == cluster_NN1', 'cluster_NN1 is not NaN')

    return dataframe


def ptcuttext(ptrange):
    return '{0}>{1} and {0}<{2}'.format('cluster_pt', *ptrange)
