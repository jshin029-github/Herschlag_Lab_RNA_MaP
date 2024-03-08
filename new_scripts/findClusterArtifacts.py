import pandas as pd
import numpy as np
import os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import skimage

# Usage: python makeClusterHeatmap.py [CPseries File with Path] [FID fluor Cutoff] [Background fluor Cutoff]

# Set Resolution within which to average (e.g. setting to 400 => avg within 400 X 400 blocks)
Resolution = 400
##################################################

Filename = sys.argv[1]
FidCutoff = float(sys.argv[2])
BgCutoff = float(sys.argv[2])

def get_mask(df, colname, clustype):
    if clustype == "FID" or clustype == "anyRNA":
        mask = (df[colname] == clustype)
    else:
        mask = (df[colname] != "FID") & (df[colname] != "anyRNA")
    return mask

fig,axs = plt.subplots(nrows=len(folders), ncols=3, sharex =True, sharey=True, figsize=(20, 35))

# Read in current tile data:
print "Reading in file: %s" % (Filename)
df = pd.read_csv(Filename, sep="\t")
df.columns = ['clusterID', 'clusterType', 'C2', 'C3', 'C4', 'C5', 'fluorescence']

## WORKING ON BACKGROUND CLUSTERS:
dfbg = df.loc[get_mask(df, 'clusterType', 'background')]
dfbg = dfbg.loc[:, ['clusterID', 'fluorescence']]
dfbg = dfbg.dropna()
# Parse clusters in current tile:
clusterIDs = dfbg['clusterID']
x = []
y = []
for clusterID in clusterIDs:
    words = clusterID.strip().split(':')
    X = int(words[5])
    Y = int(words[6])
    x.append(X)
    y.append(Y)
dfbg['X'] = x
dfbg['Y'] = y

# WORKING ON RNA CLUSTERS:
dfrna = df.loc[get_mask(df, 'clusterType', 'anyRNA')]
dfrna = dfrna.loc[:, ['clusterID', 'fluorescence']]
dfrna = dfrna.dropna()
# Parse clusters in current tile:
clusterIDs = dfrna['clusterID']
x = []
y = []
for clusterID in clusterIDs:
    words = clusterID.strip().split(':')
    X = int(words[5])
    Y = int(words[6])
    x.append(X)
    y.append(Y)
dfrna['X'] = x
dfrna['Y'] = y

# WORKING ON FID CLUSTERS:
dffid = df.loc[get_mask(df, 'clusterType', 'FID')]
dffid = dffid.loc[:, ['clusterID', 'fluorescence']]
dffid = dffid.dropna()
# Parse clusters in current tile:
clusterIDs = dffid['clusterID']
x = []
y = []
for clusterID in clusterIDs:
    words = clusterID.strip().split(':')
    X = int(words[5])
    Y = int(words[6])
    x.append(X)
    y.append(Y)
dffid['X'] = x
dffid['Y'] = y

print "Background clusters: %d", dfbg.shape[0]
print "RNA clusters: %d", dfrna.shape[0]
print "Fiducial clusters: %d", dffid.shape[0]

Xmax = max(dfbg['X'].max(), dfrna['X'].max(), dffid['X'].max())
Ymax = max(dfbg['Y'].max(), dfrna['Y'].max(), dffid['Y'].max())
N = max(Xmax, Ymax) + Resolution - max(Xmax, Ymax)%Resolution
print "Using matrix size: %d"%(N)

clusterTypes = ['Background', 'RNA', 'Fiducial']
dfs = [dfbg, dfrna, dffid]
badClusters = []
#fig,axs = plt.subplots(nrows=1, ncols=3, sharex =True, sharey=True, figsize=(20, 6))
i = 0
for df, clustype in zip(dfs, clusterTypes):
    # Create matrix from dataset with X and Y coordinates:
    M = np.empty((N,N,),)
    M[:] = np.nan
    for index, row in df.iterrows():
        M[row['X'], row['Y']] = row['fluorescence'] 

    print np.nanmean(M)
    # Block Reduce and mean this matrix:
    MReduced = skimage.measure.block_reduce(M, (Resolution, Resolution), np.nanmedian)
    

    print clustype
    print M.shape
    print MReduced.shape
    print MReduced[:5, :5]
    if i == 0 or i == 1:
        cl = "rocket"
        vmaxval=2500
    else:
        cl= "mako"
        vmaxval = 5000
    sns.heatmap(MReduced.T, ax=axs[j,i], yticklabels=False, vmin=0, vmax=vmaxval, square=True, cbar = True, cmap=cl)
    i += 1

figname = "Heatmap_plots/Tile-%d-%s.pdf"%(tileno, color)
plt.savefig(figname)
