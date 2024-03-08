import pandas as pd
import numpy as np
import os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import skimage

# Usage: python makeClusterHeatmap.py [CPseries Template] [Tile No] [color]

#### Change the following variables before running: ####
# Set folder names of green images
green_folders = [
"13_equilibrium0p91nM_green",
"14_equilibrium2p7nM_green",
"15_equilibrium8p2nM_green", 
"16_equilibrium25nM_green",
"17_equilibrium74nM_green",
"18_equilibrium222nM_green",
"19_equilibrium666nM_green",
"20_equilibrium2000nM_green"
        ]
# Set folder names of red images
red_folders = [
"13_equilibrium0p91nM_red",
"14_equilibrium2p7nM_red" ,
"15_equilibrium8p2nM_red", 
"16_equilibrium25nM_red",
"17_equilibrium74nM_red",
"18_equilibrium222nM_red",
"19_equilibrium666nM_red",
"20_equilibrium2000nM_red"
        ]
# Set flow piece concentrations corresponding to images
fpconcs = [ 2.7, 8.2, 8.2, 25, 74, 222, 666, 2000]
# Set Resolution within which to average (e.g. setting to 400 => avg within 400 X 400 blocks)
Resolution = 400
##################################################


color = sys.argv[-1]
if color == "green":
    folders = green_folders
else:
    folders = red_folders

Template = sys.argv[1]
NTile = int(sys.argv[2])

def get_mask(df, colname, clustype):
    if clustype == "FID" or clustype == "anyRNA":
        mask = (df[colname] == clustype)
    else:
        mask = (df[colname] != "FID") & (df[colname] != "anyRNA")
    return mask

for tileno in [NTile]:
    j = 0
    fig,axs = plt.subplots(nrows=len(folders), ncols=3, sharex =True, sharey=True, figsize=(20, 35))

    for folder, fpconc in zip(folders, fpconcs):
        # Read in current tile data:
        fn = Template%(tileno)
        fullfn = "%s/%s"%(folder, fn)
        if not os.path.isfile(fullfn):
            j += 1
            continue
        print "Reading in file: %s" % (fn)
        df = pd.read_csv(fullfn, sep="\t")
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
            #axs[i].set_title(clustype)
            #plt.savefig("test.png")
            i += 1
        #plt.suptitle("Tile = %d, Flow piece = %3.2f nM"%(tileno, FPconc))
        #figname = "Heatmap_plots/Colorbars.pdf"
        #plt.savefig(figname)
        j += 1
        #Fpconcname = str(fpconc).replace('.', 'p')
        #figname = "Heatmap_plots/Tile-%d-FP-%s-%s"%(tileno, Fpconcname, color)
    figname = "Heatmap_plots/Tile-%d-%s.pdf"%(tileno, color)
    plt.savefig(figname)
