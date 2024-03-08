import pandas as pd, numpy as np, gzip, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import skimage


##### Usage: python fixBadClusters.py [CPseries file] [Background Cutoff] [FID Cutoff] [output csv file] [output plot]
FidCutoff = float(sys.argv[2])
BGCutoff = float(sys.argv[3])
csvout = sys.argv[4]
plotout = sys.argv[5]

######### READ IN DATA ###################################

tile1 = pd.read_csv(sys.argv[1], sep="\t", header=None)
tile1 = tile1.iloc[:, [0, 1, -1]]
tile1.columns = ["clusterID", "clusterType", "Fluorescence"]
print "Tile 1: total number of clusters = ", tile1.shape[0]
print "Tile 1: number of clusters quantified = ", tile1.dropna(subset=["Fluorescence"]).shape[0]
df1 = tile1.loc[tile1["clusterType"] == "anyRNA"]
df2 = tile1.loc[tile1["clusterType"] == "FID"]
df3 = tile1.loc[(tile1["clusterType"] != "anyRNA") & (tile1["clusterType"] != "FID")]

##########################################################


######## DETECT OUTLIERS USING BLOCK REDUCE  #############

dfbg = df3.loc[:, ['clusterID', 'Fluorescence']]
dfbg = dfbg.dropna()
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

dffid = df2.loc[:, ['clusterID', 'Fluorescence']]
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

dfrna = df1.loc[:, ['clusterID', 'Fluorescence']]
dfrna = dfrna.dropna()
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

Resolution = 400
Xmax = max(dfbg['X'].max(), dfrna['X'].max(), dffid['X'].max())
Ymax = max(dfbg['Y'].max(), dfrna['Y'].max(), dffid['Y'].max())
N = max(Xmax, Ymax) + Resolution - max(Xmax, Ymax)%Resolution
print "Using matrix size: %d"%(N)

dfs = [dfbg, dffid]
clusterTypes = ["background", "Fiducial"]


for df, clustype in zip(dfs, clusterTypes):
    M = np.empty((N,N,),)
    M[:] = np.nan
    for index, row in df.iterrows():
        M[row['X'], row['Y']] = row['Fluorescence']
    print "M Shape : ", M.shape

    print np.nanmean(M)
    MReduced = skimage.measure.block_reduce(M, (Resolution, Resolution), np.nanmedian)
    p = MReduced.shape[0]
    print "MReduced Shape : ", MReduced.shape
    MReduced[np.isnan(MReduced)] = -50
    if clustype == 'Fiducial':
        FidOutliers = np.where((MReduced < FidCutoff) & (MReduced != -50))
        print "Fiducial Cluster Outliers: ", len(FidOutliers)
    elif clustype == 'background':
        BGOutliers = np.where((MReduced > BGCutoff) & (MReduced != -50))
        print "Background Cluster Outliers: ", len(BGOutliers)


#####################################################

######### PLOT DETECTED OUTLIERS ####################

M = np.empty((p,p,),)
M[:] = np.nan
M[FidOutliers] = 1000
M[BGOutliers] = -1000
sns.heatmap(M.T,  vmin=-1000, vmax=1000, square=True, cbar = True, cmap="rocket")
plt.savefig(plotout)


#####################################################

######### EXTRACT INDICES OF OUTLIERS ###############


BGOutliersX = BGOutliers[0]
BGOutliersY = BGOutliers[1]
FidOutliersX = FidOutliers[0]
FidOutliersY = FidOutliers[1]

BGOutliersXLowerBound = BGOutliersX * Resolution
BGOutliersXUpperBound = BGOutliersX * Resolution + Resolution - 1
BGOutliersYLowerBound = BGOutliersY * Resolution
BGOutliersYUpperBound = BGOutliersY * Resolution + Resolution - 1
FidOutliersXLowerBound = FidOutliersX * Resolution
FidOutliersXUpperBound = FidOutliersX * Resolution + Resolution - 1
FidOutliersYLowerBound = FidOutliersY * Resolution
FidOutliersYUpperBound = FidOutliersY * Resolution + Resolution - 1

filtered_dfrna = pd.DataFrame()
filtered_dfrna["clusterID"] = []
filtered_dfrna["Fluorescence"] = []
filtered_dfrna["X"] = []
filtered_dfrna["Y"] = []

BGbadClusterList = []
FidbadClusterList = []

for i in range(BGOutliersXLowerBound.shape[0]):
    #print i
    Xlb = BGOutliersXLowerBound[i]
    Xub = BGOutliersXUpperBound[i]
    #print " x: %d - %d"%(Xlb, Xub)

    Ylb = BGOutliersYLowerBound[i]
    Yub = BGOutliersYUpperBound[i]
    #print " Y: %d - %d"%(Ylb, Yub)

    maskX = (dfrna['X'] >= Xlb) & (dfrna['X'] <= Xub)
    maskY = (dfrna['Y'] >= Ylb) & (dfrna['Y'] <= Yub)
    bclist = dfrna.loc[maskX & maskY, 'clusterID']
    BGbadClusterList = BGbadClusterList + bclist.tolist()

print "Bad Clusters Detected using Background Cutoffs = %d"%(len(BGbadClusterList))

for i in range(FidOutliersXLowerBound.shape[0]):
    #print i
    Xlb = FidOutliersXLowerBound[i]
    Xub = FidOutliersXUpperBound[i]
    #print " x: %d - %d"%(Xlb, Xub)

    Ylb = FidOutliersYLowerBound[i]
    Yub = FidOutliersYUpperBound[i]
    #print " Y: %d - %d"%(Ylb, Yub)

    maskX = (dfrna['X'] >= Xlb) & (dfrna['X'] <= Xub)
    maskY = (dfrna['Y'] >= Ylb) & (dfrna['Y'] <= Yub)
    bclist = dfrna.loc[maskX & maskY, 'clusterID']
    FidbadClusterList = FidbadClusterList + bclist.tolist()

print "Bad Clusters Detected using Fiducial Cutoffs = %d"%(len(FidbadClusterList))

BadClusters = pd.DataFrame()
BadClusters["clusterID"] = (BGbadClusterList+FidbadClusterList)
BadClusters.to_csv(csvout, sep="\t", index=False)


##################################################
