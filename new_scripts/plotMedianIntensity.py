import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import gzip
import sys

greenCPseries = sys.argv[1]
redCPseries = sys.argv[2]
normCPseries = sys.argv[3]
plotdir = sys.argv[4]


Xlab = "Flow Piece Conc [nM]"
Ylab = "Median Cluster Intensity"

# Plot Green median intensities:
greenCPseries = pd.read_csv(gzip.open(greenCPseries, "rb"), sep="\t")
greenCPseries = greenCPseries.iloc[:, 1:]
print greenCPseries.head()
greenCPseries.columns = [0.91, 2.7, 8.2,  74, 222, 667, 2000]
sns.boxplot(x="variable", y="value", data=pd.melt(greenCPseries), showfliers=False,palette="Greens" )
plt.xlabel('Flow piece concentration (nM)')
plt.ylabel('Fluorescence Intensity')
plt.savefig('%s/green_median_intensity.pdf'%(plotdir))
plt.clf()

redCPseries = pd.read_csv(gzip.open(redCPseries, "rb"), sep="\t")
redCPseries = redCPseries.iloc[:, 1:]
print redCPseries.head()
redCPseries.columns = ['0.91', '2.7', '8.2', '74', '222', '667', '2000']
redCPseries.columns = [0.91, 2.7, 8.2,  74, 222, 667, 2000]
sns.boxplot(x="variable", y="value", data=pd.melt(redCPseries), showfliers=False, palette="Reds")
plt.xlabel('Flow piece concentration (nM)')
plt.ylabel('Fluorescence Intensity')
plt.savefig('%s/red_median_intensity.pdf'%(plotdir))
plt.clf()

normCPseries = pd.read_csv(gzip.open(normCPseries, "rb"), sep="\t")
normCPseries = normCPseries.iloc[:, 1:]
print normCPseries.head()
normCPseries.columns = [0.91, 2.7, 8.2,  74, 222, 667, 2000]
sns.boxplot(x="variable", y="value", data=pd.melt(normCPseries), showfliers=False, palette="Blues")
plt.xlabel('Flow piece concentration (nM)')
plt.ylabel('Fluorescence Intensity')
plt.savefig('%s/norm_median_intensity.pdf'%(plotdir))
plt.clf()

