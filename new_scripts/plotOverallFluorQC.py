# Usage: python plotOverallFLuorQC.py [all.annot.green.CPseries] [all.annot.red.CPseries] [conditions.txt] [plotdir]
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import gzip
import sys

if len(sys.argv) != 5:
    print "Usage python plotOverallFLuorQC.py [all.annot.green.CPseries] [all.annot.red.CPseries] [conditions.txt] [plotdir]"
    sys.exit(0)

greenCPseries = sys.argv[1]
redCPseries = sys.argv[2]
conditions = sys.argv[3]
plotdir = sys.argv[4]

fin = open(conditions, 'r')
conds = []
for line in fin:
    cond = line.strip()
    conds.append(cond)

print "Identified Experimental Conditions: "
print conds

print "\nReading green annotated CPseries file: %s" %(greenCPseries)
greenCPseries = pd.read_csv(gzip.open(greenCPseries, "rb"), sep="\t")
print greenCPseries.head(10)
print "\nReading red annotated CPseries file: %s" %(redCPseries)
redCPseries = pd.read_csv(gzip.open(redCPseries, "rb"), sep="\t")
print redCPseries.head(10)

# Axis Labels for future plots:
Xlab = "Exp Conditions"
Ylab = "Cluster Intensity"

########## Plot anyRNA fluorescence trends:
print "\n\n\n*********** Analyzing anyRNA Clusters: "

# GREEN:
print "Plotting Green Fluorescence trends of anyRNA clusters: "
df = greenCPseries.loc[greenCPseries['clusterFlag'] == "anyRNA"]
print "Plotting data from dataframe: "
print df.head(10)
df = df.loc[:, ~df.columns.isin(['clusterID', 'clusterFlag', 'variant_number'])]

# John Edits 6/13/21
# df.columns = conds

# sns.violinplot(x="variable", y="value", data=pd.melt(df), palette="Greens", scale="width")
# plt.xlabel(Xlab)
# plt.ylabel(Ylab)
# plt.savefig('%s/FluorQC_anyRNA_Green_ViolinPlot.pdf'%(plotdir))
# plt.clf()
# sns.boxplot(x="variable", y="value", data=pd.melt(df), palette="Greens", showfliers=False)
# plt.xlabel(Xlab)
# plt.ylabel(Ylab)
# plt.savefig('%s/FluorQC_anyRNA_Green_BoxPlot.pdf'%(plotdir))
# plt.clf()
df.columns = conds

sns.boxplot(x="variable", y="value", data=pd.melt(df), palette="Greens", showfliers=False)
plt.xlabel(Xlab)
plt.xticks(rotation = 60,ha='right')
plt.ylabel(Ylab)
ymin, ymax = plt.gca().get_ylim()
plt.tight_layout()
plt.savefig('%s/FluorQC_anyRNA_Green_BoxPlot.pdf'%(plotdir))
plt.clf()

sns.violinplot(x="variable", y="value", data=pd.melt(df.clip(upper=ymax)), palette="Greens", scale="width")
plt.xlabel(Xlab)
plt.xticks(rotation = 60,ha='right')
plt.ylabel(Ylab)
plt.ylim((ymin,ymax))
plt.tight_layout()
plt.savefig('%s/FluorQC_anyRNA_Green_ViolinPlot.pdf'%(plotdir))
plt.clf()

# End John edits

# RED:
print "Plotting Red Fluorescence trends of anyRNA clusters: "
df = redCPseries.loc[redCPseries['clusterFlag'] == "anyRNA"]
df = df.loc[:, ~df.columns.isin(['clusterID', 'clusterFlag', 'variant_number'])]

# John Edits 6/13/21
# df.columns = conds

# sns.violinplot(x="variable", y="value", data=pd.melt(df), palette="Reds", scale="width")
# plt.xlabel(Xlab)
# plt.ylabel(Ylab)
# plt.savefig('%s/FluorQC_anyRNA_Red_ViolinPlot.pdf'%(plotdir))
# plt.clf()
# sns.boxplot(x="variable", y="value", data=pd.melt(df), palette="Reds", showfliers=False)
# plt.xlabel(Xlab)
# plt.ylabel(Ylab)
# plt.savefig('%s/FluorQC_anyRNA_Red_BoxPlot.pdf'%(plotdir))
# plt.clf()

df.columns = conds

sns.boxplot(x="variable", y="value", data=pd.melt(df), palette="Reds", showfliers=False)
plt.xlabel(Xlab)
plt.xticks(rotation = 60,ha='right')
plt.ylabel(Ylab)
ymin, ymax = plt.gca().get_ylim()
plt.tight_layout()
plt.savefig('%s/FluorQC_anyRNA_Red_BoxPlot.pdf'%(plotdir))
plt.clf()

sns.violinplot(x="variable", y="value", data=pd.melt(df.clip(upper=ymax)), palette="Reds", scale="width")
plt.xlabel(Xlab)
plt.xticks(rotation = 60,ha='right')
plt.ylabel(Ylab)
plt.ylim((ymin,ymax))
plt.tight_layout()
plt.savefig('%s/FluorQC_anyRNA_Red_ViolinPlot.pdf'%(plotdir))
plt.clf()

# End John edits


###########################################
########## Plot FID fluorescence trends:
print "\n\n\n*********** Analyzing FID Clusters: "

# GREEN:
print "Plotting Green Fluorescence trends of FID clusters: "
df = greenCPseries.loc[greenCPseries['clusterFlag'] == "FID"]
df = df.loc[:, ~df.columns.isin(['clusterID', 'clusterFlag', 'variant_number'])]

# John edits 6/13/21
# df.columns = conds
#
# sns.violinplot(x="variable", y="value", data=pd.melt(df), palette="Greens", scale="width")
# plt.xlabel(Xlab)
# plt.ylabel(Ylab)
# plt.savefig('%s/FluorQC_FID_Green_ViolinPlot.pdf'%(plotdir))
# plt.clf()
# sns.boxplot(x="variable", y="value", data=pd.melt(df), palette="Greens", showfliers=False)
# plt.xlabel(Xlab)
# plt.ylabel(Ylab)
# plt.savefig('%s/FluorQC_FID_Green_BoxPlot.pdf'%(plotdir))

df.columns = conds

sns.boxplot(x="variable", y="value", data=pd.melt(df), palette="Greens", showfliers=False)
plt.xlabel(Xlab)
plt.xticks(rotation = 60,ha='right')
plt.ylabel(Ylab)
ymin, ymax = plt.gca().get_ylim()
plt.tight_layout()
plt.savefig('%s/FluorQC_FID_Green_BoxPlot.pdf'%(plotdir))
plt.clf()

sns.violinplot(x="variable", y="value", data=pd.melt(df.clip(upper=ymax)), palette="Greens", scale="width")
plt.xlabel(Xlab)
plt.xticks(rotation = 60,ha='right')
plt.ylabel(Ylab)
plt.ylim((ymin,ymax))
plt.tight_layout()
plt.savefig('%s/FluorQC_FID_Green_ViolinPlot.pdf'%(plotdir))
plt.clf()

# End John edits


# RED:
print "Plotting Red Fluorescence trends of FID clusters: "
df = redCPseries.loc[redCPseries['clusterFlag'] == "FID"]
df = df.loc[:, ~df.columns.isin(['clusterID', 'clusterFlag', 'variant_number'])]

# John edits 6/13/21
# df.columns = conds
#
# sns.violinplot(x="variable", y="value", data=pd.melt(df), palette="Reds", scale="width")
# plt.xlabel(Xlab)
# plt.ylabel(Ylab)
# plt.savefig('%s/FluorQC_FID_Red_ViolinPlot.pdf'%(plotdir))
# plt.clf()
# sns.boxplot(x="variable", y="value", data=pd.melt(df), palette="Reds", showfliers=False)
# plt.xlabel(Xlab)
# plt.ylabel(Ylab)
# plt.savefig('%s/FluorQC_FID_Red_BoxPlot.pdf'%(plotdir))

df.columns = conds

sns.boxplot(x="variable", y="value", data=pd.melt(df), palette="Reds", showfliers=False)
plt.xlabel(Xlab)
plt.xticks(rotation = 60,ha='right')
plt.ylabel(Ylab)
ymin, ymax = plt.gca().get_ylim()
plt.tight_layout()
plt.savefig('%s/FluorQC_FID_Red_BoxPlot.pdf'%(plotdir))
plt.clf()

sns.violinplot(x="variable", y="value", data=pd.melt(df.clip(upper=ymax)), palette="Reds", scale="width")
plt.xlabel(Xlab)
plt.xticks(rotation = 60,ha='right')
plt.ylabel(Ylab)
plt.ylim((ymin,ymax))
plt.tight_layout()
plt.savefig('%s/FluorQC_FID_Red_ViolinPlot.pdf'%(plotdir))
plt.clf()

# End John edits

sys.exit(0)

###########################################
########## Plot Background fluorescence trends:
print "\n\n\n*********** Analyzing Background Clusters: "

# GREEN:
print "Plotting Green Fluorescence trends of BG clusters: "
df = greenCPseries.loc[(greenCPseries['clusterFlag'] != "FID") & (greenCPseries['clusterFlag'] != "anyRNA")]
df = df.loc[:, ~df.columns.isin(['clusterID', 'clusterFlag', 'variant_number'])]
df.columns = conds

sns.violinplot(x="variable", y="value", data=pd.melt(df), palette="Greens", scale="width")
plt.xlabel(Xlab)
plt.ylabel(Ylab)
plt.savefig('%s/FluorQC_BG_Green_ViolinPlot.pdf'%(plotdir))
plt.clf()
sns.boxplot(x="variable", y="value", data=pd.melt(df), palette="Greens", showfliers=False)
plt.xlabel(Xlab)
plt.ylabel(Ylab)
plt.savefig('%s/FluorQC_BG_Green_BoxPlot.pdf'%(plotdir))

# RED:
print "Plotting Red Fluorescence trends of FID clusters: "
df = redCPseries.loc[(redCPseries['clusterFlag'] != "FID") & (redCPseries['clusterFlag'] != "anyRNA")]
df = df.loc[:, ~df.columns.isin(['clusterID', 'clusterFlag', 'variant_number'])]
df.columns = conds

sns.violinplot(x="variable", y="value", data=pd.melt(df), palette="Reds", scale="width")
plt.xlabel(Xlab)
plt.ylabel(Ylab)
plt.savefig('%s/FluorQC_BG_Red_ViolinPlot.pdf'%(plotdir))
plt.clf()
sns.boxplot(x="variable", y="value", data=pd.melt(df), palette="Reds", showfliers=False)
plt.xlabel(Xlab)
plt.ylabel(Ylab)
plt.savefig('%s/FluorQC_BG_Red_BoxPlot.pdf'%(plotdir))

sys.exit(0)
