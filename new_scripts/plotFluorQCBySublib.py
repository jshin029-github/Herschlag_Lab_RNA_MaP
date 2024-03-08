# Usage: python plotOverallFLuorQC.py [all.annot.green.CPseries] [all.annot.red.CPseries] [conditions.txt] [libchar] [plotdir]
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import gzip
import sys

if len(sys.argv) != 6:
    print "Usage python plotFluorQCBySublib.py [all.annot.green.CPseries] [all.annot.red.CPseries] [conditions.txt] [libchar] [plotdir]"
    sys.exit(0)

greenCPseries = sys.argv[1]
redCPseries = sys.argv[2]
conditions = sys.argv[3]
libchar = sys.argv[4]
plotdir = sys.argv[5]

fin = open(conditions, 'r')
conds = []
for line in fin:
    cond = line.strip()
    conds.append(cond)

print "Identified Experimental Conditions: "
print conds

lc = pd.read_csv(libchar)
sublibs = lc['body_plan'].unique()
sublib_vars = {}
for sublib in sublibs:
    lcsublib = lc.loc[lc['body_plan'] == sublib]
    var = np.array(lcsublib['unique_variant'])
    sublib_vars[sublib] = var

print "\nIdentified Sublibraries (Body Plans): "
print sublibs

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
dff = greenCPseries.loc[greenCPseries['clusterFlag'] == "anyRNA"]
#dff = greenCPseries.merge(cpannot, on="clusterID", how="left")
print "\nPlotting data from dataframe: \n\n"
print dff.head(10)


# Iterate over sublibs:
for sublib in sublibs:
    print "\n\n     Plotting sublib : %s"%(sublib)
    variant_ids = sublib_vars[sublib]
    all_ids = np.array(dff['variant_number'])
    mask = np.in1d(all_ids, variant_ids)
    df = dff.loc[mask]
    print "         No of Library Variants : %d"%(variant_ids.shape[0])
    print "         No of Detected Clusters : %d"%(df.shape[0])

    if df.shape[0] != 0:
        df = df.loc[:, ~df.columns.isin(['clusterID', 'clusterFlag', 'variant_number'])]

        # John edits 6/13/21
        # df.columns = conds
        #
        # sns.violinplot(x="variable", y="value", data=pd.melt(df), palette="Greens", scale="width")
        # plt.xlabel(Xlab)
        # plt.ylabel(Ylab)
        # print "         Saving plot: %s/SublibFluorQC_%s_Green_ViolinPlot.pdf"%(plotdir,sublib)
        # plt.savefig('%s/SublibFluorQC_%s_Green_ViolinPlot.pdf'%(plotdir,sublib))
        # plt.clf()
        # sns.boxplot(x="variable", y="value", data=pd.melt(df), palette="Greens", showfliers=False)
        # plt.xlabel(Xlab)
        # plt.ylabel(Ylab)
        # print "         Saving plot: %s/SublibFluorQC_%s_Green_BoxPlot.pdf"%(plotdir,sublib)
        # plt.savefig('%s/SublibFluorQC_%s_Green_BoxPlot.pdf'%(plotdir,sublib))
        # plt.clf()

        try:
            df.columns = conds
        except ValueError:
            break

        sns.boxplot(x="variable", y="value", data=pd.melt(df), palette="Greens", showfliers=False)
        plt.xlabel(Xlab)
        plt.xticks(rotation = 60,ha='right')
        plt.ylabel(Ylab)
        print "         Saving plot: %s/SublibFluorQC_%s_Green_BoxPlot.pdf"%(plotdir,sublib)
        ymin, ymax = plt.gca().get_ylim()
        plt.tight_layout()
        plt.savefig('%s/SublibFluorQC_%s_Green_BoxPlot.pdf'%(plotdir,sublib))
        plt.clf()

        sns.violinplot(x="variable", y="value", data=pd.melt(df.clip(upper=ymax)), palette="Greens", scale="width")
        plt.xlabel(Xlab)
        plt.xticks(rotation = 60,ha='right')
        plt.ylabel(Ylab)
        plt.ylim((ymin,ymax))
        plt.tight_layout()
        print "         Saving plot: %s/SublibFluorQC_%s_Green_ViolinPlot.pdf"%(plotdir,sublib)
        plt.savefig('%s/SublibFluorQC_%s_Green_ViolinPlot.pdf'%(plotdir,sublib))
        plt.clf()

        # End John edits


# RED:
print "\n\nPlotting Red Fluorescence trends of anyRNA clusters: \n"
dff = redCPseries.loc[redCPseries['clusterFlag'] == "anyRNA"]

# Iterate over sublibs:
for sublib in sublibs:
    print "\n\n     Plotting sublib : %s"%(sublib)
    variant_ids = sublib_vars[sublib]
    all_ids = np.array(dff['variant_number'])
    mask = np.in1d(all_ids, variant_ids)
    df = dff.loc[mask]
    print "         No of Library Variants : %d"%(variant_ids.shape[0])
    print "         No of Detected Clusters : %d"%(df.shape[0])

    if df.shape[0] != 0:
        df = df.loc[:, ~df.columns.isin(['clusterID', 'clusterFlag', 'variant_number'])]

        # John edits 6/13/21

        # df.columns = conds
        #
        # sns.violinplot(x="variable", y="value", data=pd.melt(df), palette="Reds", scale="width")
        # plt.xlabel(Xlab)
        # plt.ylabel(Ylab)
        # print "         Saving plot: %s/SublibFluorQC_%s_Red_ViolinPlot.pdf"%(plotdir,sublib)
        # plt.savefig('%s/SublibFluorQC_%s_Red_ViolinPlot.pdf'%(plotdir,sublib))
        # plt.clf()
        # sns.boxplot(x="variable", y="value", data=pd.melt(df), palette="Reds", showfliers=False)
        # plt.xlabel(Xlab)
        # plt.ylabel(Ylab)
        # print "         Saving plot: %s/SublibFluorQC_%s_Red_BoxPlot.pdf"%(plotdir,sublib)
        # plt.savefig('%s/SublibFluorQC_%s_Red_BoxPlot.pdf'%(plotdir,sublib))
        # plt.clf()

        try:
            df.columns = conds
        except ValueError:
            break

        sns.boxplot(x="variable", y="value", data=pd.melt(df), palette="Reds", showfliers=False)
        plt.xlabel(Xlab)
        plt.xticks(rotation = 60,ha='right')
        plt.ylabel(Ylab)
        print "         Saving plot: %s/SublibFluorQC_%s_Red_BoxPlot.pdf"%(plotdir,sublib)
        ymin, ymax = plt.gca().get_ylim()
        plt.tight_layout()
        plt.savefig('%s/SublibFluorQC_%s_Red_BoxPlot.pdf'%(plotdir,sublib))
        plt.clf()
        sns.violinplot(x="variable", y="value", data=pd.melt(df.clip(upper=ymax)), palette="Reds", scale="width")
        plt.xlabel(Xlab)
        plt.xticks(rotation = 60,ha='right')
        plt.ylabel(Ylab)
        plt.ylim((ymin,ymax))
        plt.tight_layout()
        print "         Saving plot: %s/SublibFluorQC_%s_Red_ViolinPlot.pdf"%(plotdir,sublib)
        plt.savefig('%s/SublibFluorQC_%s_Red_ViolinPlot.pdf'%(plotdir,sublib))
        plt.clf()
