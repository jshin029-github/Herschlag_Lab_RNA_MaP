# Usage: python plotSingleClusterFits.py [CPfitted] [CPannot] [plotdir]
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import gzip
import sys
import numpy as np

RT = 0.582
conc = 1E-9

if len(sys.argv) != 4:
    print "Usage: python plotSingleClusterFits.py [CPfitted] [CPannot] [plotdir]"
    sys.exit(0)

# Read CPfitted file - file of single cluster fits - format into median vals:
df = pd.read_csv(gzip.open(sys.argv[1], "rb"), sep="\t")
cpa = pd.read_csv(gzip.open(sys.argv[2], "rb"), sep="\t")
plotdir = sys.argv[3]

print "Read in fitted results :"
print df.head()

print "Read in CPannot file :"
print cpa.head()

#df.columns = ['clusterID', 'dG', 'fmax', 'fmin', 'dG_stde', 'fmax_stde', 'fmin_stde', 'rsq', 'exit_flag', 'rmse']
df = df.rename(columns={'Unnamed: 0':'clusterID'})
print "\n\nSCF before merging with renamed columns:"
print df.head()
dff = df.merge(cpa, on='clusterID', how="left")
#dff = dff.sort_values(by=['variant_number', 'clusterID'])
print "\n\nSCF merged with CPannot:"
print dff.head()
var = dff.groupby(by=['variant_number']).median()
var['cluster_count'] = dff.groupby(by=['variant_number']).size()
var['Kd'] = np.exp(var['dG']/RT)/conc
var['logKd'] = np.log10(var['Kd'])
print "\n\nSCF with Kds:"
print var.head()
#var.to_csv(sys.argv[1][:-3]+'.median.gz', sep='\t', header=True, index=False, compression="gzip")
df = var

# Plot histogram of cluster counts (only plot cluster counts upto 500):
x = df.loc[df['cluster_count'] <= 300]
sns.distplot(x['cluster_count'], kde=False)
#sns.distplot(np.array(x['cluster_count'], dtype=int), kde=False, hist=False, bins=100, color='blue')
#sns.histplot(data=x, x='cluster_count', stat='count', binwidth=5, binrange=[0, 300], color="blue", kde=False)
plt.xlabel('no. clusters per variant')
plt.ylabel('no. of variants')
plt.savefig('%s/SCF_ClusterCounts_hist.png'%(plotdir))
plt.xlim([0,300])
plt.clf()

# Plot histogram of median log Kds
x = df.loc[df['logKd'] <= 8]
sns.distplot(x['logKd'], kde=False)
plt.xlabel('log10 Kd (nM)')
plt.ylabel('no. of variants')
plt.savefig('%s/SCF_Kd_hist.png'%(plotdir))
plt.clf()

# Plot initial Fmax vs logKd 
sns.lmplot(fit_reg=False, data=x, x='logKd', y='fmax', markers='.', scatter_kws = {'color':'blue', 's':7, 'alpha':0.2})
plt.xlabel('initial log10 Kd (nM)')
plt.ylabel('initial Fmax')
plt.xlim([0, 6])
plt.ylim([0, 5])
plt.axvline(x=2.398, linestyle="--", color="k", linewidth=1)
plt.savefig('%s/SCF_Fmax_vs_Kd.png'%(plotdir))
plt.clf()


# Plot initial Fmin vs logKd 
sns.lmplot(fit_reg=False, data=x, x='logKd', y='fmin', markers='.', scatter_kws = {'color':'blue', 's':7, 'alpha':0.2})
plt.xlabel('initial log10 Kd (nM)')
plt.ylabel('initial Fmin')
plt.xlim([0, 6])
plt.ylim([0, 1])
plt.axvline(x=2.398, linestyle="--", color="k", linewidth=1)
plt.savefig('%s/SCF_Fmin_vs_Kd.png'%(plotdir))
plt.clf()

# Plot initial RMSE vs logKd
sns.lmplot(fit_reg=False, data=x, x='logKd', y='rmse', markers='.', scatter_kws = {'color':'blue', 's':7, 'alpha':0.2})
plt.xlabel('initial log10 Kd (nM)')
plt.ylabel('initial RMSE')
plt.xlim([0, 6])
plt.ylim([0, 10])
plt.savefig('%s/SCF_rmse_vs_Kd.png'%(plotdir))
plt.clf()

# Plot initial R^2 vs logKd
sns.lmplot(fit_reg=False, data=x, x='logKd', y='rsq', markers='.', scatter_kws = {'color':'blue', 's':7, 'alpha':0.2})
plt.xlabel('initial log10 Kd (nM)')
plt.ylabel('initial R^2')
plt.xlim([0, 6])
plt.ylim([0, 1])
plt.savefig('%s/SCF_R2_vs_Kd.png'%(plotdir))
plt.clf()

