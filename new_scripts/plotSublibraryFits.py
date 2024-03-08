# Python script to plot (sns) relevant initial %s to analyze
# the results of Kd bootstrap refinement.
# %s generated in ./%s and named: %s/BSF_*.pdf
# Usage: python plotBootstrapFits.py Bootstrap.CPvariant.gz"

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
    print "Usage python plotLibraryFits.py [cpvar] [libchar] [plotdir]"
    sys.exit()

libchar = sys.argv[2]
plotdir = sys.argv[3]
# Read CPvariant file - file of bootstrapped fits, and format it
print 'Reading in data'
df = pd.read_csv(gzip.open(sys.argv[1], "rb"), sep="\t")
df = df.rename(columns = {'Unnamed: 0': 'variant'})
df['Kd'] = np.exp(df['dG']/RT)/conc 
df['logKd'] = np.log10(df['Kd'])
df['Kd_init'] = np.exp(df['dG_init']/RT)/conc 
df['logKd_init'] = np.log10(df['Kd_init'])
df = df.loc[:, ['variant', 'dG', 'logKd']]

lc = pd.read_csv(libchar)
lc = lc.loc[:, ['variant', 'sublibrary']]

# sublibraries:
sublibs = list(lc['sublibrary'].drop_duplicates())
for sublib in sublibs:
    print "Plotting sublibrary: %s"%(sublib)
    x = lc.loc[lc['sublibrary'] == sublib]
    xx = x.merge(df, on = 'variant', how='inner')
    xx = xx.loc[xx['logKd'] <= 8]
    if xx.shape[0] == 0:
        continue
    kds = xx['logKd']
    dGs = xx['dG']
    filename = "SublibPlots_%s_Kd.pdf"%(sublib)
    sns.distplot(kds, kde=False, bins=50)
    plt.xlabel('Final log10 Kd (nM)')
    plt.ylabel('no. of variants')
    plt.xlim([0, 8])
    plt.savefig('%s/%s'%(plotdir, filename))
    plt.clf()

sys.exit()




# Plot histogram of  log Kds
x = df.loc[df['logKd'] <= 8]
sns.distplot(x['logKd'], kde=False)
plt.xlabel('Final log10 Kd (nM)')
plt.ylabel('no. of variants')
plt.xlim([0, 8])
plt.savefig('%s/BSF_Kd_hist.pdf'%(plotdir))
plt.clf()

# Plot initial Fmax vs logKd 
sns.lmplot(fit_reg=False, data=df, x='logKd', y='fmax', markers='.', scatter_kws = {'color':'blue', 's':7, 'alpha':0.2})
plt.xlabel('Final log10 Kd (nM)')
plt.ylabel('Final Fmax')
plt.xlim([0, 6])
plt.ylim([0, 5])
plt.savefig('%s/BSF_Fmax_vs_Kd.pdf'%(plotdir))
plt.clf()


# Plot initial Fmin vs logKd 
sns.lmplot(fit_reg=False, data=df, x='logKd', y='fmin', markers='.', scatter_kws = {'color':'blue', 's':7, 'alpha':0.2})
plt.xlabel('Final log10 Kd (nM)')
plt.ylabel('Final Fmin')
plt.xlim([0, 6])
plt.ylim([0, 1])
plt.savefig('%s/BSF_Fmin_vs_Kd.pdf'%(plotdir))
plt.clf()

# Plot initial RMSE vs logKd
sns.lmplot(fit_reg=False, data=df, x='logKd', y='rmse', markers='.', scatter_kws = {'color':'blue', 's':7, 'alpha':0.2})
plt.xlabel('Final log10 Kd (nM)')
plt.ylabel('Final RMSE')
plt.xlim([0, 6])
plt.ylim([0, 5])
plt.savefig('%s/BSF_rmse_vs_Kd.pdf'%(plotdir))
plt.clf()

# Plot initial R^2 vs logKd
sns.lmplot(fit_reg=False, data=df, x='logKd', y='rsq', markers='.', scatter_kws = {'color':'blue', 's':7, 'alpha':0.2})
plt.xlabel('Final log10 Kd (nM)')
plt.ylabel('Final R^2')
plt.xlim([0, 6])
plt.ylim([0, 1])
plt.savefig('%s/BSF_R2_vs_Kd.pdf'%(plotdir))
plt.clf()

