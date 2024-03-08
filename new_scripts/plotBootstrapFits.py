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

def scatterHue(dff, x, y, key, nbins, ul = None, ll = None, q1=0.01, q2=0.5, q3=0.99):
    df = dff.copy()
    xkey = x
    ykey = y

    # Calculate bins for hue mapping
    quantiles = df[key].quantile([q1, q2, q3])
    if ul is None:
        print "     %scatterHue: detected quantiles: %4.2f, %4.2f, %4.2f"%(quantiles[0.01], quantiles[0.5], quantiles[0.99])

    else:
        quantiles[0.01] = ll
        quantiles[0.99] = ul

    df = df.loc[(df[key]>=quantiles[0.01]) & (df[key]<=quantiles[0.99])]
    x = df[x]
    y = df[y]
    hue = df[key]
    flag = hue.copy()

    huemax = hue.max()
    huemin = hue.min()
    huerange_lb = np.linspace(huemin, huemax, num=nbins+1)[:-1]
    huerange_ub = np.zeros(shape=(nbins,))
    huerange_ub[0:nbins-1] = huerange_lb[1:]
    huerange_ub[-1:] = np.array([huemax])
    np.set_printoptions(precision=2, suppress=True)
    range_label = []
    print "     %scatterHue: using detected ranges for Hue (%4.2f, %4.2f)" % (huemin, huemax)
    print "     %scatterHue: using bins for Hue %d" %(nbins)
    for i in range(nbins):
        range_label.append("%4.2f-%4.2f"%(huerange_lb[i], huerange_ub[i]))
        #print "\t\t%4.2f-%4.2f"%(huerange_lb[i], huerange_ub[i])
        flag.loc[(hue >= huerange_lb[i])  & (hue <= huerange_ub[i])] = "%4.2f-%4.2f"%(huerange_lb[i], huerange_ub[i])

    df['hue'] = flag

    # Compute normalized colorbar for matplotlib:
    norm = plt.Normalize(huemax, huemin)
    sm = plt.cm.ScalarMappable(cmap="icefire", norm=norm)
    sm.set_array([])

    # Plotting using lmplot:
    ax = sns.lmplot(fit_reg=False, data=df, x=xkey, y=ykey, hue='hue', palette='icefire', scatter_kws = {'s':7, 'alpha':0.2}, legend=False)

    return sm

if len(sys.argv) != 3:
    print "Usage: python plotBootstrapFits.py Bootstrap.CPvariant.gz plotdir"
    sys.exit(0)

plotdir = sys.argv[3]
# Read CPvariant file - file of bootstrapped fits, and format it
print 'Reading in data'
df = pd.read_csv(gzip.open(sys.argv[1], "rb"), sep="\t")
df = df.rename(columns = {'Unnamed: 0': 'variant'})
df['Kd'] = np.exp(df['dG']/RT)/conc 
df['logKd'] = np.log10(df['Kd'])
df['Kd_init'] = np.exp(df['dG_init']/RT)/conc 
df['logKd_init'] = np.log10(df['Kd_init'])
df['dG_CI'] = df['dG_ub'] - df['dG_lb']

# Plot initial Kd vs final Kd hued on Fmax and Fmin
print 'plotting initial vs final Kd (Fmax hue)'
x = df.loc[(df['logKd_init'] <= 6) & (df['fmax_init'] <= 6)].copy()
sm = %scatterHue(x, 'logKd_init', 'logKd', 'fmax_init', 50)
plt.xlim([0,6])
plt.ylim([0,6])
plt.xlabel('Initial log10 Kd (nM)')
plt.ylabel('Final log10 Kd (nM)')
cbar = plt.colorbar(sm)
cbar.ax.set_ylabel('Initial Fmax')
plt.savefig('%s/BSF_Kd_final_vs_init_fmaxhue.pdf'%(plotdir))
plt.clf()
del x

print 'plotting initial vs final Kd (Fmin hue)'
x = df.loc[(df['logKd_init'] <= 6) & (df['fmin_init'] <= 0.5)].copy()
sm = %scatterHue(x, 'logKd_init', 'logKd', 'fmin_init', 50)
plt.xlim([0,6])
plt.ylim([0,6])
plt.xlabel('Initial log10 Kd (nM)')
plt.ylabel('Final log10 Kd (nM)')
cbar = plt.colorbar(sm)
cbar.ax.set_ylabel('Initial Fmin')
plt.savefig('%s/BSF_Kd_final_vs_init_fminhue.pdf'%(plotdir))
plt.clf()
del x

# Plot histogram of  log Kds
print 'plotting histogram of final log Kds'
x = df.loc[df['logKd'] <= 8]
sns.distplot(x['logKd'], kde=False)
plt.xlabel('Final log10 Kd (nM)')
plt.ylabel('no. of variants')
plt.xlim([0, 8])
plt.savefig('%s/BSF_Kd_hist.pdf'%(plotdir))
plt.clf()

# Plot histogram of  dG Error (CIs)
print 'plotting histogram of final dG Error estimates (Kd < ~5000)'
x = df.loc[df['logKd'] <= 3.7]
sns.distplot(x['dG_CI'], kde=False)
plt.xlabel('dG Error - 95% CI (nM)')
plt.ylabel('no. of variants')
#plt.xlim([0, 5])
plt.savefig('%s/BSF_dG_Error_hist.pdf'%(plotdir))
plt.clf()

# Plot initial Fmax vs logKd 
print 'plotting scatter plot of Fmax vs Kd'
sns.lmplot(fit_reg=False, data=df, x='logKd', y='fmax', markers='.', scatter_kws = {'color':'blue', 's':7, 'alpha':0.2})
plt.xlabel('Final log10 Kd (nM)')
plt.ylabel('Final Fmax')
plt.xlim([0, 6])
plt.ylim([0, 5])
plt.savefig('%s/BSF_Fmax_vs_Kd.pdf'%(plotdir))
plt.clf()

# Plot initial Fmin vs logKd sns.lmplot(fit_reg=False, data=df, x='logKd', y='fmin', markers='.', scatter_kws = {'color':'blue', 's':7, 'alpha':0.2})
print 'plotting scatter plot of Fmin vs Kd'
plt.xlabel('Final log10 Kd (nM)')
plt.ylabel('Final Fmin')
plt.xlim([0, 6])
plt.ylim([0, 1])
plt.savefig('%s/BSF_Fmin_vs_Kd.pdf'%(plotdir))
plt.clf()

# Plot initial RMSE vs logKd
print 'plotting scatterplot of RMSE vs Kd'
sns.lmplot(fit_reg=False, data=df, x='logKd', y='rmse', markers='.', scatter_kws = {'color':'blue', 's':7, 'alpha':0.2})
plt.xlabel('Final log10 Kd (nM)')
plt.ylabel('Final RMSE')
plt.xlim([0, 6])
plt.ylim([0, 5])
plt.savefig('%s/BSF_rmse_vs_Kd.pdf'%(plotdir))
plt.clf()

# Plot initial R^2 vs logKd
print 'plotting scatter plot of R^2 vs Kd'
sns.lmplot(fit_reg=False, data=df, x='logKd', y='rsq', markers='.', scatter_kws = {'color':'blue', 's':7, 'alpha':0.2})
plt.xlabel('Final log10 Kd (nM)')
plt.ylabel('Final R^2')
plt.xlim([0, 6])
plt.ylim([0, 1])
plt.savefig('%s/BSF_R2_vs_Kd.pdf'%(plotdir))
plt.clf()

