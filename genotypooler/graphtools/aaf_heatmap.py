import os, sys
import numpy as np
import pandas as pd
import json
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import Divider, Size
import argparse

import warnings
warnings.filterwarnings("ignore")

rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from genotypooler.poolSNPs import dataframe as vcfdf
from genotypooler.poolSNPs.metrics import quality as qual
from genotypooler.persotools.files import *


pop = 'impSTU'
outdir = f'/home/camille/MagicWheat/results/chrom1-uppmax/{pop}'
xdata = '/home/camille/MagicWheat/runs/poolimputeSNPs/results/data/1/STU.Chr1.SNPs.pruned.sorted.vcf.gz'
ydata = '/home/camille/MagicWheat/data/chrom1-uppmax/prophaser/STU.Chr1.SNPs.pruned.sorted.pooled.imputed.vcf.gz'
ydata2 = '/home/camille/MagicWheat/runs/poolimputeSNPs/results/data/1/STU.Chr1.SNPs.pruned.sorted.pooled.vcf.gz'
# '/home/camille/MagicWheat/runs/poolimputeSNPs/results/data/1/PNL.Chr1.SNPs.pruned.sorted.vcf.gz'
# '/home/camille/MagicWheat/data/chrom1-uppmax/prophaser/STU.Chr1.SNPs.pruned.sorted.pooled.imputed.vcf.gz'
# '/home/camille/MagicWheat/runs/poolimputeSNPs/results/data/1/STU.Chr1.SNPs.pruned.sorted.pooled.vcf.gz'
# None  # ydata2

bS = 0.01
rQ = 5

if not os.path.exists(outdir):
    os.mkdir(outdir)

dfX = vcfdf.PandasMixedVCF(xdata, format='GT', indextype='chrom:pos')
try:
    dfY = vcfdf.PandasMixedVCF(ydata, format='GL', indextype='chrom:pos')
    gY = dfY.trinary_encoding()
except KeyError:
    dfY = vcfdf.PandasMixedVCF(ydata, format='GT', indextype='chrom:pos')
    gY = dfY.trinary_encoding()

aafX = dfX.aaf.rename(columns={'aaf': 'aafX'})
aafY = dfY.aaf.rename(columns={'aaf': 'aafY'})

aafXY = aafX.join(aafY)
aafXY = aafXY.sort_values(by=['aafX', 'aafY'], ascending=True).rolling(window=rQ).mean()
print(aafXY)

dfroll = qual.QuantilesDataFrame(aafXY.aafX.to_frame(),
                                 aafXY.aafY,
                                 bins_step=bS)
rollXY = dfroll.binnedX.join(dfroll.binnedY).round(decimals=3)
rollXY.reset_index(drop=True, inplace=True)
print(rollXY.head(15))

counts_per_Xbin = rollXY.groupby(by='binned_aafX', as_index=False, dropna=False).count()\
    .rename(columns={'binned_aafY': 'counts_per_Xbin'})

counts_per_XYbin = rollXY
counts_per_XYbin['counts_per_binXY'] = 1
counts_per_XYbin = counts_per_XYbin.groupby(by=['binned_aafX', 'binned_aafY'], as_index=False, dropna=False).count()
counts_per_XYbin = counts_per_XYbin.merge(counts_per_Xbin, on='binned_aafX')
counts_per_XYbin['prop_per_binXY'] = counts_per_XYbin['counts_per_binXY'] / counts_per_XYbin['counts_per_Xbin']
counts_per_XYbin = counts_per_XYbin.dropna()\
    .reset_index(drop=True)\
    .drop(['counts_per_binXY', 'counts_per_Xbin'], axis=1)
counts_per_XYbin['binned_aafX'] = [float(b) for b in counts_per_XYbin['binned_aafX']]
counts_per_XYbin['binned_aafY'] = [float(b) for b in counts_per_XYbin['binned_aafY']]
print(counts_per_XYbin.head(15))
print(counts_per_XYbin.tail(15))

if ydata2 is not None:
    try:
        dfY2 = vcfdf.PandasMixedVCF(ydata2, format='GL', indextype='chrom:pos')
        gY2 = dfY2.trinary_encoding()
    except KeyError:
        dfY2 = vcfdf.PandasMixedVCF(ydata2, format='GT', indextype='chrom:pos')
        gY2 = dfY2.trinary_encoding()

    aafY2 = dfY2.aaf.rename(columns={'aaf': 'aafY'})
    aafXY2 = aafX.join(aafY2)
    aafXY2 = aafXY2.sort_values(by=['aafX', 'aafY'], ascending=True).rolling(window=rQ).mean()

    dfroll2 = qual.QuantilesDataFrame(aafXY2.aafX.to_frame(),
                                      aafXY2.aafY,
                                      bins_step=bS)
    rollXY2 = dfroll2.binnedX.join(dfroll2.binnedY).round(decimals=3)
    rollXY2.reset_index(drop=True, inplace=True)

    counts_per_Xbin2 = rollXY2.groupby(by='binned_aafX', as_index=False, dropna=False).count()\
        .rename(columns={'binned_aafY': 'counts_per_Xbin'})

    counts_per_XYbin2 = rollXY2
    counts_per_XYbin2['counts_per_binXY'] = 1
    counts_per_XYbin2 = counts_per_XYbin2.groupby(by=['binned_aafX', 'binned_aafY'], as_index=False, dropna=False).count()
    counts_per_XYbin2 = counts_per_XYbin2.merge(counts_per_Xbin2, on='binned_aafX')
    counts_per_XYbin2['prop_per_binXY'] = counts_per_XYbin2['counts_per_binXY'] / counts_per_XYbin2['counts_per_Xbin']
    counts_per_XYbin2 = counts_per_XYbin2.dropna()\
        .reset_index(drop=True)\
        .drop(['counts_per_binXY', 'counts_per_Xbin'], axis=1)
    counts_per_XYbin2['binned_aafX'] = [float(b) for b in counts_per_XYbin2['binned_aafX']]
    counts_per_XYbin2['binned_aafY'] = [float(b) for b in counts_per_XYbin2['binned_aafY']]
    print(counts_per_XYbin2.head(15))
    print(counts_per_XYbin2.tail(15))
    counts_per_XYbin_merged = counts_per_XYbin.merge(counts_per_XYbin2, on=['binned_aafX', 'binned_aafY'], how='outer', suffixes=(None, '2'))

try:
    heatdata = counts_per_XYbin_merged.pivot(index='binned_aafY', columns='binned_aafX', values='prop_per_binXY')
except:
    heatdata = counts_per_XYbin.pivot(index='binned_aafY', columns='binned_aafX', values='prop_per_binXY')
heatdata.index = [round(idx, 3) for idx in heatdata.index]  # weird trick for avoiding endless decimals, round does not work
heatdata.columns = [round(col, 3) for col in heatdata.columns]
heatdata = heatdata.sort_index(axis=0, ascending=False)
fxmin, fxmax = heatdata.columns[0], heatdata.columns[-1]
fymin, fymax = heatdata.index[-1], heatdata.index[0]
print('Max on y-axis: ', fymax)
print('Max on x-axis: ', fxmax)

if ydata2 is not None:
    heatdata2 = counts_per_XYbin_merged.pivot(index='binned_aafY', columns='binned_aafX', values='prop_per_binXY2')
    heatdata2.index = [round(idx, 3) for idx in
                      heatdata2.index]  # weird trick for avoiding endless decimals, round does not work
    heatdata2.columns = [round(col, 3) for col in heatdata2.columns]
    heatdata2 = heatdata2.sort_index(axis=0, ascending=False)
    fxmin2, fxmax2 = heatdata2.columns[0], heatdata2.columns[-1]
    fymin2, fymax2 = heatdata2.index[-1], heatdata2.index[0]

# Draw plot
sns.set(font_scale=2, style='white')  # increase font size of seaborn labels (overwritten then for the ticks on the axes)
sns.axes_style(style={'xtick.bottom': True,
                      'axes.spines.bottom': True,
                      'axes.facecolor': 'white',
                      'axes.edgecolor': 'black'})

fig, ax = plt.subplots(figsize=(16, 12))

sns.heatmap(heatdata, annot=False, alpha=0.7, ax=ax, xticklabels=5, yticklabels=5,
            cmap="plasma_r", cbar=True, cbar_kws={'shrink': .75,
                                                  'label': 'Imputed data'})  # viridis_r color works too
if ydata2 is not None:
    sns.heatmap(heatdata2, annot=False, alpha=0.7, ax=ax, xticklabels=5, yticklabels=5,
                cmap="Greys", cbar=True, cbar_kws={'shrink': .75,
                                                   'label': 'Pooled data'})

newax = fig.add_axes(ax.get_position(), frameon=False)
newax.patch.set_alpha(0.0)
newax.get_xaxis().set_visible(False)
newax.get_yaxis().set_visible(False)
sns.lineplot(x=[0.0, 0.5],
             y=[0.0, 0.5],
             color='green',
             linewidth=2,
             label='target',
             ax=newax)
newax.axis('equal')

ax.set_xlabel('AAF in study population', fontsize=20)
ax.set_ylabel('AAF in other population', fontsize=20)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=18)
ax.set_yticklabels(ax.get_yticklabels(), rotation=45, ha='right', fontsize=18)
cbar = ax.collections[0].colorbar
# cbar.ax.tick_params(labelsize=16)
# ax.figure.axes[-1].yaxis.label.set_size(20)
newax.legend(loc='upper left', edgecolor='white', fontsize=18)
plt.savefig(os.path.join(outdir, f'heatmap-STU-aaf-zoomin-targetline.pdf') if ydata2 is not None
            else os.path.join(outdir, f'heatmap-{pop}-aaf-zoomin-targetline.pdf'))
plt.show()
