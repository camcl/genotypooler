"""
Plots quantifying the imputation performance.
The following metrics are computed across MAF range and shown with quantiles [25%, 50%, 75%] dispersion.
* Concordance:
* Cross-entropy:

Quantiles computed on rolling windows.

Arguments are parsed from a file where they must be written in this order:
pathout <>
truegt <>
truegl <>
imp1 <path/to/file/imputed/with/Beagle>
date <>
rollwin <>
bins <>
compute <>

Beware that strings for paths should be written just as text (without quotes!) in the argsfile!

Command line usage (assuming the current directory is genotypooler/examples)
$ python3 -u quantiles_1plot.py @argsfile_example.txt
"""

import os, sys
import collections
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter, FixedLocator
import seaborn as sns
import argparse

import warnings
warnings.filterwarnings("ignore")

rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from genotypooler.poolSNPs.metrics import quality as qual
from genotypooler.persotools.files import *


### COMMAND-LINE PARSING AND PARAMETERS (arguments are parsed from a file
parser = argparse.ArgumentParser(description='Plots the imputation performance'
                                             '(concordance and cross-entropy)',
                                 fromfile_prefix_chars='@')
parser.add_argument('pathout', metavar='out', type=str, help='Results directory', default=None)
parser.add_argument('truegt', metavar='gt', type=str, help='File with true data (GT format)', default=None)
parser.add_argument('truegl', metavar='gl', type=str, help='File with true data (GL format)', default=None)
parser.add_argument('imp1', metavar='imp1', type=str, help='File with imputed data 1 (GT and/or GP formats)', default=None)
parser.add_argument('date', metavar='date', type=str, help='Date of the experiment (YYYMMDD or today)', default=None)
parser.add_argument('rollwin', metavar='wq', type=int, help='Number of markers per rolling window', default=1000)  # default option does not work
parser.add_argument('bins', metavar='bin', type=float, help='Bin size for discretizing MAF', default=0.01)  # default option does not work
parser.add_argument('compute', metavar='comp', type=int, help='If True, compute quantiles and plots, else runs plotting only', default=1)


argsin = parser.parse_args()
print('\n'.ljust(80, '*'))
print('The following arguments were parsed from file:\n')
print(argsin)
print('\n'.ljust(80, '*'))

pathout = argsin.pathout
truegt = argsin.truegt
truegl = argsin.truegl
imputed_beagle = argsin.imp1
datedir = argsin.date
rQ = argsin.rollwin
bS = argsin.bins
compute = argsin.compute

# Data parameters

x_data = 'binned_maf'

# Configure data/plots paths

outdir = os.path.join(pathout, datedir)
if not os.path.exists(outdir):
    os.makedirs(outdir)  # recursive directories creation (mkdir is not recursive)


print('\r\nData written to {}'.format(outdir))

# Plot styling

# General parameters
# stylesheet

# Specific to this plotting script
sns.set(rc={'figure.figsize': (10, 8)})  # specific to this plotting sripts
sns.set_style('whitegrid')
titlesz = 24
axlabsz= 20
axticksz = 16
legsz = 20
yscale = {
    'concordance': (0.0, 1.0),
    'cross_entropy': (0.0, 12.0)
}
x_bins = [0.0, 0.02, 0.04, 0.06, 0.1, 0.2, 0.4, 0.5]
lab_bins = [0.01, 0.03, 0.05, 0.08, 0.15, 0.3, 0.45]
lab_fmt = ['{:.2f}-{:.2f}'.format(i, j) for i, j in zip(x_bins[:-1], x_bins[1:])]

# Function/Tools

def rollquants(dX: pd.DataFrame, dS1: pd.Series) -> pd.DataFrame:
    pdf1 = qual.QuantilesDataFrame(dX,
                                   dS1,
                                   bins_step=bS)
    pctY1 = pdf1.binnedX_rolling_quantilY(rollwin=rQ)
    pctY1['dataset'] = ['beagle'] * pctY1.shape[0]

    rollquants = pctY1  # pd.concat([pctY1, pctY2])

    return rollquants


# Load data and check

qbeaglegt = qual.QualityGT(truegt, imputed_beagle, 0, idx='chrom:pos')
qbeaglegl = qual.QualityGL(truegl, imputed_beagle, 0, idx='chrom:pos')

print('\r\n{} variants from {} samples read from {}'.format(len(qbeaglegt.trueobj.variants),
                                                            len(qbeaglegt.trueobj.samples),
                                                            os.path.basename(truegt)))
print('\r\n{} variants from {} samples read from {}'.format(len(qbeaglegt.imputedobj.variants),
                                                            len(qbeaglegt.imputedobj.samples),
                                                            os.path.basename(imputed_beagle)))
if compute:
    bgldiff = qbeaglegt.diff()


print('\r\n{} variants from {} samples read from {}'.format(len(qbeaglegl.trueobj.variants),
                                                            len(qbeaglegl.trueobj.samples),
                                                            os.path.basename(truegl)))

mafS = qbeaglegt.trueobj.maf  # maf_info

# Create bins for barplot

binS = pd.cut(mafS.values.squeeze(), x_bins, labels=lab_fmt, include_lowest=True)
binDF = pd.DataFrame(data=binS, index=mafS.index, columns=['maf_bin']).reset_index()
countDict = collections.Counter(binDF['maf_bin'])
binDF['bin_counts'] = binDF['maf_bin'].apply(lambda x: countDict[x])
print(binDF)
print(countDict)

# Compute metrics

if compute:
    metrics = {'precision_score': qbeaglegt.precision,
               'recall_score': qbeaglegt.recall,
               'f1_score':  qbeaglegt.f1_score,
               'concordance': qbeaglegt.concordance(),
               'allelic_dos': None,
               'cross_entropy': qbeaglegl.cross_entropy
               }

dataquants = {'precision_score': None,  # os.path.join(outdir, 'rolling_quantiles_precision_score.json'),
              'recall_score': None,  # os.path.join(outdir, 'rolling_quantiles_recall_score.json'),
              'f1_score': None,  # os.path.join(outdir, 'rolling_quantiles_f1_score.json'),
              'concordance': os.path.join(outdir, 'rolling_quantiles_concordance.json'),
              'allelic_dos': None,
              'cross_entropy': os.path.join(outdir, 'rolling_quantiles_cross_entropy.json')
              }

# Process and write data

if compute:
    for metric, d in metrics.items():
        if d is not None:
            yS_beagle = d
            print(yS_beagle)
            # Compute quantiles
            print('Computing quantiles for {}'.format(metric).ljust(80, '.'))
            pctY_comp = rollquants(mafS, yS_beagle)
            # Compute mean over all markers
            print('Computing means for {}'.format(metric).ljust(80, '.'))
            pctY_comp['mean'] = yS_beagle.mean()
            jsonf = dataquants[metric]
            pctY_comp.to_json(jsonf,
                              orient='records')


# Read processed reshaped data for plotting and draw figures

sns.set(font_scale=1.75)  # multiplication factor!

# Histogram of markers count per MAF-bin

plt.bar(x=binDF['maf_bin'], height=binDF['bin_counts'].astype(int, copy=False), color='silver', alpha=0.5)

plt.xticks(rotation=45)
plt.xlabel('True minor allele frequency in {} population'.format('study' if x_data == 'binned_maf'
                                                                      else 'main'),
                fontsize=axlabsz)
plt.ylabel('Counts', fontsize=axlabsz)
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'histogram-counts.pdf'))
plt.show()
# plt.close()

# Dispersion of metrics

for dquant, f in dataquants.items():
    # break
    if f is not None:
        dataf = pd.read_json(f, orient='records')
        print(dataf)
        meanf = {}

        gY = sns.lineplot(data=dataf[dataf.quantiles == 0.5], x=x_data, y=dquant,
                          hue='dataset', palette="deep", linewidth=1)

        for i, dset in enumerate(['beagle']):
            df = dataf[dataf['dataset'] == dset]
            meanf[dset] = df['mean'].mean()
            gY.fill_between(df[df.quantiles == 1.0][x_data],
                            df[df.quantiles == 0.0][dquant],
                            df[df.quantiles == 1.0][dquant],
                            color=sns.color_palette('deep')[i],
                            alpha=0.1)
            gY.fill_between(df[df.quantiles == 0.99][x_data],
                            df[df.quantiles == 0.01][dquant],
                            df[df.quantiles == 0.99][dquant],
                            color=sns.color_palette('deep')[i],
                            alpha=0.25)
            gY.fill_between(df[df.quantiles == 0.75][x_data],
                            df[df.quantiles == 0.25][dquant],
                            df[df.quantiles == 0.75][dquant],
                            color=sns.color_palette('deep')[i],
                            alpha=0.40)

        gY.set_xlabel('True minor allele frequency in {} population'.format('study' if x_data == 'binned_maf'
                                                                            else 'main'),
                      fontsize=axlabsz)
        gY.set_ylabel(str.capitalize(dataf.columns[2].replace('_', ' ')), fontsize=axlabsz)

        gY.set(ylim=yscale[dquant])
        handles, labels = gY.get_legend_handles_labels()
        labels[-1] = '{} (mean = {:.5f})'.format(labels[-1], meanf['beagle'])
        gY.legend(handles, labels, loc='best', fontsize=legsz)
        plt.savefig(os.path.join(outdir, '{}_percentiles_rQ={}_bS={}_xdata={}.pdf'.format(dquant, rQ, bS, x_data.lstrip('binned_'))))
        plt.show()
