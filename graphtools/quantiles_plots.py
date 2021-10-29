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
imp2 <path/to/file/imputed/with/Phaser>
date <>
rollwin <>
bins <>
compute <>

Beware that strings for paths should be written just as text (without quotes!) in the argsfile!

Command line usage (assuming the current directory is genotypooler/manus)
$ python3 -u quantiles_plots.py @/home/camcl609/PoolImpHuman/results/20210320/argsfile20210320.txt
"""

import os, sys
import numpy as np
import pandas as pd
import json
import seaborn as sns
import matplotlib.pyplot as plt
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
parser.add_argument('imp2', metavar='imp2', type=str, help='File with imputed data 2 (GT and/or GP formats)', default=None)
parser.add_argument('date', metavar='date', type=str, help='Date of the experiment (YYYMMDD)', default=None)
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
imputed_phaser = argsin.imp2
datedir = argsin.date
rQ = argsin.rollwin
bS = argsin.bins
compute = argsin.compute

# Data parameters

x_data = 'binned_maf'  # 'binned_maf_info'
# rQ = 1000
# bS = 0.01
# compute = False

# Configure data/plots paths

outdir = os.path.join(pathout, datedir)
if not os.path.exists(outdir):
    os.mkdir(outdir)


print('\r\nData written to {}'.format(outdir))

# Plot styling

# General parameters
plt.style.use('./manus-style.mplstyle')

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


# Function/Tools

def rollquants(dX: pd.DataFrame, dS1: pd.Series, dS2: pd.Series) -> pd.DataFrame:
    pdf1 = qual.QuantilesDataFrame(dX,
                                   dS1,
                                   bins_step=bS)
    pctY1 = pdf1.binnedX_rolling_quantilY(rollwin=rQ)
    pctY1['dataset'] = ['beagle'] * pctY1.shape[0]

    pdf2 = qual.QuantilesDataFrame(dX,
                                   dS2,
                                   bins_step=bS)
    pctY2 = pdf2.binnedX_rolling_quantilY(rollwin=rQ)
    pctY2['dataset'] = ['phaser'] * pctY2.shape[0]

    rollquants = pd.concat([pctY1, pctY2])

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

qphasergt = qual.QualityGT(truegt, imputed_phaser, 0, idx='chrom:pos')
qphasergl = qual.QualityGL(truegl, imputed_phaser, 0, idx='chrom:pos')

print('\r\n{} variants from {} samples read from {}'.format(len(qbeaglegl.trueobj.variants),
                                                            len(qbeaglegl.trueobj.samples),
                                                            os.path.basename(truegl)))
print('\r\n{} variants from {} samples read from {}'.format(len(qphasergl.imputedobj.variants),
                                                            len(qphasergl.imputedobj.samples),
                                                            os.path.basename(imputed_phaser)))

mafS = qbeaglegt.trueobj.maf

if compute:
    metrics = {
        'concordance': {'beagle': qbeaglegt.concordance(),
                        'phaser': qphasergt.concordance()},
        'cross_entropy': {'beagle': qbeaglegl.cross_entropy,
                          'phaser': qphasergl.cross_entropy}
    }

dataquants = {
    'concordance': os.path.join(outdir, 'rolling_quantiles_concordance.json'),
    'cross_entropy': os.path.join(outdir, 'rolling_quantiles_cross_entropy.json')
}

# Process and write data

if compute:
    for metric, d in metrics.items():
        if d is not None:
            yS_beagle, yS_phaser = list(d.values())
            # Compute quantiles
            print('Computing quantiles for {}'.format(metric).ljust(80, '.'))
            pctY_comp = rollquants(mafS, yS_beagle, yS_phaser)
            # Compute mean over all markers
            print('Computing means for {}'.format(metric).ljust(80, '.'))
            pctY_comp['mean'] = pctY_comp['dataset'].apply(lambda x: yS_beagle.mean() if x == 'beagle' else yS_phaser.mean())
            jsonf = dataquants[metric]
            pctY_comp.to_json(jsonf,
                              orient='records')


# Read processed reshaped data for plotting and draw figures

sns.set(font_scale=1.75)  # multiplication factor!

for dquant, f in dataquants.items():
    if f is not None:
        dataf = pd.read_json(f, orient='records')
        print(dataf)
        meanf = {}

        gY = sns.lineplot(data=dataf[dataf.quantiles == 0.5], x=x_data, y=dquant,
                          hue='dataset', palette="deep", linewidth=1)
        for i, dset in enumerate(['beagle', 'phaser']):
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
        # if dquant == 'cross_entropy':
        #     gY.set(yscale="log")
        handles, labels = gY.get_legend_handles_labels()
        labels[-2] = '{} (mean = {:.5f})'.format(labels[-2], meanf['beagle'])
        labels[-1] = '{} (mean = {:.5f})'.format(labels[-1], meanf['phaser'])
        gY.legend(handles, labels, loc='lower right' if dquant == 'concordance' else 'upper right', fontsize=legsz)
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, '{}_percentiles_rQ={}_bS={}_xdata={}.pdf'.format(dquant, rQ, bS, x_data.lstrip('binned_'))))
        # plt.show()
        plt.close()
