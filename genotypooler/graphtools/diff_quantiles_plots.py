"""
Plots comparing the accuracy of imputation between consecutive runs.
The following metrics are computed across MAF range and shown with quantiles [25%, 50%, 75%] dispersion.
* Concordance: marker-wise difference cycleX - cycleX+1
* Cross-entropy: marker-wise difference cycleX - cycleX+1

Quantiles computed on rolling windows.

Arguments are parsed from a file where they must be written in this order:
pathout <>
truegt <>
truegl <>
imp1 <path/to/file/imputed/with/data1>
imp2 <path/to/file/imputed/with/data2>
date <>
rollwin <>
bins <>
compute <>

Beware that strings for paths should be written just as text (without quotes!) in the argsfile!

Command line usage (assuming the current directory is genotypooler/graphtools)
$ python3 -u diff_quantiles_plots.py @/home/camille/IterDecodeImpute/runs/poolimputeSNPs/results/data/1/cycle2/quantilesargsfile.txt
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
imputed_1 = argsin.imp1
imputed_2 = argsin.imp2
datedir = argsin.date
rQ = argsin.rollwin * 4  # try some integers in the range [2, 5] and see what is the smoothest
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

def diff_rollquants(dX: pd.DataFrame, dS1: pd.Series) -> pd.DataFrame: # , dS2: pd.Series
    pdf1 = qual.QuantilesDataFrame(dX,
                                   dS1,
                                   bins_step=bS)
    pctY1 = pdf1.binnedX_rolling_quantilY(rollwin=rQ)
    pctY1['dataset'] = ['cycleX - cycleX+1'] * pctY1.shape[0]
    diff_rollquants = pctY1

    # pdf2 = qual.QuantilesDataFrame(dX,
    #                                dS2,
    #                                bins_step=bS)
    # pctY2 = pdf2.binnedX_rolling_quantilY(rollwin=rQ)
    # pctY2['dataset'] = ['2'] * pctY2.shape[0]
    #
    # diff_rollquants = pd.concat([pctY1, pctY2])

    return diff_rollquants


# Load data and check

q1gt = qual.QualityGT(truegt, imputed_1, 0, idx='chrom:pos')
q1gl = qual.QualityGL(truegl, imputed_1, 0, idx='chrom:pos')

print('\r\n{} variants from {} samples read from {}'.format(len(q1gt.trueobj.variants),
                                                            len(q1gt.trueobj.samples),
                                                            os.path.basename(truegt)))
print('\r\n{} variants from {} samples read from {}'.format(len(q1gt.imputedobj.variants),
                                                            len(q1gt.imputedobj.samples),
                                                            os.path.basename(imputed_1)))
if compute:
    bgldiff = q1gt.diff()

q2gt = qual.QualityGT(truegt, imputed_2, 0, idx='chrom:pos')
q2gl = qual.QualityGL(truegl, imputed_2, 0, idx='chrom:pos')

print('\r\n{} variants from {} samples read from {}'.format(len(q1gl.trueobj.variants),
                                                            len(q1gl.trueobj.samples),
                                                            os.path.basename(truegl)))
print('\r\n{} variants from {} samples read from {}'.format(len(q2gl.imputedobj.variants),
                                                            len(q2gl.imputedobj.samples),
                                                            os.path.basename(imputed_2)))

mafS = q1gt.trueobj.maf

if compute:
    metrics = {
        'concordance': {'1': q1gt.marker_wise_concordance(),
                        '2': q2gt.marker_wise_concordance()},
        'cross_entropy': {'1': q1gl.marker_wise_cross_entropy,
                          '2': q2gl.marker_wise_cross_entropy}
    }

dataquants = {
    'concordance': os.path.join(outdir, 'diff_rolling_quantiles_concordance.json'),
    'cross_entropy': os.path.join(outdir, 'diff_rolling_quantiles_cross_entropy.json')
}

# Process and write data

if compute:
    for metric, d in metrics.items():
        if d is not None:
            yDF_1, yDF_2 = list(d.values())
            yS1_2 = yDF_1.subtract(yDF_2).mean(axis=1).rename(metric)

            # Compute quantiles
            print('Computing quantiles for {}'.format(metric).ljust(80, '.'))
            pctY_comp = diff_rollquants(mafS, yS1_2)

            # Compute mean over all markers
            print('Computing means for {}'.format(metric).ljust(80, '.'))
            pctY_comp['mean'] = pctY_comp['dataset'].apply(lambda x: yS1_2.mean())
            jsonf = dataquants[metric]
            pctY_comp.to_json(jsonf,
                              orient='records')
if True:
    # Read processed reshaped data for plotting and draw figures

    sns.set(font_scale=1.75)  # multiplication factor!

    for dquant, f in dataquants.items():
        if f is not None:
            dataf = pd.read_json(f, orient='records')
            print(dataf)
            meanf = {}

            gY = sns.lineplot(data=dataf[dataf.quantiles == 0.5], x=x_data, y=dquant,
                              hue='dataset', palette="husl", linewidth=1)
            for i, dset in enumerate(['cycleX - cycleX+1']):
                df = dataf[dataf['dataset'] == dset]
                meanf[dset] = df['mean'].mean()
                gY.fill_between(df[df.quantiles == 1.0][x_data],
                                df[df.quantiles == 0.0][dquant],
                                df[df.quantiles == 1.0][dquant],
                                color=sns.color_palette('husl')[i],
                                alpha=0.1)
                gY.fill_between(df[df.quantiles == 0.99][x_data],
                                df[df.quantiles == 0.01][dquant],
                                df[df.quantiles == 0.99][dquant],
                                color=sns.color_palette('husl')[i],
                                alpha=0.25)
                gY.fill_between(df[df.quantiles == 0.75][x_data],
                                df[df.quantiles == 0.25][dquant],
                                df[df.quantiles == 0.75][dquant],
                                color=sns.color_palette('husl')[i],
                                alpha=0.40)
            gY.set_xlabel('True minor allele frequency in {} population'.format('study' if x_data == 'binned_maf'
                                                                                else 'main'),
                          fontsize=axlabsz)
            gY.set_ylabel(str.capitalize(dataf.columns[2].replace('_', ' ')), fontsize=axlabsz)

            # gY.set(ylim=yscale[dquant])
            # if dquant == 'cross_entropy':
            #     gY.set(yscale="log")
            handles, labels = gY.get_legend_handles_labels()
            # labels[-2] = '{} (mean = {:.5f})'.format(labels[-2], meanf['1'])
            labels[-1] = '{} (mean = {:.5f})'.format(labels[-1], meanf['cycleX - cycleX+1'])
            gY.legend(handles, labels, loc='lower right' if dquant == 'concordance' else 'upper right', fontsize=legsz)
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, 'diff_{}_percentiles_rQ={}_bS={}_xdata={}.pdf'.format(dquant, rQ, bS, x_data.lstrip('binned_'))))
            # plt.show()
            plt.close()