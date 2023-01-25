import os, sys
import numpy as np
import pandas as pd
import seaborn as sns
import timeit
import multiprocessing as mp
import matplotlib.pyplot as plt

from typing import *

rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from genotypooler.poolSNPs.metrics import quality as qual
from genotypooler.poolSNPs import dataframe as vcfdf
from genotypooler.persotools.files import *

'''
This script addresses the following review for the manuscript:
"How many more markers are correctly imputed in pooled data in addition to the ones that are fully decoded?"
'''

compute = True

# Plot styling

# General parameters
plt.style.use('/home/camille/1000Genomes/src/genotypooler/manus/manus-style.mplstyle')

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
dash_styles = [
    (1, 1),
    (3, 1, 1.5, 1),
    (5, 1, 1, 1),
    (5, 1, 2, 1, 2, 1),
    (2, 2, 3, 1.5),
    (1, 2.5, 3, 1.2),
    "",
    (4, 1.5),
]

# Data parameters

bins_step = 0.01
rQ = 1000
x_data = 'maf_info'
x_bins = np.arange(0.0, 0.5 + bins_step, bins_step) if x_data in ['maf_info', 'maf'] \
            else np.arange(0.0, 1.0 + bins_step, bins_step)
# Custom bins for exact matches counts
x2_bins = [0.0, 0.02, 0.04, 0.06, 0.1, 0.2, 0.4, 0.5]  # MAF
lab2_bins = [0.01, 0.03, 0.05, 0.08, 0.15, 0.3, 0.45]  # MAF
lab2_fmt = ['{:.2f}-{:.2f}'.format(i, j) for i, j in zip(x2_bins[:-1], x2_bins[1:])]


# Configure data/plots paths

datedir = '20200827'
outdir = os.path.join('/home/camille/PoolImpHuman/results', datedir)
if not os.path.exists(outdir):
    os.mkdir(outdir)
    
    
# Function/Tools

def rollquants(dX: pd.DataFrame, dS1: pd.Series, dS2: pd.Series) -> pd.DataFrame:
    pdf1 = qual.QuantilesDataFrame(dX,
                                   dS1,
                                   bins_step=bins_step)
    pctY1 = pdf1.binnedX_rolling_quantilY(rollwin=rQ)
    pctY1['dataset'] = ['beagle'] * pctY1.shape[0]
    print(pctY1)

    pdf2 = qual.QuantilesDataFrame(dX,
                                   dS2,
                                   bins_step=bins_step)
    pctY2 = pdf2.binnedX_rolling_quantilY(rollwin=rQ)
    pctY2['dataset'] = ['phaser'] * pctY2.shape[0]
    print(pctY2)

    rollquants = pd.concat([pctY1, pctY2])

    return rollquants


# Coordinates of HDonly and LDonly variants chr20 x Illumina (35,682 and 17,015 variants)

ld_vars = pd.read_csv('/home/camille/PoolImpHuman/data/omniexpress-isec-LDHD/LDonly.coords',
                      sep='\t',
                      header=None,
                      names=['CHROM', 'POS'])
ld_vars['variants'] = pd.Series(['20:{}'.format(pos) for pos in ld_vars['POS']], dtype=str)
ld_vars.set_index('variants', inplace=True)

hd_vars = pd.read_csv('/home/camille/PoolImpHuman/data/omniexpress-isec-LDHD/HDonly.coords',
                      sep='\t',
                      header=None,
                      names=['CHROM', 'POS'])
hd_vars['variants'] = pd.Series(['20:{}'.format(pos) for pos in hd_vars['POS']], dtype=str)
hd_vars.set_index('variants', inplace=True)


# Read files

truegt = '/home/camille/PoolImpHuman/data/20200827/IMP.chr20.snps.gt.vcf.gz'
truegl = '/home/camille/PoolImpHuman/data/20200827/IMP.chr20.snps.gl.vcf.gz'
# pooled can also be the file with full LD and missing HD
pooledgt = None  # '/home/camille/PoolImpHuman/data/20201029/sHG01063.IMP.chr20.missingHD.fullLD.snps.gt.vcf.gz'
pooledgl = '/home/camille/PoolImpHuman/data/20210607/IMP.chr20.pooled.snps.gl.vcf.gz'  # Caution: older files might have min(log-GL) = -5 and not -12
# imputation with Beagle or with Phaser
imputed_gtgp1 = '/home/camille/PoolImpHuman/data/20200827/IMP.chr20.pooled.imputed.vcf.gz'  # Phaser
imputed_gtgp2 = '/home/camille/PoolImpHuman/data/20210607/IMP.chr20.pooled.imputed.vcf.gz'  # Beagle


# Build Quality and DataFrame objects for analysis

startT = timeit.default_timer()

quality1gt = qual.QualityGT(truegt, imputed_gtgp1, 0, idx='chrom:pos')
quality1gl = qual.QualityGL(truegl, imputed_gtgp1, 0, idx='chrom:pos')

if imputed_gtgp2 is not None:
    quality2gt = qual.QualityGT(truegt, imputed_gtgp2, 0, idx='chrom:pos')
    quality2gl = qual.QualityGL(truegl, imputed_gtgp2, 0, idx='chrom:pos')

if pooledgt is not None:
    dfpooled = vcfdf.PandasMixedVCF(pooledgt, format='GT', indextype='chrom:pos')

if pooledgl is not None:
    dfpooled = vcfdf.PandasMixedVCF(pooledgl, format='GL', indextype='chrom:pos')


if compute:

    # Filter/mask markers on mismatches and decoding status
    disc1 = quality1gt.diff().dropna()
    miss1 = dfpooled.trinary_encoding()
    entro1 = quality1gl.trueobj.genotypes().combine(quality1gl.imputedobj.genotypes(), quality1gl.intergl_entropy)

    disc2 = quality2gt.diff().dropna()  # Non NaN i.e. everything imputed
    miss2 = dfpooled.trinary_encoding()
    entro2 = quality2gl.trueobj.genotypes().combine(quality2gl.imputedobj.genotypes(), quality2gl.intergl_entropy)

    unassayed_disc1 = 0.5 * disc1.where(miss1 == -1, other=np.nan)  # keep (mis)matches of only not decoded markers
    unassayed_entro1 = entro1.where(miss1 == -1, other=np.nan)

    unassayed_disc2 = 0.5 * disc2.where(miss2 == -1, other=np.nan)  # discordance is 0, 1, or 2
    unassayed_entro2 = entro2.where(miss2 == -1, other=np.nan)

    # Additional exact matches

    # Count not decoded and correctly imputed
    unassayed_impok1 = unassayed_disc1.where(unassayed_disc1 == 0.0, other=np.nan)  # keep not decoded and full matches after imputation
    counts_unassayed_impok1 = unassayed_impok1.notnull().sum(axis=1)
    counts_unassayed_impok1.name = 'counts'

    unassayed_impok2 = unassayed_disc2.where(unassayed_disc2 == 0.0, other=np.nan)  # keep not decoded and full matches after imputation
    counts_unassayed_impok2 = unassayed_impok2.notnull().sum(axis=1)
    counts_unassayed_impok2.name = 'counts'

    # Count not decoded per variant
    not_decoded1 = miss1.where(miss1 == -1, other=np.nan)
    counts_not_decoded1 = not_decoded1.notnull().sum(axis=1)

    not_decoded2 = miss2.where(miss2 == -1, other=np.nan)
    counts_not_decoded2 = not_decoded2.notnull().sum(axis=1)

    # TODO: exact matches counts per bin
    # Proportion of full matches after imputation for not decoded markers
    dX = quality1gt.trueobj.maf_info
    matchbins = pd.cut(dX.values.squeeze(), x2_bins, labels=lab2_fmt, include_lowest=True)  # lab2_bins
    dBinsX = pd.DataFrame(matchbins, index=dX.index, columns=['binned_' + x_data])

    sExact1 = counts_unassayed_impok1 / counts_not_decoded1
    sExact1.name = 'exact_matches'  # counts exactly imputed only
    # exact_matches1 = qual.QuantilesDataFrame(dX, sExact1)
    dfmatch1 = dBinsX.join(sExact1).groupby(by='binned_maf_info').mean()
    dfmatch1['dataset'] = 'phaser'

    sExact2 = counts_unassayed_impok2 / counts_not_decoded2
    sExact2.name = 'exact_matches'  # counts exactly imputed only
    # exact_matches2 = qual.QuantilesDataFrame(dX, sExact2)
    dfmatch2 =  dBinsX.join(sExact2).groupby(by='binned_maf_info').mean()
    dfmatch2['dataset'] = 'beagle'

    dfmatches = pd.concat([dfmatch1, dfmatch2])

    tab1 = dfmatch1.reset_index().pivot(index='dataset', columns='binned_maf_info', values='exact_matches')
    tab2 = dfmatch2.reset_index().pivot(index='dataset', columns='binned_maf_info', values='exact_matches')
    tabmatches = pd.concat([tab1, tab2])
    tabmatches.index.name = ''
    tabmatches.index = [algo.capitalize() for algo in tabmatches.index]
    tabmatches.columns.name = ''
    tabmatches.to_latex(buf=os.path.join(outdir, 'pooled-not-decoded-exact-matches-bin.tex'),
                        sparsify=True,
                        multirow=True,
                        caption='Exact matches after imputation for not decoded genotypes in pooled HD per data {} bin: Columns labels are the central values of each interval'.format(
                            'MAF' if (x_data == 'maf' or x_data == 'maf_info') else 'AAF'),
                        label='tab:pooled-not-decoded-exact-matches-bin')

    # Continuous quality metrics with quantiles

    # Concordance after imputation for not decoded markers
    unassayed_concordance1 = (1 - unassayed_disc1.mean(axis=1))
    unassayed_concordance1.dropna(inplace=True)
    unassayed_concordance1.name = 'concordance'

    unassayed_concordance2 = (1 - unassayed_disc2.mean(axis=1))
    unassayed_concordance2.dropna(inplace=True)
    unassayed_concordance2.name = 'concordance'

    # Cross-entropy after imputation for not decoded markers
    unassayed_crossentro1 = unassayed_entro1.mean(axis=1)
    unassayed_crossentro1.dropna(inplace=True)
    unassayed_crossentro1.name = 'cross_entropy'

    unassayed_crossentro2 = unassayed_entro2.mean(axis=1)
    unassayed_crossentro2.dropna(inplace=True)
    unassayed_crossentro2.name = 'cross_entropy'

    # unassayed_concordance2.notnull().sum()
    # Out[99]: 51395
    # unassayed_concordance1.notnull().sum()
    # Out[100]: 51395

    # Metrics to compute and save

    dX = quality1gt.trueobj.maf_info
    qualDF11 = qual.QuantilesDataFrame(dX, unassayed_concordance1)
    qualDF12 = qual.QuantilesDataFrame(dX, unassayed_crossentro1)

    print('Example of quantiles for concordance:')
    print(qualDF11.binnedX_rolling_quantilY())
    print('Example of quantiles for cross-entropy:')
    print(qualDF12.binnedX_rolling_quantilY())

    qualDF21 = qual.QuantilesDataFrame(dX, unassayed_concordance2)


    metrics = {
        'concordance': {'beagle': unassayed_concordance2,
                        'phaser': unassayed_concordance1},
        'cross_entropy': {'beagle': unassayed_crossentro2,
                          'phaser': unassayed_crossentro1}
    }

dataquants = {
              'concordance': os.path.join(outdir, 'rolling_quantiles_not_decoded_concordance.json'),
              'cross_entropy': os.path.join(outdir, 'rolling_quantiles_not_decoded_cross_entropy.json')
              }

if compute:

    dX = quality1gt.trueobj.maf_info
    # Save as table

    for metric, d in metrics.items():
        break
        if d is not None:
            yS_beagle, yS_phaser = list(d.values())
            # Compute quantiles
            print('Computing quantiles for {}'.format(metric).ljust(80, '.'))
            pctY_comp = rollquants(dX, yS_beagle, yS_phaser)
            # Compute mean over all markers
            print('Computing means for {}'.format(metric).ljust(80, '.'))
            pctY_comp['mean'] = pctY_comp['dataset'].apply(lambda x: yS_beagle.mean() if x == 'beagle' else yS_phaser.mean())
            jsonf = dataquants[metric]
            pctY_comp.to_json(jsonf,
                              orient='records')

# jsonf = os.path.join(outdir, 'exact_matches-not-decoded.json')
# jsonfroll = os.path.join(outdir, 'exact_matches-not-decoded-rolling.json')
#
# df_agg.to_json(jsonf, orient='records')
# dataf = pd.read_json(jsonf, orient='records')
#
# df_roll.to_json(jsonfroll, orient='records')
# rolldataf = pd.read_json(jsonfroll, orient='records')

# Plot
# dataf.plot()
# plt.show()

stopT = timeit.default_timer()
print('Time elapsed for computing and building data to plot *= {}'.format(stopT-startT).ljust(80, '.'))

# Read processed reshaped data for plotting and draw figures

sns.set(font_scale=1.75)  # multiplication factor!

for dquant, f in dataquants.items():
    break
    if f is not None:
        dataf = pd.read_json(f, orient='records')
        meanf = {}

        gY = sns.lineplot(data=dataf[dataf.quantiles == 0.5], x='binned_' + x_data, y=dquant,
                          hue='dataset', palette="deep", linewidth=1)
        for i, dset in enumerate(['beagle', 'phaser']):
            df = dataf[dataf['dataset'] == dset]
            meanf[dset] = df['mean'].mean()
            gY.fill_between(df[df.quantiles == 1.0]['binned_' + x_data],
                            df[df.quantiles == 0.0][dquant],
                            df[df.quantiles == 1.0][dquant],
                            color=sns.color_palette('deep')[i],
                            alpha=0.1)
            gY.fill_between(df[df.quantiles == 0.99]['binned_' + x_data],
                            df[df.quantiles == 0.01][dquant],
                            df[df.quantiles == 0.99][dquant],
                            color=sns.color_palette('deep')[i],
                            alpha=0.25)
            gY.fill_between(df[df.quantiles == 0.75]['binned_' + x_data],
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
        labels[-2] = '{} (mean = {:.5f})'.format(labels[-2], meanf['beagle'])
        labels[-1] = '{} (mean = {:.5f})'.format(labels[-1], meanf['phaser'])
        gY.legend(handles, labels, fontsize=legsz)
        plt.savefig(os.path.join(outdir, '{}_percentiles_rQ={}_bS={}_xdata={}_not_decoded.pdf'.format(dquant, rQ, bins_step, x_data.lstrip('binned_'))))
        plt.show()


