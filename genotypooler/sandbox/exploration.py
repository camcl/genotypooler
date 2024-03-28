"""
Script for exploring the genotype data from the MagicWheat data set.
VCF files are created with PLINK.
"""

import numpy as np
import pandas as pd
from sklearn.metrics import multilabel_confusion_matrix, confusion_matrix
from sklearn.preprocessing import *
from scipy.stats import *
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, to_rgba_array, to_rgba
import seaborn as sns
import os, sys
import argparse

rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from genotypooler.poolSNPs import dataframe as vcfdf

BASIC_GWAS = False
FOUNDERS = True

if BASIC_GWAS:
    truef = '/home/camille/MagicWheat/data/BASIC_GWAS/DATA/MAGIC_IMPUTED_PRUNED/MAGIC_imputed.pruned.vcf.gz'
    print(f"\nParsing file {truef}")

    dftrue = vcfdf.PandasMixedVCF(truef, format='GT')

    # truegt = dftrue.genotypes()

    # trueaaf = dftrue.aaf

    truemiss = dftrue.missing_rate

    print(truemiss)
    print(f"\nStatistics for missing rate:")
    print(truemiss.describe())
    print(f"Number of missing markers:{len(truemiss[truemiss['missing_rate'] == 0.0])}")

    truehet = dftrue.het_rate

    print(truehet)
    print(f"\nStatistics for heterozygosity rate:")
    print(truehet.describe())

    # What if I remove SNPs with missing data?

if FOUNDERS:
    foundersf = '/home/camille/MagicWheat/data/FOUNDERS/1/Chr1.Founders.vcf.gz'
    print(f"\nParsing file {foundersf}")
    inbredf = '/home/camille/MagicWheat/data/20220611/1/Chr1.SNPs.pruned.nomiss.vcf.gz'
    print(f"\nParsing file {inbredf}")

    dffounders = vcfdf.PandasMixedVCF(foundersf, format='GT')
    dfinbred = vcfdf.PandasMixedVCF(inbredf, format='GT')

    print(dffounders.pos)
    print(dfinbred.pos.values.flatten())

    dfsharedpos = dffounders.pos.join(dfinbred.pos, how='inner', lsuffix='_founders', rsuffix='_inbred')
    print(dfsharedpos)

    plt.rcParams["figure.figsize"] = [16, 4]
    fig, axes = plt.subplots(nrows=1, ncols=1)#, sharex=True)
    axes.scatter(dfinbred.pos.values.flatten(), np.ones_like(dfinbred.pos.values.flatten()), s=1, c='g')
    axes.scatter(dffounders.pos.values.flatten(), np.zeros_like(dffounders.pos.values.flatten()), s=1, c='b')
    plt.show()

