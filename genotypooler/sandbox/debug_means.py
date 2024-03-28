import sys, os
import argparse
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

# force PYTHONPATH to look into the project directory for modules
rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from genotypooler.poolSNPs.metrics import quality as qual
from genotypooler.poolSNPs import dataframe as vcfdf
from genotypooler.persotools.files import *


datadir = '/home/camille/IterDecodeImpute/runs/poolimputeSNPs/results/data/1/'
datapaths = {
    'alltrue': 'STU.Chr1.SNPs.pruned.sorted.vcf.gz',
    'allimp': 'cycle2/STU.Chr1.SNPs.pruned.sorted.pooled.imputed.vcf.gz',
    'updated': 'cycle2/STU.Chr1.SNPs.pruned.sorted.pooled.changed_priors.vcf.gz'
}

dfmask = vcfdf.PandasMixedVCF(datadir + datapaths['updated'], format='GP', indextype='chrom:pos', mask=None)
gdata = dfmask.genotypes()
boolmask = gdata.applymap(lambda x: True if x[0] is None else False).values  # masks unchanged priors (which are None)

dftrueall = vcfdf.PandasMixedVCF(datadir + datapaths['alltrue'], format='GT', indextype='chrom:pos', mask=None)
dfimpall = vcfdf.PandasMixedVCF(datadir + datapaths['allimp'], format='GP', indextype='chrom:pos', mask=None)
dfimpupdat = vcfdf.PandasMixedVCF(datadir + datapaths['allimp'], format='GP', indextype='chrom:pos', mask=boolmask)
dfimpnochg = vcfdf.PandasMixedVCF(datadir + datapaths['allimp'], format='GP', indextype='chrom:pos', mask=~boolmask)

allgenos = dfimpall.genotypes()
updatgenos = dfimpupdat.genotypes()
unchggenos = dfimpnochg.genotypes()

notnaSall = allgenos.notna().sum(axis=1)
notnaSupdat = updatgenos.notna().sum(axis=1)
notnaSunchg = unchggenos.notna().sum(axis=1)

print(notnaSall, '\n')
print(notnaSupdat, '\n')
print(notnaSunchg, '\n')


print(f'''\nAll data retrieved from updated set and unchanged set? 
--> {notnaSall.all() == notnaSupdat.all() 
     + notnaSunchg.all()}\n''')

q1gt = qual.QualityGT(datadir + datapaths['alltrue'],
                      datadir + datapaths['allimp'],
                      0, idx='chrom:pos',
                      mask=None)
q1gl = qual.QualityGL(datadir + datapaths['alltrue'],
                      datadir + datapaths['allimp'],
                      0, idx='chrom:pos',
                      mask=None)

q2gt = qual.QualityGT(datadir + datapaths['alltrue'],
                      datadir + datapaths['allimp'],
                      0, idx='chrom:pos',
                      mask=boolmask)
q2gl = qual.QualityGL(datadir + datapaths['alltrue'],
                      datadir + datapaths['allimp'],
                      0, idx='chrom:pos',
                      mask=boolmask)

q3gt = qual.QualityGT(datadir + datapaths['alltrue'],
                      datadir + datapaths['allimp'],
                      0, idx='chrom:pos',
                      mask=~boolmask)
q3gl = qual.QualityGL(datadir + datapaths['alltrue'],
                      datadir + datapaths['allimp'],
                      0, idx='chrom:pos',
                      mask=~boolmask)

print('CONCORDANCE: all, updated data, unchanged data')
print(q1gt.concordance(), '\n')
print(q2gt.concordance(), '\n')
print(q3gt.concordance(), '\n')
