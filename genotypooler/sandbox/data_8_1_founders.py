import os, sys
import pandas as pd

rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from genotypooler.poolSNPs import dataframe as vcfdf
from genotypooler.persotools.files import *

'''
This script creates the founder population 8 (heterozygotes) + 1 (minor allele homozygote) samples 
to be used as reference panel.
Expected imputation result: nearly perfect. Using half of the number of founders might speed up the compuation time 
and memory complexity for imputation.
'''

CHROM_TABLE = {'1A': '1',
               '1B': '2',
               '1D': '3',
               '2A': '4',
               '2B': '5',
               '2D': '6',
               '3A': '7',
               '3B': '8',
               '3D': '9',
               '4A': '10',
               '4B': '11',
               '4D': '12',
               '5A': '13',
               '5B': '14',
               '5D': '15',
               '6A': '16',
               '6B': '17',
               '6D': '18',
               '7A': '19',
               '7B': '20',
               '7D': '21'}

filein = '/home/camille/MagicWheat/src/genotypooler/examples/Chr1.Founders.vcf.gz'
print('\n'.ljust(80, '*'))
print('The following arguments were parsed from file:\n')
print(filein)
print('\n'.ljust(80, '*'))

fileout = '/home/camille/MagicWheat/src/genotypooler/examples/8+1-beagle/Chr1.Founders.vcf'

truegt = os.path.expanduser(filein)

dftruegt = vcfdf.PandasMixedVCF(truegt, format='GT', indextype='id')

dfgenin = dftruegt.genotypes()

founders = dfgenin.columns
dfgenleft = dfgenin[founders[::2]].applymap(lambda x: sum(x) // 2)
dfgenright = dfgenin[founders[1::2]].applymap(lambda x: sum(x) // 2)

dictgenout = {}
for col_left, col_right in zip(dfgenleft.columns, dfgenright.columns):
    dictgenout[col_left.split('_')[0]
               + '_'
               + col_right.split('_')[0]] = list([(gleft, gright)
                                                  for gleft, gright in zip(dfgenleft[col_left],
                                                                           dfgenright[col_right])])
dfminorin = dftruegt.minor_allele
print('What ID for variants?', dfminorin.index)
dictgenout['Dummy_Dummy'] = list([(gmin, gmin) for gmin in dfminorin.values.flatten()])

dfgenout = pd.DataFrame.from_dict(dictgenout)
dfgenout.index = dfgenin.index

print(dfgenout)

VariantFileOut = vcfdf.PandasToVCF(truegt, fileout, dfgenout, 'GT', True)
VariantFileOut.write()

