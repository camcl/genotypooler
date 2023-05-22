import sys, os
import argparse
import timeit

# force PYTHONPATH to look into the project directory for modules
rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)
sys.path.insert(0, rootdir + '/genotypooler')

from genotypooler.poolSNPs import poolvcf
from genotypooler.poolSNPs import pybcf

'''
Applies pooling simulation to a VCF file

* the input VCF file contains variants of type SNP only, with genotype formatted as GT,
* the output VCF file contains the same variants, formatted as GL,
* the number of samples is a multiple of 16 (size of the pooling blok in our design).
The decoding step of pooling is adaptive with respect to the pooling pattern observed (see README.md)
* the samples are assumed to be sorted in row-order flattened blocks order e.g. the 16 first columns in the VCF file
correspond  to the samples assigned to the first block. 
Samples 1-4 form the first pool in the block, samples 5-8 the second pool, and so on.
* the output file has to be written to unbgzipped format (.vcf) and then compressed to 
bgzipped format (.vcf.gz) with bcftools.

For VCF-file bigger than some dozen of thousands of variants, pooling can be parallelized.

Command line usage (assuming the current directory is VCFPooling/examples)
$ python3 -u pooling-ex.py <path-to-file-in> <path-to-file-out> <decoding-format>

Examples for GP decoding format (assuming current directory is /examples):
python3 -u pooling-ex.py Chr1.SNPs.pruned.nomiss.vcf.gz Chr1.SNPs.pruned.nomiss.pooled-example.vcf.gz GP
'''

### COMMAND-LINE PARSING AND PARAMETERS
parser = argparse.ArgumentParser(description='Run pooling simulation'
                                             'on the whole set of samples')
parser.add_argument('pathin', metavar='in', type=str, help='File to pool', default=None)
parser.add_argument('pathout', metavar='out', type=str, help='Pooled file', default=None)
parser.add_argument('formatto', metavar='fmt', type=str, help='Genootype format to decode to (GP or GT)', default='GP')  # default option does not work


argsin = parser.parse_args()
filin = argsin.pathin
filout = argsin.pathout
fmt = argsin.formatto.upper()
plookup = os.path.join(os.getcwd(), 'adaptive_gls.csv')  # look-up table for converting pooled GT to GL

print('\n'.ljust(80, '*'))
print('Input file = {}'.format(os.path.expanduser(argsin.pathin)))
print('Output file = {}'.format(os.path.expanduser(argsin.pathout)))
print('\n'.rjust(80, '*'))

# make sure to write to .vcf
if filout.endswith('.gz'):
    vcfout = filout[:-3]
if filout.endswith('.vcf'):
    vcfout = filout

### SIMULATE POOLING
start = timeit.default_timer()
if fmt == 'GP':
    poolvcf.pysam_pooler_gp(filin, vcfout, plookup, os.getcwd())
if fmt == 'GT':
    poolvcf.pysam_pooler_gt(filin, vcfout, plookup, os.getcwd())

print('\r\nTime elapsed --> ', timeit.default_timer() - start)
