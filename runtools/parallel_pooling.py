import sys, os
import numpy as np
from itertools import starmap, repeat
import shutil
import multiprocessing as mp
import argparse
import timeit

home_dir = os.path.expanduser("~")
proj_dir = os.path.join(home_dir, '../')

# force PYTHONPATH to look into the project directory for modules
rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from genotypooler.poolSNPs import poolvcf
from genotypooler.poolSNPs import pybcf

from genotypooler.persotools.files import delete_file, mkdir

'''
Parallelized file processing simulating pooling with adaptive genotype likelihoods.
For large VCF-files with some dozen of thousands or a million of variants.
If million of variants, it is recommened to run this script
on a cluster of servers e.g. UPPMAX resources (https://www.uppmax.uu.se/) .

Steps:
* Read main VCF file and write chunks: bash script based on bcftools
* Simulate pooling on chunks and decode into adaptive GP (chunks have to be in a temporary directory inside the working directory
* Merge pooled chunks back to read-for-use compressed VCF file

Usage (from genotypooler/runtools folder): 
$ python3 parallel_pooling.py /home/camille/MagicWheat/runs/devtests/data/1 STU.Chr1.SNPs.pruned.sorted.vcf.gz tmp/ STU.Chr1.SNPs.pruned.sorted.pooled.vcf.gz 4
'''

### COMMAND-LINE PARSING AND PARAMETERS
parser = argparse.ArgumentParser(description='Run parallelized pooling simulation'
                                             'on the whole set of samples')
parser.add_argument('pathwd', metavar='wd', type=str, help='Path to working directory (files to process)', default=None)
parser.add_argument('infile', metavar='in', type=str, help='File to pool', default=None)
parser.add_argument('tempdir', metavar='tmp', type=str, help='File to pool', default=None)
parser.add_argument('outfile', metavar='out', type=str, help='Pooled file', default=None)
parser.add_argument('cores', metavar='nc', type=int, nargs='?', help='Number of cores to use', default=None)

argsin = parser.parse_args()

nb_cores = os.cpu_count() if argsin.cores is None else argsin.cores
tmp_path = os.path.join(argsin.pathwd, os.path.basename(argsin.tempdir))
if not os.path.exists(tmp_path):
    os.makedirs(tmp_path)

print('\n'.ljust(80, '*'))
print('Number of cpu to be used = {}'.format(nb_cores))
print('Input file = {}'.format(os.path.join(argsin.pathwd, argsin.infile)))
print('Output file = {}'.format(os.path.join(argsin.pathwd, argsin.outfile)))
print('\n'.rjust(80, '*'))

### SIMULATE POOLING ON PACKED DATA CHUNKS
# data_dir = snic_proj
# os.chdir(data_dir)
fingz = argsin.infile  # os.path.expanduser(argsin.pathin)  # /home/usr
basin = os.path.basename(fingz).rstrip('.gz').replace('.gt', '.gl')
basout = os.path.basename(argsin.outfile)


start = timeit.default_timer()
prename = 'pack'

# Clean tmp directory: keep only GLs values and chunks to process (index the files)
files0 = []
for f0 in os.listdir(tmp_path):
    if f0.startswith(prename) and f0.endswith('.vcf.gz'):
        print('Indexing ', f0)
        pybcf.index(f0, tmp_path)
        files0.append(f0)
    elif f0 != 'adaptive_gls.csv':
        print('Removing ', f0)
        os.remove(os.path.join(tmp_path, f0))

indices = np.arange(len(files0))

print('\r\n{} files found will be pooled'.format(len(files0)).ljust(80, '.'))
print([os.path.join(tmp_path, f0) for f0 in files0])

# Create names for the pooled files
files1 = [os.path.join(tmp_path,
                       'pooled.{}'.format(f0).rstrip('.gz').replace('.gt', '.gl'))
          for f0 in files0]
print(files1)

args1 = list(zip([f0 for f0 in files0],
                 [f1 for f1 in files1],
                 repeat(os.path.join(tmp_path, 'adaptive_gls.csv'), len(indices)),  # path to lookup table
                 repeat(tmp_path, len(indices))))  # directory for temporary output

print(args1)

with mp.Pool(processes=nb_cores) as mpool:
    _ = all(mpool.starmap(poolvcf.pysam_pooler_gp, args1))

print('\r\nTime elapsed --> ', timeit.default_timer() - start)

### CONCATENATE RESULTS
files2 = ['{}.gz'.format(f1) for f1 in files1]
pybcf.concat(files2, basin, tmp_path)
pybcf.sort(basin, tmp_path)
pybcf.bgzip(basin, os.path.join(argsin.pathwd, basout), tmp_path)
pybcf.index(basout, argsin.pathwd)
delete_file(os.path.join(tmp_path, basin))

print('\r\nTime elapsed --> ', timeit.default_timer() - start)
