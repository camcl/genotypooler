"""
Git: Commit this script to a new branch?

Compare datasets 1 and 2: what genotypes (variants x samples) got changed priors from 1 to 2?
Epsilon for change = 1e-05, unchanged if delta <= epsilon.

Some kind of masked array?
Find a way to write the masked arrays to VCF? (just reuses the quantiles scripts on the "masked" VCF then)
How to represent a masked value in a VCF-friendly format? Missing entry?
How to compute the metrics with missing entries in the data set?
"""

import numpy as np
import pandas as pd
import pysam
import timeit
import os, sys
import argparse

rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from genotypooler.poolSNPs import dataframe as vcfdf
from genotypooler.poolSNPs import pybcf

changes = True  # False --> unchanged priors (not working?)
prev_cycle = 1
parent_dir = "/home/camille/IterDecodeImpute/runs/poolimputeSNPs/results/data/1/"
vcf1 = parent_dir + f"cycle{prev_cycle}/STU.Chr1.SNPs.pruned.sorted.pooled.vcf.gz"  # (log) GL
vcf2 = parent_dir + f"cycle{prev_cycle + 1}/STU.Chr1.SNPs.pruned.sorted.pooled.vcf.gz"  # (log) GL
vcfout = parent_dir + f"cycle{prev_cycle + 1}/STU.Chr1.SNPs.pruned.sorted.pooled.changed_priors.vcf"  # write as text i.e. uncompressed


class VariantRecordDiff(object):
    """
    Computes indicator for masking values if (un)changed priors between two cycles.
    """
    def __init__(self, var: pysam.VariantRecord,
                 format_to: str = 'GP', wd: str = os.getcwd()):
        """
        ...
        """
        self.var = var.copy()
        assert ('GL' in dict(self.var.format)),  'Reading from other format than GP not implemented'
        self.fmt_to = format_to  # write to GP
        self.wd = wd

    @property
    def samples(self):
        """Samples in the input VCF file"""
        return list(self.var.samples.keys())

    @property
    def genotypes(self):
        """Estimated Genotype Likelihoods GL only"""
        return np.asarray([v['GL'] for v in self.var.samples.values()])

    def _diff_prior_null(self, prior1: np.ndarray, prior2: np.ndarray) -> bool:
        """Test whether two priors represented as GL(RR, RA, AA) are equal"""
        # if np.equal(prior1, prior2).all():
        #     return True
        # else:
        #     return False
        return np.allclose(prior1, prior2, rtol=0.0, atol=1e-05)

    def new_var(self, var1, var2, changed) -> str:
        """
        Outputs a string representation of a pooled variant
        since pysam.VariantRecord objects are not writable
        """
        if self.fmt_to == 'GP':
            _genotypes = self.genotypes  # .reshape((self.genotypes.shape[0], 3))  # 3 --> GP(RR, RA, AA)
            info_fields = ['='.join([str(k), str(np.asarray(v).flatten()[0])]) for k, v in self.var.info.items()]
            info = ';'.join([kv for kv in info_fields])
            str_var = '''{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'''.format(self.var.chrom,
                                                                      self.var.pos,
                                                                      self.var.id,
                                                                      self.var.ref,
                                                                      self.var.alts[0],
                                                                      int(self.var.qual) if self.var.qual is not None else '.',
                                                                      self.var.filter.keys()[0] if self.var.filter.keys() != [] else '.',
                                                                      info,
                                                                      'GP')
            for _i, v in enumerate(self.var.samples.values()):
                pr1, pr2 = np.array(var1.samples.values()[_i]['GL']), np.array(var2.samples.values()[_i]['GL'])
                if not (changed and self._diff_prior_null(pr1, pr2)):
                    str_var = str_var + '\t' + ','.join(np.power(10, _genotypes[_i]).astype(str))
                else:
                    str_var = str_var + '\t' + ','.join('.')

            str_var = str_var + '\n'

        return str_var


class VariantFileDiff(object):
    """
    Writes a new VariantFile with GP (possibly missing data '.').
    """
    def __init__(self, path_in: str, path_out: str,
                 format_to: str = 'GP', wd: str = os.getcwd()):
        """
        ...
        """
        self.vcf_in = pysam.VariantFile(path_in)
        self.path_in = path_in
        self.path_out = path_out
        self.fmt_to = format_to.upper()
        assert (self.fmt_to == 'GP'), 'Other formats than GP not handled'
        self.wd = wd
        self.header = None
        self.data = None
        self.n_variants = 0

    def _new_header(self):
        """Modifies VCF header in-place and replace GL format by GP"""
        if self.fmt_to == 'GP':
            for hrec in list(self.vcf_in.header.records):
                drec = dict(hrec)
                if drec != {} and drec['ID'] == 'GL':
                    hrec.remove()
            self.vcf_in.header.add_line(
                '##FORMAT=<ID=GP,Number=G,Type=Float,Description="Estimated Genotype Probability">')
        self.header = iter([hrec for hrec in self.vcf_in.header.records])

    def _new_data(self, data1, data2, changed):
        """Masks data"""
        print('Processing data in {}'.format(self.path_in).ljust(80, '.'))
        pvar = []
        tm = timeit.default_timer()
        for n, (var1, var2) in enumerate(zip(data1.fetch(), data2.fetch())):
            varout = VariantRecordDiff(var2, self.fmt_to)
            pvar.append(varout.new_var(var1, var2, changed))
            self.n_variants += 1
            if n % 1000 == 0:
                print('{} variants processed in {:06.2f} sec'.format(n + 1, timeit.default_timer() - tm).ljust(80, '.'))
        self.data = iter(pvar)

    def write(self, data1, data2, changed) -> None:
        """Writes masked data into an output file"""
        # load header and data to write
        self._new_header()
        self._new_data(data1, data2, changed)
        # write as text-like file
        print('\r\nWriting metadata in {}'.format(self.path_out).ljust(80, '.'))
        vcf_out = pysam.VariantFile(self.path_out, 'w', header=self.vcf_in.header)
        vcf_out.close()
        print('\r\nWriting data in {}'.format(self.path_out).ljust(80, '.'))
        with open(self.path_out, 'a') as f_out:
            f_out.writelines(self.data)
        print('Writing data in {}: Done'.format(self.path_out).rjust(80, '.'))
        # compress and index the VCF file
        pybcf.bgzip(self.path_out, self.path_out + '.gz', wd=self.wd)
        pybcf.index(self.path_out + '.gz', wd=self.wd)
        os.remove(self.path_out)


vf1 = pysam.VariantFile(vcf1)
vf2 = pysam.VariantFile(vcf2)
vfout = VariantFileDiff(vcf2, vcfout)
vfout.write(vf1, vf2, changes)

# Test reading the file that was written
newvf = pysam.VariantFile(vcfout + '.gz')
for n, record in enumerate(newvf.fetch()):
    if n < 10:
        print(n + 1, end='\t')
    for samp, genos in record.samples.items():
        if n < 10:
            print(samp, genos['GP'], end='\t')
    if n < 10:
        print('\n')

# Create a mask from the data (test)
dfobj = vcfdf.PandasMixedVCF(vcfout + '.gz', format='GP', indextype='id', mask=None)

gdata = dfobj.genotypes()
gdf = gdata.applymap(lambda x: None if x[0] is None else x)
garr = gdf.values
gbool = gdf.applymap(lambda x: True if x is None else False)  # np.full_like(garr, True)
gmasked = np.ma.array(gdata, mask=gbool)

print('\n')
print(garr)
print('\n')
print(gbool)
print('\n')
print(gmasked)
print('\n')

a = np.random.random(garr.shape)
ama = np.ma.masked_array(a, mask=gbool)
print(ama)
print()
print(a.mean(axis=0)[:5])
print()
print(pd.DataFrame(ama).mean(axis=0, skipna=True).head(5).values)
print(ama.mean(axis=0)[:5])  # masked values are ignored for computing the mean
