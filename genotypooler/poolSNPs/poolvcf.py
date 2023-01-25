import sys, os

# force PYTHONPATH to look into the project directory for modules
rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from genotypooler.poolSNPs.pooler import *
from genotypooler.poolSNPs import pybcf
from genotypooler.poolSNPs import utils

import numpy as np
import timeit
import pysam

"""
Classes for applying pooling simulation to a VCF file using Pysam
"""

np.random.seed(123)  # fix the shuffling order


class ShuffleSplitVCF(object):
    """
    Based on sklearn ShuffleSplit. Shuffle samples in a VCF file and split them into
    a reference panel population and a study (target) population.
    Adapted for small VCF files only (e.g 10-20,000 variants).
    """
    def __init__(self, design_matrix: np.ndarray, vcf_in: str, stu_size: float = 0.1, wd: str = os.getcwd()):
        """
        Per default, the number of splits is chosen to be 2 (reference panel + study population)
        such that the study population size is:
        * a multiple of the block's size chosen in the design
        * approximately 1/10 of the reference panel
        Samples that are not assigned to the study population are included in the reference panel.
        Shuffled samples are assigned to populations and blocks, and written back to VCF files
        in block-order (read/write in row-major order for a block).
        :param design_matrix: pooling design chosen expressed as a pooler.SingleBlockDecoder.matrix object
        :param vcf_in: path to the VCF file to process
        :param stu_size: relative size of the study population vs. reference panel
        :param wd: directory for writing the output files
        """
        self._n_splits = 2  # 1 split to reference panel + 1 split to study population
        self.design = design_matrix
        self.filein = vcf_in
        self.stu_size = stu_size
        self.ref_pan = os.path.join(wd, 'reference.panel')
        self.stu_pop = os.path.join(wd, 'study.population')
        self.ref_prefix = 'REF'
        self.stu_prefix = 'IMP'
        self.wd = wd

    @property
    def samples(self):
        """Samples in the input VCF file"""
        return pysam.VariantFile(self.filein).header.samples

    def _split_samples(self):
        """Assign shuffled samples to reference/study by writing their ID in separate files"""
        # Calculate the study population size
        n_idv = len(self.samples)
        block_size = self.design.shape[1]
        n_blocks, left_over = divmod(n_idv, block_size)  # design constraint
        n_samples = n_blocks * block_size
        stu_blocks = utils.ppcm(n_samples, block_size) * self.stu_size // block_size  # number of blocks for the target
        n_stu = int(stu_blocks * block_size)  # number of samples for the target set
        # Shuffle samples order
        samples = list(self.samples)
        np.random.seed(123)  # fix the shuffling order
        rng = np.random.default_rng()
        rng.shuffle(samples)
        # Assign samples to populations
        with open(self.ref_pan, 'w') as fpan:
            fpan.writelines([s + '\n' for s in samples[n_stu:]])
        with open(self.stu_pop, 'w') as fstu:
            fstu.writelines([s + '\n' for s in samples[:n_stu]])

    def split_file(self, base_name_out: str):
        """Write reference an target populations to files"""
        self._split_samples()
        pybcf.sampling(self.filein,
                       self.ref_prefix + '.' + base_name_out,
                       self.ref_pan,
                       wd=self.wd)
        pybcf.index(self.ref_prefix + '.' + base_name_out, wd=self.wd)
        pybcf.sampling(self.filein,
                       self.stu_prefix + '.' + base_name_out,
                       self.stu_pop,
                       wd=self.wd)
        pybcf.index(self.stu_prefix + '.' + base_name_out, wd=self.wd)


class VariantRecordPooler(object):
    """
    Applies pooling simulation to samples' genotypes at a variant.
    """
    def __init__(self, design_matrix: np.ndarray, var: pysam.VariantRecord,
                 dict_lookup: dict, format_to: str, wd: str = os.getcwd()):
        """
        The NonOverlapping Repeated Block pooling design applied is provided with the design matrix.
        """
        self.dm = design_matrix
        self.var = var.copy()
        assert ('GT' in dict(self.var.format)),  'Pooling from other format than GT not implemented'
        self.lookup = dict_lookup
        self.fmt_to = format_to.upper()  # decode into GT or GP
        self.wd = wd

    @property
    def samples(self):
        """Samples in the input VCF file"""
        return list(self.var.samples.keys())

    @property
    def genotypes(self):
        """True genotypes GT only"""
        return np.asarray([v['GT'] for v in self.var.samples.values()])

    @property
    def n_blocks(self) -> int:
        """Number of pooling blocks from the sampples"""
        assert len(self.samples) % self.dm.shape[1] == 0  # check all samples fit in some block
        return len(self.samples) // self.dm.shape[1]

    def _encode(self) -> np.ndarray:
        dse = Design(blocks=self.n_blocks)
        dme = dse.matrix
        enc = Encoder(dme)
        return enc.encode(self.genotypes.sum(axis=-1))

    def _decode(self) -> np.ndarray:
        x_shift = self.dm.shape[1]  # 16
        y_shift = self.dm.shape[0]  # 8
        varp = self._encode()
        return dict_blocks_decoder(self.n_blocks, varp.sum(axis=-1), y_shift, self.lookup, self.fmt_to.lower())

    def new_var(self) -> str:
        """
        Outputs a string representation of a pooled variant
        since pysam.VariantRecord objects are not writable
        """
        #_genotypes = self._decode().reshape((self.genotypes.shape[0], 3))  # 3 --> GP(RR, RA, AA)
        if self.fmt_to == 'GP':
            _genotypes = self._decode().reshape((self.genotypes.shape[0], 3))  # 3 --> GP(RR, RA, AA)
            # GP must be written as GL (literaly) for compatibility with Beagle
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
                                                                      'GL')  # GL for compatibility with Beagle
        # var.__str__()
        # '20\t61098\trs6078030\tC\tT\t100\tPASS\tAC=4;AF=0.287141;AN=32;NS=2504;DP=16880;EAS_AF=0.4415;AMR_AF=0.2709;AFR_AF=0.1823;EUR_AF=0.2167;SAS_AF=0.3538;AA=.|||;VT=SNP\tGT\t0|0\t0|0\t0|0\t0|0\t0|1\t0|0\t1|0\t0|0\t0|0\t0|0\t0|0\t0|1\t0|0\t0|1\t0|0\t0|0\n'
        # self.var.format = 'GP'
        # AttributeError: attribute 'format' of 'pysam.libcbcf.VariantRecord' objects is not writable
        elif self.fmt_to == 'GT':
            _genotypes = self._decode().reshape((self.genotypes.shape[0], 2))  # 3 --> Gt(allele 0, allele 1)
        for _i, v in enumerate(self.var.samples.values()):
            if self.fmt_to == 'GP':
                str_var = str_var + '\t' + ','.join(_genotypes[_i].astype(str))
                # self.var.format.clear() # Process finished with exit code 139 (interrupted by signal 11: SIGSEGV)
            elif self.fmt_to == 'GT':
                try:
                    v[self.fmt_to] = tuple(_genotypes[_i])  # in-place modification of the genotype values
                except ValueError:  # pysam does not accept -1 but set None for missing alleles
                    if _genotypes[_i].sum() == 0:  # case [1, -1]
                        v[self.fmt_to] = (1, None)
                    elif _genotypes[_i].sum() == -1:  # case [0, -1]
                        v[self.fmt_to] = (0, None)
                    else:  # case [-1, -1]
                        v[self.fmt_to] = (None, None)
                finally:
                    v.phased = False  # set phase to unphased after pooling
                str_var = self.var.__str__()
        if self.fmt_to == 'GP':
            str_var = str_var + '\n'
        return str_var


class VariantFilePooler(object):
    """
    Writes a new VariantFile.
    GP are converted to GL format at writing for compatibility with the imputation methods.
    Add GL format to the header if necessary.
    """
    def __init__(self, design_matrix: np.ndarray, vcf_in: str, vcf_out: str,
                 dict_lookup: dict, format_to: str, wd: str = os.getcwd()):
        """
        The NonOverlapping Repeated Block pooling design applied is provided with the design matrix.
        Pooling from only GT genotype format to only GT or GP format implemented.
        """
        self.design = design_matrix
        self.vcf_in = pysam.VariantFile(vcf_in)
        self.path_in = vcf_in
        self.path_out = vcf_out
        self.lookup = dict_lookup
        self.fmt_to = format_to.upper()
        assert (self.fmt_to == 'GT' or self.fmt_to == 'GP'), 'Pooling to other formats than GT or GP not implemented'
        self.wd = wd
        self.header = None
        self.data = None
        self.n_variants = 0

    def _new_header(self):
        """Modifies VCF header in-place if necessary (new GL format from pooling)"""
        # GP header record
        # ##FORMAT=<ID=GP,Number=G,Type=Float,Description="Estimated Genotype Probability">
        if self.fmt_to == 'GP':
            for hrec in list(self.vcf_in.header.records):
                drec = dict(hrec)
                if drec != {} and drec['ID'] == 'GT':
                    hrec.remove()
            self.vcf_in.header.add_line(
                '##FORMAT=<ID=GL,Number=G,Type=Float,Description="Estimated Genotype Likelihood">')
            # GP must be written as GL (literaly) for compatibility with Beagle
        self.header = iter([hrec for hrec in self.vcf_in.header.records])

    def _new_data(self):
        """Simulates NORB pooling on the variants in the input file"""
        print('Pooling data in {}'.format(self.path_in).ljust(80, '.'))
        pvar = []
        tm = timeit.default_timer()
        for n, var in enumerate(self.vcf_in.fetch()):
            prec = VariantRecordPooler(self.design, var, self.lookup, self.fmt_to)
            pvar.append(prec.new_var())
            self.n_variants += 1
            if n % 1000 == 0:
                print('{} variants processed in {:06.2f} sec'.format(n + 1, timeit.default_timer() - tm).ljust(80, '.'))
        self.data = iter(pvar)

    def write(self) -> None:
        """Writes pooling simulation result into an output file"""
        # load header and data to write
        self._new_header()
        self._new_data()
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


class VariantRecordConverter(pysam.VariantRecord):
    """Converts format and genotypes"""
    def __new__(cls, var: pysam.VariantRecord, format_to: str, genotypes=None):
        cls.var = var.copy
        cls.fmt_to = format_to
        cls.genotypes = genotypes
        return super(VariantRecordConverter, cls).__new__(cls)

    def do(self):
        new_var = dict(zip(self.var.__dir__(), [None] * len(self.var.__dir__())))
        new_var['alleles'] = self.var.alleles
        new_var['alts'] = self.var.alts
        new_var['chrom'] = self.var.chrom
        new_var['contig'] = self.var.contig
        new_var['filter'] = self.var.filter
        new_var['format'] = 'GP' if self.fmt_to == 'GP' else 'GT'
        new_var['id'] = self.var.id
        new_var['info'] = self.var.info
        new_var['pos'] = self.var.pos
        new_var['qual'] = self.var.qual
        new_var['ref'] = self.var.ref
        new_var['rid'] = self.var.rid
        new_var['rlen'] = self.var.rlen
        new_var['samples'] = self.var.samples
        new_var['start'] = self.var.start
        new_var['stop'] = self.var.stop
        new_var['__class__'] = pysam.libcbcf.VariantRecord


def pysam_pooler_gp(file_in: str, file_out: str, path_to_lookup: str, wd: str):
    """
    Process a VCF file with NORB pooling simulation. Decode into GP based on a look up table.
    :param file_in: name of the file to be processed (.vcf.gz or .vcf only)
    :param file_out: name of the file to output (NO .gz)
    :param path_to_lookup: lookup table to use for GP decoding
    :param wd: path to the data directory
    """
    design = Design()
    dm = design.matrix

    #path_to_lookup = '/home/camille/PoolImpHuman/data/main'
    dict_gl = load_lookup_dict(path_to_lookup)

    poolf = VariantFilePooler(dm,
                              os.path.join(wd, file_in),
                              os.path.join(wd, file_out),
                              dict_gl,
                              'GP')

    tstart = timeit.default_timer()
    poolf.write()
    tstop = timeit.default_timer()
    print('Time for pooling {} variants = {} sec'.format(poolf.n_variants, tstop - tstart))


def pysam_pooler_gt(file_in: str, file_out: str, path_to_lookup: str, wd: str):
    """
    Process a VCF file with NORB pooling simulation. Decode into GP based on a look up table.
    :param file_in: name of the file to be processed (.vcf.gz or .vcf only)
    :param file_out: name of the file to output (NO .gz)
    :param path_to_lookup: lookup table to use for GP decoding
    :param wd: path to the data directory
    """
    design = Design()
    dm = design.matrix

    #path_to_lookup = '/home/camille/PoolImpHuman/data/main'
    dict_gl = load_lookup_dict(path_to_lookup)

    poolf = VariantFilePooler(dm,
                              os.path.join(wd, file_in),
                              os.path.join(wd, file_out),
                              dict_gl,
                              'GT')

    tstart = timeit.default_timer()
    poolf.write()
    tstop = timeit.default_timer()
    print('Time for pooling {} variants = {} sec'.format(poolf.n_variants, tstop - tstart))

