import sys, os
import pysam

rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

nb_cores = os.cpu_count()

from genotypooler.persotools.files import delete_file, FilePath

'''
Classes for parallelized file processing: Read main VCF-file and write chunks from it
1) Implementation based on cyvcf2 file handling (e.g. pooling)
2) Implementation based on pysam file handling (SNPsPool to reimplement to rewrite)

Steps:
* Read main VCF-file: iterate over variants and repack them into chunks
* Write packed chunks to independent VCF-files
* Apply pooling simulation to a population at one marker
* Scale up pooling to a population over all variants in a VCF-file
'''


class PysamVariantCallGenerator(object):
    """
    Generates single-type formatted calls of variants
    """

    def __init__(self, vcfpath: FilePath, format: str = None):
        """
        :param vcfpath:
        :param indextype: identifier for variants: 'id', 'chrom:pos'.
        Must be 'chrom:pos' if the input has been generated by Phaser
        """
        self.path = vcfpath
        self.fmt = format

    def __iter__(self):
        vcfobj = pysam.VariantFile(self.path)
        for var in vcfobj:
            yield [g[self.fmt] for g in var.samples.values()]


class PysamVariantChunkGenerator(object):
    """
    Generates chunks of single-type formatted calls of variants
    """

    def __init__(self, vcfpath: FilePath, format: str = None, chunksize: int = None):
        """
        :param vcfpath:
        :param indextype: identifier for variants: 'id', 'chrom:pos'.
        Must be 'chrom:pos' if the input has been generated by Phaser
        """
        self.path = vcfpath
        self.fmt = format
        self.chksz = chunksize
        self.chrom = [*pysam.VariantFile(self.path).header.contigs][0]
        # extract chrom, works if only 1 chrom in the file
        self.pack = True
        self.newpos = 1  # for valid self.newpos - 1 at the start of the first chunk

    def chunk(self, chunksize: int, newpos: int):
        """Build generators of variants calls"""
        iterator = pysam.VariantFile(self.path)
        try:
            for i, v in enumerate(iterator.fetch(contig=self.chrom, start=newpos - 1, reopen=False)):
                # newpos - 1: avoids first variant truncation in the next chunk
                var = v
                if i == chunksize:
                    break
                yield var

        except StopIteration:
            print('Could not build chunk')

    def incrementer(self, chunksize: int):
        """update position and packing bool"""
        iterator = pysam.VariantFile(self.path)
        try:
            for i, v in enumerate(iterator.fetch(contig=self.chrom, start=self.newpos - 1, reopen=False)):
                # self.newpos - 1: avoids first variant truncation in the next chunk
                var = v
                if i == chunksize:
                    self.newpos = var.pos
                    break
            if var.pos != self.newpos:  # reached EOF
                self.pack = False

        except StopIteration:
            self.newpos = None
            self.pack = False

        finally:
            return self.newpos, self.pack

    def chunkpacker(self):
        while self.pack:
            chk = self.chunk(self.chksz, self.newpos)
            # function output and included sttributes updates NOT unpacked hence NOT updated
            self.newpos, self.pack = self.incrementer(self.chksz)
            yield chk
            print(self.newpos, self.pack)



