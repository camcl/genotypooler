import subprocess
import numpy as np
import os

from genotypooler.poolSNPs import parameters as prm
from genotypooler.persotools.files import *

"""
Bash commands for bcftools manipulations written as Python-functions.
Increase readability.
(pysam should provide samtools commands but it does not work)
"""


def bgzip(f_vcf: str, f_gz: str, wd: str) -> None:
    """
    Bgzip a vcf file into a vcf.gz and checks if creation succeeded.
    :param f_vcf: input vcf file name
    :param f_gz: output vcf.gz file name
    :param wd: path to working directory
    :return: None
    """
    cmd = ' '.join(['bcftools',
                    'view -Oz -o',
                    f_gz,
                    '--threads {}'.format(os.cpu_count()),
                    f_vcf
                    ])
    subprocess.run(cmd, shell=True, cwd=wd)
    print('{}:\r\n File created? -> {}'.format(os.path.join(wd, f_gz),
                                               check_file_creation(wd, f_gz)))


def sort(f_gz: str, wd: str) -> None:
    """
    Sorts a VCF-file per increasing marher position on chromosom.
    :param f_gz: input file name (vcf or vcf.gz)
    :param wd: path to working directory
    :return: None
    """
    cmd = ' '.join(['bcftools',
                    'sort -Oz -o',
                    f_gz,
                    f_gz
                    ])
    subprocess.run(cmd, shell=True, cwd=wd)
    print('{}:\r\n File sorted? -> {}'.format(os.path.join(wd, f_gz),
                                              check_file_creation(wd, f_gz)))


def index(f_gz: str, wd: str) -> None:
    """
    Index a VCF-file and replace old .csi file if it exixts.
    :param f_gz: input file name (vcf or vcf.gz)
    :param wd: path to working directory
    :return: None
    """
    cmd = ' '.join(['bcftools',
                    'index -f',
                    f_gz
                    ])
    subprocess.run(cmd, shell=True, cwd=wd)
    print('{}:\r\n File indexed? -> {}'.format(os.path.join(wd, f_gz),
                                               check_file_creation(wd, f_gz + '.csi')))


def sampling(f_gz: str, f_out: str, f_samp: str, wd: str) -> None:
    """
    Creates a VCF-file (vcf.gz) for a chosen set of samples.
    Samples should be stored in a text-file, one sample name per line.
    Checks if sampled file was created.
    :param f_gz: input VCF-file name
    :param f_out: output VCF-file name (vcf.gz)
    :param f_samp: samples file name (text file)
    :param wd: path to working directory
    :return: None
    """
    cmd = ' '.join(['bcftools',
                    'view -Oz -o',
                    f_out,
                    '-S {}'.format(f_samp),
                    '--threads {}'.format(os.cpu_count()),
                    f_gz
                    ])
    subprocess.run(cmd, shell=True, cwd=wd)
    print('{}:\r\n File created? -> {}'.format(os.path.join(wd, f_out),
                                               check_file_creation(wd, f_out)))


def pick_sample(indir, filein, sampleid):
    """
    Creates a sample-wise VCF file from a main one.
    Output file is written in a 'tmp' directory.
    :param indir: path to directory with main file
    :param filein: name of the main file
    :param sampleid: ID of sample to pick

    Usage example with multiprocessing:
    ```
    import multiprocessing as mp

    study_pop = pd.read_csv('/home/camille/PoolImpHuman/data/20200827/study.population', header=None, names=['sampleID'])

    sampleID = study_pop.iloc[0]['sampleID']
    listID = study_pop['sampleID'].values.tolist()

    # Args to pass to multiprocessing
    args_pick_sample = [['/home/camille/PoolImpHuman/data/20210622',
                         'IMP.chr20.snps.gt.vcf.gz',
                         s] for s in listID]
    workers_pick_sample = mp.pool.Pool(processes=3)
    _ = workers_pick_sample.starmap(pick_sample, args_pick_sample)
        ```
    """
    # create tmp folder
    if not os.path.exists(os.path.join(indir, 'tmp')):
        os.mkdir(os.path.join(indir, 'tmp'))
    fileout = os.path.join(indir, 'tmp', 's{}.{}'.format(sampleid, filein))
    cmd1 = ' '.join(['bcftools view',
                     '-s',
                     sampleid,
                     '-Oz -o',
                     fileout,
                     filein
                     ])
    cmd2 = ' '.join(['bcftools index',
                     '-f',
                     fileout
                     ])
    subprocess.run(cmd1, shell=True, cwd=indir)
    subprocess.run(cmd2, shell=True, cwd=indir)


def stratified_aaf_sampling(f_gz: str, wd: str, binning: bool = False) -> None:
    """
    Samples markers from a VCF-file by selecting an equal number of
    markers for each 0.1-MAF interval. Total number of markers depends
    on the global chunk size parameter.
    :param f_gz: input VCF-file with all markers
    :param wd: path to working directory
    :param binning: boolean for splitting the original data set into AAF-binned layers
    :return: None
    """
    os.chdir(wd)
    print(wd)
    extract_header(os.path.join(prm.DATA_PATH, 'gt', 'ALL.chr20.snps.gt.vcf.gz'),
                   'headers.ALL.chr20.snps.gt.chunk{}.strat.vcf'.format(prm.CHK_SZ),
                   wd)
    bins = np.arange(0.0, 1.0, 0.1)
    for b in bins:
        print('interval for AAF: [{}, {})'.format(b, b+0.99))  # 0.99 and not 0.1: avoid intervals q/Q overlapping
        cmd1 = ' '.join(['bcftools',
                         'view -Oz -o',
                         f_gz.replace('.vcf.gz', '.maf_{}_{}.vcf.gz'.format(b, b+0.99)),
                         '-q {} -Q {}'.format(b, b+0.99),
                         f_gz
                         ])
        #TODO: what is done here?
        tmp = ' '.join(['cat ./bins_for_chunk10000/chunk{}.vcf'.format(b),
                        '| sort -R | head -{}'.format(prm.CHK_SZ // len(bins)),
                        '> chunk{}.vcf'.format(b)
                        ])  # for generating strat_chunk1000 from strat_chunk10000

        cmd2 = ' '.join(['bcftools view -H',
                         f_gz.replace('.vcf.gz', '.maf_{}_{}.vcf.gz'.format(b, b+0.99)),
                         '| sort -R | head -{}'.format(prm.CHK_SZ // len(bins)),
                         '> chunk{}.vcf'.format(b)
                         ])

        if binning:
            subprocess.run(cmd1, shell=True, cwd=wd)
        subprocess.run(tmp, shell=True, cwd=wd)
        # subprocess.run(cmd2, shell=True, cwd=wd)

    subprocess.run(' '.join(['cat headers.ALL.chr20.snps.gt.chunk{}.strat.vcf '.format(prm.CHK_SZ),
                             ' '.join(['chunk{}.vcf'.format(i) for i in bins]),
                             '> TMP.chr20.snps.gt.strat.vcf']),
                   shell=True,
                   cwd=wd)

    delete_file('ALL.chr20.snps.gt.chunk{}.vcf.gz.csi'.format(prm.CHK_SZ))

    bgzip('TMP.chr20.snps.gt.strat.vcf'.format(prm.CHK_SZ),
          'ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(prm.CHK_SZ),
          wd)
    sort('ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(prm.CHK_SZ),
         wd)
    index('ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(prm.CHK_SZ),
          wd)
    delete_file('headers.ALL.chr20.snps.gt.chunk{}.strat.vcf'.format(prm.CHK_SZ))
    print('Stratification ended\n\n'.ljust(80, '.'))
    # delete_file('TMP.chr20.snps.gt.strat.vcf')
    # for i in bins:
    #     delete_file('chunk{}.vcf'.format(i))


def rename_samples(file_in: str, file_out:str, wd: str, suffix:str) -> None:
    """
    Rename samples of a file by adding suffix to the old names.
    :param file_in:
    :param file_out:
    :param suffix:
    :return:
    """
    cmd0 = ' '.join(['bcftools query -l',
                     '-l',
                     file_in,
                     '> tmp.samples.get_names.txt'])
    subprocess.run(cmd0, shell=True, cwd=wd)

    with open(os.path.join(wd, 'tmp.samples.get_names.txt'), 'r') as get_names:
        names = get_names.readlines()

    with open(os.path.join(wd, 'tmp.samples.set_names.txt'), 'w') as set_names:
        lines = ['{} {}'.format(n.strip('\n\r'), n.strip('\n\r') + suffix + '\n') for n in names]
        set_names.writelines(lines)

    cmd1 = ' '.join(['bcftools reheader',
                    '-s tmp.samples.set_names.txt',
                     '-o',
                     file_out,
                     file_in])
    subprocess.run(cmd1, shell=True, cwd=wd)

    # delete_file(os.path.join(wd, 'tmp.samples.get_names.txt'))
    # delete_file(os.path.join(wd, 'tmp.samples.set_names.txt'))


def extract_header(f_gz: str, f_head: str, wd: str) -> None:
    """
    Extract header from VCF to uncompressed VCF format (text file)
    :param f_gz: input vcf.gz file name
    :param f_head: output vcf file name
    :param wd: path to working directory
    :return: None
    """
    cmd = ' '.join(['bcftools',
                    'view -h -Ov -o',
                    f_head,
                    f_gz
                    ])
    subprocess.run(cmd, shell=True, cwd=wd)
    print('{}:\r\n File created? -> {}'.format(os.path.join(wd, f_head),
                                               check_file_creation(wd, f_head)))


def chunk_markers(f_gz: str, chk_sz: int, wd: str) -> None:
    """

    :param f_gz:
    :param chk_sz:
    :param wd:
    :return:
    """
    cmd = ' '.join(['bcftools',
                    'view -H',
                    f_gz,
                    '| sort -R | head -{}'.format(str(chk_sz)),
                    '> chunk_{}.vcf'.format(str(chk_sz))
                    ])
    subprocess.run(cmd, shell=True, cwd=wd)
    print('{}:\r\n File created? -> {}'.format(os.path.join(wd,
                                                            'chunk_{}.vcf'.format(str(chk_sz)),
                                               check_file_creation(wd,
                                                                   'chunk_{}.vcf'.format(str(chk_sz))))))


def cat(flist_in: list, f_out: FilePath, wd: str):
    """
    Merge text-readable files
    Merge header.vcf and body.vcf parts.
    :param flist_in:
    :param f_out:
    :return:
    """
    cmd = ' '.join(['cat',
                    ' '.join([f for f in flist_in]),
                    '>',
                    f_out
                    ])
    subprocess.run(cmd, shell=True, cwd=wd)
    print('{}:\r\n File created? -> {}'.format(os.path.join(wd, f_out),
                                               check_file_creation(wd, f_out)))


def concat(flist_in: list, f_out: FilePath, wd: str) -> None:
    """
    Concatenate or combine VCF/BCF files.
    All source files must have the same sample columns appearing in the same order.
    Can be used, to concatenate chromosome VCFs into one VCF, or combine a SNP VCF and an indel VCF into one.
    The input files must be sorted by chr and position. The files must be given in the correct order to produce
    sorted VCF on output unless the -a, --allow-overlaps option is specified. With the --naive option, the files
    are concatenated without being recompressed, which is very fast..
    :param flist_in: list of vcf files names
    :param f_out: output vcf file name
    :param wd: path to working directory
    :return: None
    """
    cmd = ' '.join(['bcftools',
                    'concat -aOv -o',
                    f_out,
                    ' '.join([f for f in flist_in])
                    ])
    subprocess.run(cmd, shell=True, cwd=wd)
    print('{}:\r\n File created? -> {}'.format(os.path.join(wd, f_out),
                                               check_file_creation(wd, f_out)))


def get_first_pos(f: FilePath, wd: str) -> None:
    """
    Return position of the first variant in the file
    :param f:
    :param wd:
    :return:
    """
    cmd = ' '.join(['bcftools query -f',
                    '"%CHROM:%POS\n"',
                    f,
                    '| ghead -1'
                    ])
    process = subprocess.run(cmd, shell=True, cwd=wd, capture_output=True)
    return process.stdout


def view_one_variant(f: FilePath, pos: str, chrom: str = '20', wd=os.getcwd()) -> None:
    """
        Extract a given variant on CHROM:POS identifier and create a VCF file from it.
        :param f:
        :param wd:
        :return:
        """
    filout = 'snp.{}:{}.{}'.format(chrom, pos, os.path.basename(f))
    cmd = ' '.join(['bcftools view ',
                    '-r {}:{}'.format(chrom, pos),
                    '-Oz -o {}'.format(os.path.join(wd, filout)),
                    f
                    ])
    process = subprocess.run(cmd, shell=True, cwd=wd)
    return os.path.join(wd, filout)