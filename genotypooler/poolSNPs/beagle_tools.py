import subprocess
import shutil
import os

from typing import *

from genotypooler.poolSNPs import parameters as prm
from genotypooler.poolSNPs import pybcf
from genotypooler.persotools.files import *

'''
Functions for running the steps necessary to Beagle execution
Return success status for each step
'''


def bgzip_working_files(dic: Dict[str, str], path_gt_files: str, path_gl_files: str, cd: str) -> None:
    """
    Bgzip, sort and index files created from pooloing-decoding (text files vcf formatted).
    :param dic:
    :param path_gt_files:
    :param path_gl_files:
    :param cd:
    :return:
    """
    print('\n\nBGZIP in {}'.format(path_gt_files).ljust(80, '.'))
    # process prm.GTGL == 'GT' format anyway: enables to start imputing with GT or GL indifferently
    print('{} compressed to {}'.format(os.path.join(path_gt_files, dic['vcf']),
                                       dic['gz']))
    pybcf.bgzip(os.path.join(path_gt_files, dic['vcf'].replace('.gl', '.gt')),
                os.path.join(path_gt_files, dic['gz'].replace('.gl', '.gt')),
                path_gt_files)  # bgzip the file in the corresponding GT folder for missing values
    # Delete the index file for avaoiding Warning
    # "[W::hts_idx_load2] The index file is older than the data file" when sorting
    delete_file(os.path.join(path_gt_files,
                             dic['gz'].replace('.gl', '.gt') + '.csi'))
    pybcf.sort(dic['gz'].replace('.gl', '.gt'), path_gt_files)
    pybcf.index(dic['gz'].replace('.gl', '.gt'), path_gt_files)
    #delete_file(path_gt_files + dic['vcf'])

    if prm.GTGL == 'GL':
        print('\n\nBGZIP in {}'.format(os.getcwd()).ljust(80, '.'))
        print('{} compressed to {}'.format(os.path.join(path_gl_files, dic['vcf']),
                                           dic['gz']))
        pybcf.bgzip(os.path.join(path_gl_files, dic['vcf']),
                    dic['gz'],
                    cd)  # bgzip the file in the corresponding GL folder for unknown genotypes
        delete_file(dic['gz'] + '.csi')
        pybcf.sort(dic['gz'], cd)
        pybcf.index(dic['gz'], cd)
        #delete_file(path_gl_files + dic['vcf'])


def create_ref_imp_lists(cd: str, sizeref: int = prm.NB_REF, sizeimp: int = prm.NB_IMP) -> None:
    """
    Creates text files with samples ID for the reference panel and the study population.
    1 sample ID per line.
    Text files saved in the data GT folder.
    :param cd: does not affect the process directly
    :param sizeref:
    :param sizeimp:
    :return:
    """
    print('Set size for REF (reference panel): ', sizeref)
    print('Set size for IMP (study population): ', sizeimp)
    samples_files = ['cat {}/ALL.chr20.snps.allID.txt '.format(prm.WD + '/gt')
                     + '| head -{} > {}/ALL.chr20.snps.impID.txt'.format(sizeimp, prm.WD + '/gt'),
                     'cat {}/ALL.chr20.snps.allID.txt '.format(prm.WD + '/gt')
                     + '| head -{} | tail -{} > {}/ALL.chr20.snps.refID.txt'.format(sizeimp + sizeref,
                                                                                    sizeref,
                                                                                    prm.WD + '/gt'),
                     'dos2unix {}/ALL.chr20.snps.refID.txt'.format(prm.WD + '/gt'),
                     'dos2unix {}/ALL.chr20.snps.impID.txt'.format(prm.WD + '/gt')]
    for f in samples_files:
        subprocess.run(f, shell=True, cwd=cd)


def partition_imp(tpl: Tuple[str, Dict[str, str]], total_ref: bool = False) -> None:
    print('\n\nREF/IMP SAMPLING'.ljust(80, '.'))
    folder, dic = tpl
    delete_file(os.path.join(folder, dic['imp']))
    delete_file(os.path.join(folder, dic['imp'] + '.csi'))
    pybcf.sampling(dic['gz'],
                   dic['imp'],
                      '{}/ALL.chr20.snps.impID.txt'.format(prm.WD + '/gt'),
                   folder)

    if total_ref and dic.name != 'raw':
        pybcf.rename_samples(dic['imp'], dic['imp'][:-3], folder, '_IMP')
        pybcf.bgzip(dic['imp'][:-3], dic['imp'], folder)
        delete_file(dic['imp'][:-3])

    pybcf.index(dic['imp'], folder)


def partition_ref(dic: dict,  path: str) -> None:
    pybcf.sampling(dic['gz'].replace('.gl', '.gt'),
                   dic['ref'].replace('.gl', '.gt'),
                      '{}/ALL.chr20.snps.refID.txt'.format(prm.WD + '/gt'),
                   path)
    pybcf.index(dic['ref'].replace('.gl', '.gt'), path)


def beagle_haplo_to_geno() -> None:
    if prm.GTGL == 'GL':
        pass
    # elif reference panel is GL encoded:
    # bgl1gtgl = ' '.join(['java -Xss5m -jar {}'.format(prm.BEAGLE_JAR),
    #                      '{}='.format('gtgl') + raw['imp'],
    #                      'impute=false',
    #                      'gprobs=true',
    #                      'out=' + 'temp.' + raw['imp'][:-7],
    #                      '&',
    #                      'java -Xmx4000m -jar {}'.format(prm.BEAGLE_JAR),
    #                      '{}='.format('gtgl') + raw['ref'],
    #                      'impute=false',
    #                      'gprobs=true',
    #                      'out=' + 'temp.' + raw['ref'][:-7]
    #                      ])
    #
    # bgl1gt = ' '.join(['java -Xss5m -jar {}'.format(prm.BEAGLE_JAR),
    #                    '{}='.format('gt')
    #                    + '{}'.format('temp.' if GTGL == 'GL' else '')
    #                    + raw['imp'],
    #                    'impute=false',
    #                    'gprobs=true',
    #                    'out=' + raw['b1i'],
    #                    '&',
    #                    'java -Xmx4000m -jar {}'.format(prm.BEAGLE_JAR),
    #                    '{}='.format('gt')
    #                    + '{}'.format('temp.' if GTGL == 'GL' else '')
    #                    + raw['ref'],
    #                    'impute=false',
    #                    'gprobs=true',
    #                    'out=' + raw['b1r']
    #                    ])
    # if GTGL == 'GL':
    #     subprocess.run(bgl1gtgl, shell=True, cwd=cd)


def beagle_phasing(dic: dict, path_gt_files: str, cd: str) -> None:
    print('\n\nBEAGLE ROUND#1'.ljust(80, '.'))
    os.chdir(cd)
    print('Directory: ', cd)

    if dic.name == 'raw':
        delete_file(dic['b1r'] + '.vcf.gz')
        delete_file(dic['b1i'] + '.vcf.gz')

        bgl1gt = ' '.join(['java -Xss5m -jar {}'.format(prm.BEAGLE_JAR),
                           # -Xmx4000m replaced by -Xss5m because of Exception in thread "main"
                           # java.lang.StackOverflowError at dag.MergeableDag.similar(MergeableDag.java:374)
                           '{}='.format('gt')
                           + os.path.join(path_gt_files, dic['imp'].replace('.gl', '.gt')),
                           'impute=false',
                           'gprobs=true',
                           'out=' + dic['b1i'],
                           '&',
                           'java -Xss5m -jar {}'.format(prm.BEAGLE_JAR),
                           '{}='.format('gt')
                           + os.path.join(path_gt_files, dic['ref'].replace('.gl', '.gt')),
                           'impute=false',
                           'gprobs=true',
                           'out=' + dic['b1r'],
                           # 'map=' + os.path.join(os.path.expanduser('~'),
                           #                       '1000Genomes/data/plink.GRCh37.map/plink.chr20.GRCh37.map')
                           ])

        subprocess.run(bgl1gt, shell=True, cwd=cd)
        pybcf.index(dic['b1r'] + '.vcf.gz', cd)
        pybcf.index(dic['b1i'] + '.vcf.gz', cd)
        delete_file('temp.' + dic['imp'])
        delete_file('temp.' + dic['ref'])

    else:  # dic.name == 'pooled' or 'missing'
        delete_file(dic['b1'] + '.vcf.gz')

        bgl1gtgl = ' '.join(['java -Xss5m -jar {}'.format(prm.BEAGLE_JAR),
                             '{}='.format('gtgl') + dic['imp'],
                             'impute=false',
                             'gprobs=true',
                             'out=' + 'temp.b1',
                             # 'map=' + os.path.join(os.path.expanduser('~'),
                             #                       '1000Genomes/data/plink.GRCh37.map/plink.chr20.GRCh37.map')
                             ])

        bgl1gt = ' '.join(['java -Xss5m -jar {}'.format(prm.BEAGLE_JAR),
                           '{}='.format('gt')
                           + '{}'.format('temp.b1.vcf.gz' if prm.GTGL == 'GL' else dic['imp']),
                           'impute=false',
                           'gprobs=true',
                           'out=' + dic['b1'],
                           # 'map=' + os.path.join(os.path.expanduser('~'),
                           #                       '1000Genomes/data/plink.GRCh37.map/plink.chr20.GRCh37.map')
                           ])

        if prm.GTGL == 'GL':
            subprocess.run(bgl1gtgl, shell=True, cwd=cd)
        subprocess.run(bgl1gt, shell=True, cwd=cd)
        pybcf.index(dic['b1'] + '.vcf.gz', cd)
        delete_file('temp.b1' + '.vcf.gz')


def switch_off_markers(dic: Dict[str, str], cd: str, rm_snp: str) -> str:
    # remove snp for simulating it as ungenotyped
    mkdir('./rm_20:{}'.format(rm_snp))
    path_out = os.path.join(cd, 'rm_20:{}'.format(rm_snp))
    cmd = ' '.join(['bcftools view -e POS={}'.format(rm_snp),
                    '-Oz -o',
                    # os.path.join(path_out, dic['b1'] + '.vcf.gz'),
                    # dic['b1'] + '.vcf.gz'
                    os.path.join(path_out, dic['imp']),
                    dic['imp']
                    ])

    subprocess.run(cmd, shell=True, cwd=cd)

    # pybcf.index(dic['b1'] + '.vcf.gz', path_out)
    pybcf.index(dic['imp'], path_out)

    return path_out


def keep_single_sample(dic: Dict[str, str], cd: str, sample_name: str) -> str:
    mkdir('./keeponly_{}'.format(sample_name))
    path_out = os.path.join(cd, 'keeponly_{}'.format(sample_name))
    cmd = ' '.join(['bcftools view -s {}'.format(sample_name),
                    '-Oz -o',
                    # os.path.join(path_out, dic['b1'] + '.vcf.gz'),
                    # dic['b1'] + '.vcf.gz'
                    os.path.join(path_out, dic['imp']),
                    dic['imp']
                    ])

    subprocess.run(cmd, shell=True, cwd=cd)

    # pybcf.index(dic['b1'] + '.vcf.gz', path_out)
    pybcf.index(dic['imp'], path_out)

    return path_out


def all_snps_all_samples(dic: Dict[str, str], cd: str) -> str:
    mkdir('./all_snps_all_samples')
    path_out = os.path.join(cd, 'all_snps_all_samples')
    cmd = ' '.join(['bcftools view',
                    '-Oz -o',
                    # os.path.join(path_out, dic['b1'] + '.vcf.gz'),
                    # dic['b1'] + '.vcf.gz'
                    os.path.join(path_out, dic['imp']),
                    dic['imp']
                    ])

    subprocess.run(cmd, shell=True, cwd=cd)

    # pybcf.index(dic['b1'] + '.vcf.gz', path_out)
    pybcf.index(dic['imp'], path_out)

    return path_out


def conform_gt(dic: dict, dicraw: dict, cd: str) -> bool:
    print('\n\nCONFORM-GT'.ljust(80, '.'))
    # GT for reference files, even when working with GLs
    os.chdir(cd)
    print(os.getcwd())
    delete_file(dic['cfgt'] + '.vcf.gz')

    cfgt = ' '.join(['java -jar {}'.format(prm.CFGT_JAR),
                     '{}='.format('gt') + dic['b1'] + '.vcf.gz',
                     'chrom=20:60343-62965354',
                     'ref={}'.format(os.path.join(cd,
                                                  dicraw['b1r'] + '.vcf.gz')),
                     'out=' + dic['cfgt']
                     ])
    try:
        subprocess.run(cfgt, shell=True, cwd=cd)
        assert os.path.exists(dic['cfgt'] + '.vcf.gz')
    except AssertionError:  # if duplicated markers, just copy phased file
        shutil.copy(dic['b1'] + '.vcf.gz', dic['cfgt'] + '.vcf.gz')
    pybcf.index(dic['cfgt'] + '.vcf.gz', cd)

    return check_file_creation(cd, dic['cfgt'] + '.vcf.gz')


def beagle_imputing(dic_study: dict, dicref: dict, cd: str) -> bool:
    print('\n\nBEAGLE (ROUND#2)'.ljust(80, '.'))
    os.chdir(cd)
    delete_file(dic_study['b2'] + '.vcf.gz')

    bgl2 = ' '.join(['java -Xss5m -jar {}'.format(prm.BEAGLE_JAR),  # -Xss5m option: fix StackOverFlow error of Java
                     'gt=' + dic_study['cfgt'] + '.vcf.gz',
                     'ref={}'.format(os.path.join(cd,
                                                  dicref['b1r'] + '.vcf.gz')),
                     'impute=true',
                     'gprobs=true',
                     'out=' + dic_study['b2'],
                     # 'map=' + os.path.join(os.path.expanduser('~'),
                     #                       '1000Genomes/data/plink.GRCh37.map/plink.chr20.GRCh37.map')
                     ])

    subprocess.run(bgl2, shell=True, cwd=cd)
    pybcf.index(dic_study['b2'] + '.vcf.gz', cd)

    return check_file_creation(cd, dic_study['b2'] + '.vcf.gz')


def reformat_fields(dic_study: dict, cd: str) -> bool:
    print('\n\nREFORMATTING GP AND DS FIELDS'.ljust(80, '.'))
    os.chdir(cd)
    delete_file(dic_study['corr'] + '.vcf.gz')

    refmt = ' '.join(["bcftools view {}.vcf.gz".format(dic_study['b2']),
                      "| sed 's/##FORMAT=<ID=DS,Number=A,Type=Float/##FORMAT=<ID=DS,Number=1,Type=String/'",
                      "| sed 's/##FORMAT=<ID=GP,Number=G,Type=Float/##FORMAT=<ID=GP,Number=3,Type=String/'",
                      "| bcftools view -Oz -o {}.vcf.gz".format(dic_study['corr'])
                      ])

    subprocess.run(refmt, shell=True, cwd=cd)
    pybcf.index(dic_study['corr'] + '.vcf.gz', cd)

    gtonly = ' '.join(["bcftools annotate -x 'FORMAT'",
                       dic_study['corr'] + '.vcf.gz',
                       "| bgzip >",
                       dic_study['gtonly'] + '.vcf.gz'])

    subprocess.run(gtonly, shell=True, cwd=cd)
    pybcf.index(dic_study['gtonly'] + '.vcf.gz', cd)

    return check_file_creation(cd, dic_study['gtonly'] + '.vcf.gz')


def clean_imputed_directory(cd: str) -> bool:
    os.chdir(cd)
    print('\n\nCLEANING DIRECTORY {}'.format(os.getcwd()).ljust(80, '.'))
    for f in os.scandir(cd):
        if f.is_file()and (f.path.endswith('.log') or '.cfgt.' in f.path):
            delete_file(f)

    return True


def move_file_to(f_in: str, f_out: str, cd: str):
    # cd: path to out directory
    cmd = ' '.join(['bcftools view',
                    '-Oz -o',
                    f_out,
                    f_in
                    ])
    print(cmd)
    subprocess.run(cmd, shell=True, cwd=cd)
    pybcf.index(f_out, cd)

    return check_file_creation(os.getcwd(), f_out)


def merge_files(pattern: str, f_out: str, cd: str):
    mkdir(os.path.join(cd, 'single_samples_merged'))
    os.chdir(os.path.join(cd, 'single_samples_merged'))
    # with open('files2merge.txt', mode='w+', encoding='utf-8') as f:
    #     for line in flist:
    #         f.write(line)
    #         f.write('\n')
    # files = ' '.join(flist)
    cmd = ' '.join(['bcftools merge',
                    pattern,
                    '-Oz -o',
                    f_out,
                    #files
                    ])
    print(cmd)
    subprocess.run(cmd, shell=True, cwd=os.getcwd())
    pybcf.index(f_out, os.getcwd())

    return check_file_creation(os.getcwd(), f_out)


def concat_files(flist: list, f_out: str, dic: dict, cd: str):
    rechr20pos = re.compile(r'(20\:)(\d+)')
    poslist: list = []
    for f in flist:
        dirpath = os.path.dirname(f)
        os.chdir(os.path.join(cd, dirpath))
        pos = re.search(rechr20pos, f)[0]
        poslist.append(pos)

        delete_file(f + '.csi')
        variant = ' '.join(['bcftools view',
                            '-i POS={}'.format(pos[3:]),  # trim 20: from the regex
                            '-Oz -o',
                            'marker_{}.vcf.gz'.format(pos),
                            dic['gtonly'] + '.vcf.gz'
                            ])
        print(variant)
        subprocess.run(variant, shell=True, cwd=os.getcwd())
        pybcf.index(f, os.getcwd())

    mkdir(os.path.join(cd, 'ko_markers_merged'))
    os.chdir(os.path.join(cd, 'ko_markers_merged'))

    f_init = ' '.join(['bcftools view',
                       '-t ^' + ','.join(poslist),
                       '-Oz -o',
                       'ko_markers_core.vcf.gz',
                       os.path.join(os.path.dirname(flist[0]),
                                    dic['gtonly'] + '.vcf.gz')
                       ])
    print(f_init)
    subprocess.run(f_init, shell=True, cwd=os.getcwd())
    pybcf.index('ko_markers_core.vcf.gz', os.getcwd())

    files = ' '.join(flist)
    cmd = ' '.join(['bcftools concat -a',
                    '-Oz -o',
                    f_out,
                    'ko_markers_core.vcf.gz',
                    files
                    ])
    print(cmd)
    subprocess.run(cmd, shell=True, cwd=os.getcwd())
    pybcf.index(f_out, os.getcwd())

    return check_file_creation(os.getcwd(), f_out)
