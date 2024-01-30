#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os, sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import pysam

from typing import *

rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from genotypooler.poolSNPs import pybcf
from genotypooler.poolSNPs import dataframe as vcfdf
from genotypooler.persotools.files import *

import img2pdf
from PIL import Image, ImageDraw, ImageFont, ImageColor

"""
For a given variant, plot the pooling blocks as coloured square matrices with 
genotypes displayed each with RGB color proportional to GL(RR|RA|AA).
Useful to get a picture of how genotypes are processed from the true ones into pooled and finally imputed ones   
"""

### COMMAND-LINE PARSING AND PARAMETERS (arguments are parsed from a file
parser = argparse.ArgumentParser(description='Plots the pooling blocks as colored squares')
parser.add_argument('cycleN', metavar='cycleN', type=str, help='Cycle', default=None)
argsin = parser.parse_args()

# TODO: docstrings and annotations

'''
Human data (15 blocks 4*4 per variant):
20:264365    0.512181  --> pools #1, #2, #4 --> row 1
20:62915126  0.006190  --> pools #3, #4, #14 relevant with ALT carrier --> row 1

Wheat data (31 blocks 4*4 per variant):
1:11814388 (same as in TAG and PLoS manuscripts
'''

# chr, pos = '20', '264365'  # human
chr, pos = '1', '11814388'
myvar = '{}:{}'.format(chr, pos)
cycle = argsin.cycleN

print(f'\nCycle {cycle}\r'.ljust(80, '.'))


# Paths to files to read genotypes from

# outdir = '/home/camille/PoolImpHuman/results/20200827' # human
outdir = f'/home/camille/IterDecodeImpute/runs/poolimputeSNPs/results/data/1/20230608/cycle{cycle}'  # wheat
if not os.path.exists(outdir):
    os.mkdir(outdir)

# truef = '/home/camille/PoolImpHuman/data/20200827/IMP.chr20.snps.gt.vcf.gz'  # human
# pooledf = '/home/camille/PoolImpHuman/data/20200827/IMP.chr20.pooled.snps.gl.vcf.gz'  # human
# imputedf = '/home/camille/PoolImpHuman/data/20200827/IMP.chr20.pooled.imputed.vcf.gz'  # human
truef = '/home/camille/IterDecodeImpute/runs/poolimputeSNPs/results/data/1/STU.Chr1.SNPs.pruned.sorted.vcf.gz'
pooledf = f'/home/camille/IterDecodeImpute/runs/poolimputeSNPs/results/data/1/20230608/cycle{cycle}/STU.Chr1.SNPs.pruned.sorted.pooled.vcf.gz'
imputedf = f'/home/camille/IterDecodeImpute/runs/poolimputeSNPs/results/data/1/20230608/cycle{cycle}/STU.Chr1.SNPs.pruned.sorted.pooled.imputed.vcf.gz'

for f in [truef, pooledf, imputedf]:
    pybcf.index(f, os.path.dirname(f))


# Create files for the chosen SNP only (speed up processing then

truevar = pybcf.view_one_variant(truef, pos, chr, wd=outdir)
pooledvar = pybcf.view_one_variant(pooledf, pos, chr, wd=outdir)
imputedvar = pybcf.view_one_variant(imputedf, pos, chr, wd=outdir)


# Colouring and labelling functions

def gl_colors(gen: np.ndarray, min_log_gl: float = -12.) -> np.ndarray:
    """
    Convert log-GL values to not logged ones, used as a RGB color triplet
    :param gen: genotype (log-)likelihoods. Array-like.
    :param min_log_gl: float (negative) indicating the minimum value used instead of log(GP=0.0)
    :return: triplet coding RGB channels
    """
    gen = np.array([g for g in gen])  # forces GL-tuple to GL-array conversion
    tenpow = np.vectorize(lambda x: 0.0 if x == min_log_gl else pow(10.0, x))
    rgb = tenpow(gen)
    return rgb


def map_gt_gl(gt_val):
    # gt_val from trinary encoding, true or imputed VCF file
    if gt_val == 2.0:
        gl_val = np.array([0., 0., 1.])
    elif gt_val == 1.0:
        gl_val = np.array([0., 1., 0.])
    elif gt_val == 0.0:
        gl_val = np.array([1., 0., 0.])
    else:  # missing data
        gl_val = np.array([1/3, 1/3, 1/3])

    return gl_val


def gt_colors(gen):
    """

    :param gen: raw genotypes with phase. Array-like.
    :return: triplet coding rgb
    """

    rgb = np.apply_along_axis(map_gt_gl, axis=-1, arr=gen)
    return rgb


def gl_labels(gen):
    """
    Returns the most likely genotype (slots labels of the pooling blocks).
    :param gen: triplet for genotype likelihoods
    :return:
    """
    if gen[0] == 1.0:
        gl = '0/0'
    elif gen[2] == 1.0:
        gl = '1/1'
    elif gen[1] == 1.0:
        gl = '1/0'
    elif round(gen[1], 1) == 0.5 and round(gen[2], 1) == 0.5:
        gl = '1/.'
    elif round(gen[0], 1) == 0.5 and round(gen[1], 1) == 0.5:
        gl = '0/.'
    else:
        gl = './.'
    return gl


def g_gl_labels(gen):
    """
    Returns the most likely genotype (slots labels of the pooling blocks).
    :param gen: triplet for genotype likelihoods
    :return:
    """
    if gen[0] == 1.0:
        g = '0'
    elif gen[2] == 1.0:
        g = '2'
    elif gen[1] == 1.0:
        g = '1'
    else:
        g = '-1'
    return g


def gt_labels(gen):
    """
    Returns the label genotype (slots labels of the pooling blocks).
    :param gen: hexa or trinary encoding
    :return: str: GT representation
    """
    if gen == 0.0:
        gt = '0/0'
    elif gen == 2.0:
        gt = '1/1'
    elif gen == 1.0:
        gt = '1/0'
    elif gen == 0.5:
        gt = '1/.'
    elif gen == -0.5:
        gt = '0/.'
    else:
        gt = './.'
    return gt


def g_gt_labels(gen):
    """
    Returns the label genotype (slots labels of the pooling blocks).
    :param gen: hexa or trinary encoding
    :return: str: G representation
    """
    if gen == 0.0:
        g = '0'
    elif gen == 2.0:
        g = '2'
    elif gen == 1.0:
        g = '1'
    else:
        g = '-1'
    return g


def gt_gl_to_rgb(gtgl, var):
    """
    Compute custom RGB color-channels and labels for genotypes.
    """
    n_pop_blocks = len(var) // (4 * 4)  # /!\ 4 ' 4 is hardcoded!
    print('\nn_pop_blocks = ', n_pop_blocks)
    pool_rgb = np.zeros((var.shape[0], 3), dtype=float)  # 3 color channels
    p_labs = np.empty(var.shape[0], dtype='<U3')  # 'x/x' is 3-character long
    if gtgl.upper() == 'GL':
        pool_rgb = np.apply_along_axis(gl_colors, -1, var)  # gl_colors(var[:, np.newaxis])
        for idx, g in enumerate(pool_rgb):  # unlogged GL colors
            p_labs[idx] = g_gl_labels(g)
        pool_rgb = pool_rgb.reshape(n_pop_blocks, 16, 3)  # block-wise reshaping
        p_labs = p_labs.reshape(n_pop_blocks, 16)
    if gtgl.upper() == 'GT':
        pool_rgb = gt_colors(var.values[:, np.newaxis])
        for idx, g in enumerate(var):
            p_labs[idx] = g_gt_labels(g)
        pool_rgb = pool_rgb.reshape(n_pop_blocks, 16, 3)  # block-wise reshaping
        p_labs = p_labs.reshape(n_pop_blocks, 16)

    pool_colors = np.multiply(pool_rgb, np.broadcast_to([1.0, 1.0, 1.0], pool_rgb.shape))  # [0.5, 0.8, 0.5]

    return pool_rgb, pool_colors, p_labs


def plot_n_first_blocks(pool_rgb, pool_colors, pool_labels, step, snp, nb_blocks, af_info=None, aaf=None):
    """
    First nb_blocks pools only
    :param var: trinary encoded genotypes or log GL
    :return:
    """
    if af_info is None:
        af_info = np.nan
    if aaf is None:
        aaf = np.nan

    plt.rcParams['axes.titlesize'] = 12
    plots_sz = (1, nb_blocks)  # plot nb_blocks first blocks as 1 row
    fig, axes = plt.subplots(plots_sz[0], plots_sz[1],
                             figsize=(8, 2.5),
                             subplot_kw={'xticks': [], 'yticks': []})
    fig.suptitle('{0} data (SNP {1}; AAF = {2:.5f})'.format(step.capitalize(),
                                                            snp,
                                                            aaf), fontsize=16, fontweight='bold', va='top')

    k = 0
    for j in range(1, plots_sz[1] + 1):
        k += 1
        colors = np.ceil(pool_colors[k-1, :, :].reshape((4, 4, 3)))
        alphas = 1.0 + 2 * np.log10(np.max(pool_colors[k-1, :, :].reshape((4, 4, 3)), axis=-1))  # custom transparency value
        rgba = np.dstack([colors, alphas])
        if step == 'true':
            rgba = np.ones_like(rgba)  # white color everywhere, no transparency
        axes[j - 1].imshow(rgba,
                           cmap='plasma')
        axes[j - 1].set_title('Block #{}'.format(k), pad=2, fontsize=16)
        axes[j - 1].set_xticks(np.arange(4) + 0.5, minor=True)
        axes[j - 1].set_yticks(np.arange(4) + 0.5, minor=True)
        axes[j - 1].tick_params(which="both",
                                axis='both',
                                bottom=False,
                                left=False,
                                labelsize=8,
                                length=1,
                                pad=0.3)
        axes[j - 1].grid(which='minor', axis='both',
                         color="k" if step == 'true' else "w",
                         linestyle='-', linewidth=0.25)
        # remnove borders
        if step != 'true':
            axes[j - 1].spines['top'].set_visible(False)
            axes[j - 1].spines['right'].set_visible(False)
            axes[j - 1].spines['bottom'].set_visible(False)
            axes[j - 1].spines['left'].set_visible(False)
        tx_i_j = pool_labels[k - 1, :].reshape((4, 4))
        for m in range(4):
            for n in range(4):
                axes[j - 1].text(n, m, tx_i_j[m, n], ha="center", va="center",
                                 color="k" if step == 'true' else "k", fontsize=16)
    plt.autoscale()
    fig.tight_layout()
    plt.savefig(os.path.join(outdir, 'pools-patterns-{}-snp{}.jpg'.format(step, snp)),
                dpi=500)


def gt_gl_to_rg(genos_gt, genos_gl):
    """
    Color in green pooled or imputed predictions that match the true genotype, otherwise color in red.
    The color shade is proportional to Pr(G).
    """
    genos_rg = []
    dim3 = genos_gt.shape
    for gt, gl in zip(genos_gt.reshape((dim3[0]*dim3[1], 3)), genos_gl.reshape((dim3[0]*dim3[1], 3))):
        if np.argmax(gt) == np.argmax(gl):
            genos_rg.append(np.array([0.0, np.max(gl), 0.0]))
        else:
            genos_rg.append(np.array([np.max(gl), 0.0, 0.0]))

    return np.array(genos_rg).reshape(dim3)


def get_concat_v(img_list, color=(256, 256, 256)):
    """
    For tet: see examples at https://holypython.com/python-pil-tutorial/how-to-add-text-to-images-in-python-via-pil-library/
    """
    offset_top = 500
    img_files = [Image.open(img_path) for img_path in img_list]
    dst = Image.new('RGB', (img_files[0].width, offset_top + sum([img.height for img in img_files])), color)
    dst.paste(img_files[0], (0, offset_top))
    for i in range(1, len(img_files)):
        dst.paste(img_files[i], (0, offset_top + sum([img.height for img in img_files[:i]])))
    # add title for the concatenated jpg images (List of fonts in Ubuntu: $ fc-list
    font = ImageFont.truetype('./liberation/LiberationSans-Regular.ttf', 182)
    draw = ImageDraw.Draw(dst)
    draw.text((dst.width//2 - 400, offset_top//4), f'Cycle {cycle}', fill='black', font=font)

    return dst


def before_after_pooling(snp):
    """
    Combine pools overview in a PDF document
    :param snp:
    :param sz:
    :return:
    """
    imglist = [os.path.join(outdir, 'pools-patterns-true-snp{}.jpg'.format(snp)),
               os.path.join(outdir, 'pools-patterns-pooled-snp{}.jpg'.format(snp)),
               os.path.join(outdir, 'pools-patterns-imputed-snp{}.jpg'.format(snp))] # imputed data at top (more visible)
    with open(os.path.join(outdir, 'pools-patterns-snp{}.pdf'.format(snp.replace(':', '-'))), "wb") as f:
        f.write(img2pdf.convert([i for i in imglist], dpi=300, title='SNP {}'.format(snp)))
    get_concat_v(imglist).save(os.path.join(outdir, 'pools-patterns-snp{}-cycle{}.jpg'.format(snp.replace(':', '-'),
                                                                                              cycle.zfill(2))))


# Read files, extract relevant informations for plotting and plot pooling blocks

dftrue = vcfdf.PandasMixedVCF(truevar, format='GT', indextype='chrom:pos')
dfpooled = vcfdf.PandasMixedVCF(pooledvar, format='GL', indextype='chrom:pos')
dfimputed = vcfdf.PandasMixedVCF(imputedvar, format='GT', indextype='chrom:pos')

genos_true = dftrue.trinary_encoding().loc[myvar]
# af_info = float(dftrue.af_info.loc[myvar].values)
aaf_true = float(dftrue.aaf.loc[myvar].values)

genos_pooled = dfpooled.genotypes().loc[myvar].values  # .trinary_encoding().loc[myvar]
max_pooled = [np.where(gl_colors(gen) > 0.5, gl_colors(gen), 0.0) for gen in genos_pooled]
# print(np.where(np.sum(max_pooled, axis=-1) > 0.0, max_pooled, np.array([[1.0, 1.0, 1.0]])))
aaf_pooled = float(dfpooled.aaf.loc[myvar].values)

genos_imputed = dfimputed.trinary_encoding().loc[myvar]
aaf_imputed = float(dfimputed.aaf.loc[myvar].values)

true_rgb, true_colors, true_labels = gt_gl_to_rgb('GT', genos_true)
plot_n_first_blocks(true_rgb, true_colors, true_labels, 'true', myvar, 5, af_info=None, aaf=aaf_true)

pool_rgb, pool_colors, pool_labels = gt_gl_to_rgb('GL', genos_pooled)
pooled_rg = gt_gl_to_rg(true_colors, pool_colors)
plot_n_first_blocks(pool_rgb, pooled_rg, pool_labels, 'pooled', myvar, 5, af_info=None, aaf=aaf_pooled)

# imputed_colors = plot_n_first_blocks(genos_imputed, 'GT', 'imputed', myvar, 5, af_info=None, aaf=aaf_imputed)
imputed_rgb, imputed_colors, imputed_labels = gt_gl_to_rgb('GT', genos_imputed)
imputed_rg = gt_gl_to_rg(true_colors, imputed_colors)
plot_n_first_blocks(imputed_rgb, imputed_rg, imputed_labels, 'imputed', myvar, 5, af_info=None, aaf=aaf_imputed)
#print(gt_gl_to_rg(true_colors, imputed_colors))
before_after_pooling(myvar)
