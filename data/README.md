## Description

Marker coordinates from the LD and HD Illumina bead chips intersecting the 1KGP chr20 markers:
* LD only (Illumina Infinium OmniExpress, 35682 markers),
* HD only (Illumina Infinium Omni2.5, 17015 markers).


Lookup file for decoding the genotype of the pools into individual GP format:
* adaptive_gls.csv


Marker map for the chromosome 20 based on the GRCh37 assembly:
* plink.ch20.GRCh37.map
Ths map is to be used for imputation.


ID of samples in the reference panel (REF) and study population (IMP):
* reference.panel (2264 samples, 1 sample per line)
* study.population (240 samples, 1 sample per line)

VCF files in GT format for the reference and the study population with genotypes for all HD markers (i.e. LD only + HD only:
* IMP.chr20.snps.gt.vcf.gz (and .csi index file)
* REF.chr20.snps.gt.vcf.gz (and .csi index file)
The files are created from the 1KGP chromosome 20 phase 3 v5 file.