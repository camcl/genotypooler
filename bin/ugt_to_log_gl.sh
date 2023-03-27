#!/bin/bash

# Converts a VCF file written with unphased GT to a VCF file with log GL
# Usage example: 
# $ bash /home/camille/MagicWheat/src/genotypooler/bin/ugt_to_log_gl.sh file.gt.vcf.gz file.gl.vcf.gz

fin=$1
fout=$2
fbname=$(basename "$fout" .vcf.gz)

bcftools view $fin | sed 's/##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">/##FORMAT=<ID=GL,Number=G,Type=Float,Description="Estimated Genotype Probability">/' | sed 's/GT/GL/g' | sed -e 's/1\/1/-12.0,-12.0,0.0/g' -e 's/0\/0/0.0,-12.0,-12.0/g' -e 's/1\/0/-12.0,0.0,-12.0/g' -e 's/0\/1/-12.0,0.0,-12.0/g' -e 's/\.\/\./-0.48148606,-0.48148606,-0.48148606/g' > $fbname.vcf
bcftools view -Oz -o $fout $fbname.vcf
bcftools index -f $fout
rm $fbname.vcf


