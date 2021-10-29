#!/bin/bash

# this scripts sorts the sample in input file accordingly to the samples ID in text file

infile=$1
samplesId=$2

bcftools view -S $samplesId -Oz -o tmp.sorted.vcf.gz $infile
bcftools view -Oz -o $infile tmp.sorted.vcf.gz
bcftools index -f $infile
rm tmp.sorted.vcf.gz 

