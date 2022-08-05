#!/bin/bash

beaglejar=../bin/beagle.11Mar19.69c.jar
cfgtjar=../bin/conform-gt.jar

chrom=$( bcftools query -f '%CHROM\n' REF.chr20.snps.gt.vcf.gz | head -1 )
startpos=$( bcftools query -f '%POS\n' REF.chr20.snps.gt.vcf.gz | head -1 )
endpos=$( bcftools query -f '%POS\n' REF.chr20.snps.gt.vcf.gz | tail -1 )

echo 'Contigs in the reference file'
echo '.................................................................................'
echo 'Chromosome ' $chrom '   Startpos =' $startpos '   Endpos =' $endpos 
echo ''

echo ''
echo 'Check FORMAT field in files for imputation'
echo '.................................................................................'
reftruefmt=$( bcftools query -f '%LINE\n' REF.chr20.snps.gt.vcf.gz | head -1 | cut -f9 )
echo 'FORMAT in reference panel: ' $reftruefmt
imppoolfmt=$( bcftools query -f '%LINE\n' IMP.chr20.pooled.snps.gl.vcf.gz | head -1 | cut -f9 )
echo 'FORMAT in target: ' $imppoolfmt
echo ''

echo ''
echo 'Check number of samples and number of markers in files for imputation'
echo '.................................................................................'
echo 'reference:'
bcftools query -l REF.chr20.snps.gt.vcf.gz | wc -l
#bcftools view -H  REF.chr20.snps.gt.vcf.gz | wc -l
echo ''
echo 'target:'
bcftools query -l IMP.chr20.pooled.snps.gl.vcf.gz | wc -l
#bcftools view -H IMP.chr20.pooled.snps.gl.vcf.gz | wc -l
echo ''

echo ''
echo 'Phase reference and target with BEAGLE'
echo '.................................................................................'
echo 'Beagle .jar file used at:' $beaglejar
echo ''
#java -Xss5m -jar $beaglejar impute=false gtgl=REF.chr20.snps.gt.vcf.gz out=REF.chr20.phased
cp REF.chr20.snps.gt.vcf.gz REF.chr20.phased.vcf.gz
bcftools index -f REF.chr20.phased.vcf.gz
refphasfmt=$( bcftools query -f '%LINE\n' REF.chr20.phased.vcf.gz | head -1 | cut -f9 )
echo 'FORMAT in the phased ref file:' $refphasfmt 

### two runs for phasing from GL: 1st run from GL to GT:DP:GP (unphased GT)
java -Xss5m -jar $beaglejar impute=false gtgl=IMP.chr20.pooled.snps.gl.vcf.gz out=IMP.chr20.pooled.unphased
### two runs for phasing from GL: 2nd run from GT:DP:GP to phased GT
### with gt argument, all genotypes in the output file will be phased and non-missing  (Beagle4.1 documentation)
echo ''
java -Xss5m -jar $beaglejar impute=false gt=IMP.chr20.pooled.unphased.vcf.gz out=IMP.chr20.pooled.phased
bcftools index -f IMP.chr20.pooled.phased.vcf.gz
impphasfmt=$( bcftools query -f '%LINE\n' IMP.chr20.pooled.phased.vcf.gz | head -1 | cut -f9 )
echo 'FORMAT in the phased target file:' $impphasfmt
echo ''

echo ''
echo 'Deduplicate possibly duplicated markers'
echo '.................................................................................'
bcftools norm --rm-dup all -Oz -o IMP.chr20.pooled.phased.dedup.vcf.gz IMP.chr20.pooled.phased.vcf.gz
bcftools index -f IMP.chr20.pooled.phased.dedup.vcf.gz
bcftools norm --rm-dup all -Oz -o REF.chr20.phased.dedup.vcf.gz REF.chr20.phased.vcf.gz
bcftools index -f REF.chr20.phased.dedup.vcf.gz
echo ''

echo ''
echo 'Unify reference and target markers with CONFORM-GT'
echo '.................................................................................'
echo 'conform-gt .jar file used at:' $cfgtjar
### necessary to get proper imputation results and GT:DS:GP fields with gprobs=true
echo ''
java -jar $cfgtjar ref=REF.chr20.phased.dedup.vcf.gz gt=IMP.chr20.pooled.phased.dedup.vcf.gz chrom=$chrom:$startpos-$endpos out=IMP.chr20.pooled.cfgt
bcftools index -f IMP.chr20.pooled.cfgt.vcf.gz
echo ''

echo ''
echo 'Impute target from reference with BEAGLE'
echo '.................................................................................'
echo 'Beagle .jar file used at:' $beaglejar
### impute=true must be used with gt= and ref= (Beagle4.1 documentation)
echo ''
java -Xss5m -jar $beaglejar gt=IMP.chr20.pooled.cfgt.vcf.gz ref=REF.chr20.phased.vcf.gz impute=true gprobs=true out=IMP.chr20.pooled.imputed
bcftools index -f IMP.chr20.pooled.imputed.vcf.gz
impimpfmt=$( bcftools query -f '%LINE\n' IMP.chr20.pooled.imputed.vcf.gz | head -1 | cut -f9 )
echo 'FORMAT in the imputed target file:' $impimpfmt
echo''

echo ''
echo 'Cleaning directory from log files'
echo '.................................................................................'
rm ./*.log
rm ./*cfgt.vcf.gz* # files from conform-gt cannot be overwritten
echo 'Done.'
echo ''
