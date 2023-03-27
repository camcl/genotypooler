#!/bin/bash

beaglejar=../bin/beagle.11Mar19.69c.jar
cfgtjar=../bin/conform-gt.jar

modelscale=0.8

chrom=$( bcftools query -f '%CHROM\n' Chr1.Founders.vcf.gz | head -1 )
startpos=$( bcftools query -f '%POS\n' Chr1.Founders.vcf.gz | head -1 )
endpos=$( bcftools query -f '%POS\n' Chr1.Founders.vcf.gz | tail -1 )

echo 'Contigs in the reference file'
echo '.................................................................................'
echo 'Chromosome ' $chrom '   Startpos =' $startpos '   Endpos =' $endpos 
echo ''

echo ''
echo 'Check FORMAT field in files for imputation'
echo '.................................................................................'
reftruefmt=$( bcftools query -f '%LINE\n' Chr$chrom.Founders.vcf.gz | head -1 | cut -f9 )
echo 'FORMAT in reference panel: ' $reftruefmt
imppoolfmt=$( bcftools query -f '%LINE\n' Chr$chrom.SNPs.pruned.nomiss.pooled.vcf.gz | head -1 | cut -f9 )
echo 'FORMAT in target: ' $imppoolfmt
echo ''

echo ''
echo 'Check number of samples and number of markers in files for imputation'
echo '.................................................................................'
echo 'reference:'
bcftools query -l Chr$chrom.Founders.vcf.gz | wc -l
#bcftools view -H  Chr$chrom.Founders.vcf.gz | wc -l
echo ''
echo 'target:'
bcftools query -l Chr$chrom.SNPs.pruned.nomiss.pooled.vcf.gz | wc -l
#bcftools view -H Chr$chrom.SNPs.pruned.nomiss.pooled.vcf.gz | wc -l
echo ''

echo ''
echo 'Phase reference and target with BEAGLE'
echo '.................................................................................'
echo 'Beagle .jar file used at:' $beaglejar
echo ''
#java -Xss5m -jar $beaglejar impute=false gtgl=Chr$chrom.Founders.vcf.gz out=Chr$chrom.Founders.phased
cp Chr$chrom.Founders.vcf.gz Chr$chrom.Founders.phased.vcf.gz
bcftools index -f Chr$chrom.Founders.phased.vcf.gz
refphasfmt=$( bcftools query -f '%LINE\n' Chr$chrom.Founders.phased.vcf.gz | head -1 | cut -f9 )
echo 'FORMAT in the phased ref file:' $refphasfmt 

### two runs for phasing from GL: 1st run from GL to GT:DP:GP (unphased GT)
java -Xss5m -jar $beaglejar impute=false gtgl=Chr$chrom.SNPs.pruned.nomiss.pooled.vcf.gz map=1_interpolated_wheat_map_plink.map modelscale=$modelscale out=Chr$chrom.SNPs.pruned.nomiss.pooled.unphased
### two runs for phasing from GL: 2nd run from GT:DP:GP to phased GT
### with gt argument, all genotypes in the output file will be phased and non-missing  (Beagle4.1 documentation)
echo ''
java -Xss5m -jar $beaglejar impute=false gt=Chr$chrom.SNPs.pruned.nomiss.pooled.unphased.vcf.gz map=1_interpolated_wheat_map_plink.map modelscale=$modelscale out=Chr$chrom.SNPs.pruned.nomiss.pooled.phased
bcftools index -f Chr$chrom.SNPs.pruned.nomiss.pooled.phased.vcf.gz
impphasfmt=$( bcftools query -f '%LINE\n' Chr$chrom.SNPs.pruned.nomiss.pooled.phased.vcf.gz | head -1 | cut -f9 )
echo 'FORMAT in the phased target file:' $impphasfmt
echo ''

echo ''
echo 'Deduplicate possibly duplicated markers'
echo '.................................................................................'
bcftools norm --rm-dup all -Oz -o Chr$chrom.SNPs.pruned.nomiss.pooled.phased.dedup.vcf.gz Chr$chrom.SNPs.pruned.nomiss.pooled.phased.vcf.gz
bcftools index -f Chr$chrom.SNPs.pruned.nomiss.pooled.phased.dedup.vcf.gz
bcftools norm --rm-dup all -Oz -o Chr$chrom.Founders.phased.dedup.vcf.gz Chr$chrom.Founders.phased.vcf.gz
bcftools index -f Chr$chrom.Founders.phased.dedup.vcf.gz
echo ''

echo ''
echo 'Unify reference and target markers with CONFORM-GT'
echo '.................................................................................'
echo 'conform-gt .jar file used at:' $cfgtjar
### necessary to get proper imputation results and GT:DS:GP fields with gprobs=true
echo ''
java -jar $cfgtjar ref=Chr$chrom.Founders.phased.dedup.vcf.gz gt=Chr$chrom.SNPs.pruned.nomiss.pooled.phased.dedup.vcf.gz chrom=$chrom:$startpos-$endpos out=Chr$chrom.SNPs.pruned.nomiss.pooled.cfgt
bcftools index -f Chr$chrom.SNPs.pruned.nomiss.pooled.cfgt.vcf.gz
echo ''

echo ''
echo 'Impute target from reference with BEAGLE'
echo '.................................................................................'
echo 'Beagle .jar file used at:' $beaglejar
### impute=true must be used with gt= and ref= (Beagle4.1 documentation)
echo ''
### map=$chrom_interpolated_wheat_map_plink.map does not work (why???), looks like it has to be hard cpded
java -Xss5m -jar $beaglejar gt=Chr$chrom.SNPs.pruned.nomiss.pooled.cfgt.vcf.gz ref=Chr$chrom.Founders.phased.vcf.gz map=1_interpolated_wheat_map_plink.map modelscale=$modelscale impute=true gprobs=true out=Chr$chrom.SNPs.pruned.nomiss.pooled.imputed
bcftools index -f Chr$chrom.SNPs.pruned.nomiss.pooled.imputed.vcf.gz
impimpfmt=$( bcftools query -f '%LINE\n' Chr$chrom.SNPs.pruned.nomiss.pooled.imputed.vcf.gz | head -1 | cut -f9 )
echo 'FORMAT in the imputed target file:' $impimpfmt
echo''

echo ''
echo 'Cleaning directory from log files'
echo '.................................................................................'
rm ./*.log
rm ./*cfgt.vcf.gz* # files from conform-gt cannot be overwritten
echo 'Done.'
echo ''
