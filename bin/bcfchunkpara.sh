#!/bin/bash

# Usage example:
# $ cd genotypooler/data
# $ bash ../bin/bcfchunkpara.sh IMP.chr20.snps.gt.vcf.gz ./tmp 1000
# NB: file_in parameter cannot have path prefix, it must be a file name only

filein=$1
pathout=$2
echo 'Sorting file' $filein
bcftools view -Oz -o tmp.$filein $filein
rm $filein*
bcftools index -f tmp.$filein
bcftools sort -Oz -o $filein tmp.$filein
rm tmp.$filein*
# Remove PR INFO field that causes error
bcftools view -Oz -o tmp.$filein $filein
bcftools annotate -x INFO/PR -Oz -o $filein tmp.$filein
bcftools index -f $filein
rm tmp.$filein*
echo 'Counting lines in' $filein
nblin=$(bcftools query -f '%ID\n' $filein | wc -l)
chksz=$3
div=$(( $nblin/$chksz )) 
floor=$(echo $div | cut -f1 -d ".")
nbchk=$(( $floor + 1 ))
echo 'Number of files to pack = ' $nbchk

mkdir tmp
chrom=$(tabix -l $filein)

# array of assets, assuming at least 1 item exists
listChunks=( { 0..$nbchk } ) # example: 0 1 2 3

# replace with your task
task() { # $1 = idWorker, $2 = chunki
i=$2
startlin=$(( $nblin-$i*$chksz ))
if (( $startlin >= 0 ))
then
	startpos=$(bcftools query -f '%POS\n' $filein | tail -$startlin | head -1)
	endpos=$(bcftools query -f '%POS\n' $filein | tail -$startlin | head -$chksz | tail -1)
	if [[ ($endpos != '' && $startpos != '') ]]
	then
		echo "Worker $1: Packing and writing chunk" $i
		echo 'Starts at POS' $startpos 'and ends at POS' $endpos

		bcftools view -r $chrom:$startpos-$endpos -Oz -o $pathout/pack$i.$filein $filein

		echo "    Worker $1: Chunk '$i' OK!"
	fi
fi
}

nVirtualCores=$(nproc --all)
nWorkers=$(( $nVirtualCores * 1 )) # I want 1 process per core

worker() { # $1 = idWorker
  echo "Worker $1 GO!"
  idAsset=0
  for (( asset=0; asset<=nbchk; asset++ )); do 
    # split assets among workers (using modulo); each worker will go through
    # the list and select the asset only if it belongs to that worker
    (( idAsset % nWorkers == $1 )) && task $1 $asset
    (( idAsset++ ))
  done
  echo "    Worker $1 ALL DONE!"
}

for (( idWorker=0; idWorker<nWorkers; idWorker++ )); do
  # start workers in parallel, use 1 process for each
  worker $idWorker &
done
wait # until all workers are done
