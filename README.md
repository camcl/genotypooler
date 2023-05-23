# genotypooler

## Description
This project implements a simulation of SNP genotype pooling with a simple shifted transversal design.
The block size chosen for the pooling design is 4*4, with 8 pools and a design weight of 2.
The encoding and decoding part of the pooling procedure can be represented as follows: ![Pooling simulation on genotypes in the DNA Sudoku style](pooling-sim-gtgl.png)

where {0, 2, -1} are the allelic dosages from the true genotypes values at one SNP of any sample in (a). {0, 2, -1} stand for 
homozygote reference allele, homozygote alternate allele, missing genotype. As we represent a simulation with the genotype data of inbred lines, we assume that the samples are full homozygotes and there is therefore no heterozygous genotype (genotype value 1).

## Set up
* a Python 3.6 environment with packages listed in [requirements.txt](requirements.txt), e.g. for a Linux-based OS from the **genotypooler** folder:

(if `venv` for Python 3.6 is not installed: `apt install libpython3.6-dev python3.6-venv`)

`/usr/bin/python3.6 -m venv venv3.6`

`source venv3.6/bin/activate`

`pip install --upgrade pip`

`pip install -r requirements.txt`

(see [official venv documentation](https://docs.python.org/3/library/venv.html))

* **bcftools** installed on the OS. See [official page](https://samtools.github.io/bcftools/).

* **tabix**

## Usage
Some data and scripts are provided as use cases in [/examples](/examples). 
In particular, the following files can be found:
* *adaptive_gls.csv*: posterior genotypes probabilities of pooled individuals
* *Chr1.Founders.vcf.gz* and its index *.csi*: genotypes of the founder individuals in unphased GT format
* *Chr1.SNPs.pruned.nomiss.vcf.gz* and its index *.csi*: genotypes of the study samples in unphased GT format
* *Chr1.SNPs.pruned.nomiss.gl.vcf.gz* and its index *.csi*: genotypes of the study samples in GL format
* *pooling-ex.py*: a minimalistic command-line program for simulating SNPs genotypes pooling from VCF files
* *pooling-imputing-ex.ipynb*: a pipeline showing pooling simulation, imputation in pooled data with Beagle and impuatation quality visualization.
* directory *prophaser*: example of file imputed with prophaser
* directory *beagle-with-map*: example of file imputed with Beagle 4.1
* *argsfile_example.txt*: file with arguments required for plotting the concordance and the cross-entropy in the imputed study populations

## References
* [DNA Sudoku pooling designs](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2704425/pdf/1243.pdf/?tool=EBI)
* Beagle 4.1 articles for [phasing](https://linkinghub.elsevier.com/retrieve/pii/S0002929707638828) and [imputation](https://www.cell.com/ajhg/fulltext/S0002-9297(15)00491-7) 
* Beagle 4.1 [documentation and binaries](https://faculty.washington.edu/browning/beagle/b4_1.html)
* The [1000 Genomes Project](https://www.internationalgenome.org/) and its [VCF phase 3 data release](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/).
* Our paper in BMC Bioinformatics: ["A joint use of pooling and imputation for genotyping SNPs"](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04974-7)