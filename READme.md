# phASE-Stitcher

***A python program to segregate and stitch the ReadBackPhased genotypes in F1 hybrids to prepare
a genome wide haplotype using first order markov chain and transition probabilities.\
This tool can be used as a companion tool along with
[`phase-Extender`](https://github.com/everestial/phase-Extender) or as a standalone tool.***

Developed by [Bishwa K. Giri](mailto:kirannbishwa01@gmail.com) in 
the [Remington Lab](https://biology.uncg.edu/people/david-remington/) at the 
University of North Carolina at Greensboro, Biology department.

## Citation

Giri, B. K., Remington D. L. Haplotype phase extension and preparation of 
diploid genome using phase-Extender and phase-Stitcher. biorxiv (2018) [not uploaded yet].

## AUTHOR/SUPPORT

Bishwa K. Giri (bkgiri@uncg.edu; kirannbishwa01@gmail.com) \
Support @ https://groups.google.com/d/forum/phase-extender

## Intro to ReadBackPhasing

**Check these links for details on readbackphasing**

- <https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_phasing_ReadBackedPhasing.php>
- <https://github.com/secastel/phaser/tree/master/phaser>

## BACKGROUND

Haplotype phasing is a second "go to" problem in bioinformatics after read alignment.
The importance of haplotype phasing applies directly to the analyses of ASE (allele specific expression),
preparation of extended haplotype for EHH (extended haplotype homozygosity) test,
and preparation of dipolid genome which will soon be a new standard in bioinformatics in coming years, etc.
The necessity for haplotype phasing (and eventually diploid genome) increases with the increase in heterozygosity
in the genome, because higher hetetogeneity leads to larger alignment bias and complicates the reliability of
the variants that are called using that alignment data (SAM, BAM files).

Gene expression quantification using sequence reads is occupying a dominant and standard sphere of functional
analyses of putative genetic loci.
Quantification of gene expression as a test of phenotypic differences comes in two flavor:

  - DE (differential expression) > gene expression differences quantified between two individuals or 
  groups categorized by population, treatments, space or time.
  - ASE (allele specific expression) > gene expression differences quantified between the two alleles of 
  the same locus in a diploid individual, categorized mainly by haplotypes but may be further 
  grouped by population, treatments, space or time.

Quantification of RNAseq reads using reference genome based approach mainly involves haploid genome.
ASE (allele specific expression) quantification using alignment of RNAseq reads from F1 hybrids or outbred individuals
on haploid reference genomes is more susceptible to biased ASE observation considering following factors:
  - alignment to reference genome will likely trigger higher mapping of the reads from the 
  population closer to the reference genome.
  - using allele masking based approach on haploid reference to address alignment bias is again more likely to attract
  paralogous alignment in the masked region creating futher biases.
  
Several new approaches have been created for better estimation of ASE analyses. Most optimal approach to 
ASE however involves **1)** preparation of phased genome of the hybrids or outbred individual, 
**2)** then preparation of diploid genome and/or transcriptome and then **3)** competitive alignment of the reads 
on diploid genome using conserved strategy.
The first step involving haplotype phasing has mostly concentrated around fixing the phase state in 
humans who have highly homogenous genome, in inbred lines of mice and in other model systems that have lots of 
haplotype reference panels available.
Existing phase methods of F1 hybrids involves `mom, dad, child` trio, which is not optimal when parental information
are missing, or in the natural hybrids where parental identification is not possible. 
Also, existing trio methods are mainly genotype based which take allele at single position at once.

ASE (allele specific expression) analyses which aims to identify cis regulatory factors and mechanism underlying gene expression
differences are heading toward more genomic and rnaseq based approaches.
The full resolution of ASE therefore relies on the quality of the phased dipoid genome.


**`phase-Stitcher`** is designed to utilize RBphase data with population based haplotype phasing 
to solve phase states in F1s. The approach is to take RBphased haplotype blocks of a single F1 hybrid individual 
and several haplotype from the two different parental background of the F1, then segregate the haplotype of F1 
by computing likelihood of each haplotype against two parental background. **The advantages of using `phase-Stitcher` is 
explained below:**
  - With increase in the size of sequence reads (mainly PE i.e paired end reads) we are able to generate larger 
  RBphased fragments. These fragments are again considerably larger when a heterogenous population is sequenced.
  F1 hyrbids of these heterogenous population have even larger RBphased fragments.
  Thus, haplotype phasing using RBphased data with population based likelihood estimates provides more optimal approach
  to solving phase state.
  - This tool doesn't require exact `maternal, parental` genotype data to solve phase state in F1.
  Rather phasing can be casually approached by supplying genotype data from `maternal vs. parental` background.


## Data Requirements

**phASE-Stitcher** can be used with the multi-sample vcf files produced by GATK pipeline or other tools that generate
readbackphased haplotype blocks in the output VCF.
A HAPLOTYPE file is created using the RBphased VCF and then piped into **phase-Stitcher**. 
Use [VCF-Simplify](https://github.com/everestial/VCF-simplify) to prepare HAPLOTYPE file from multisample VCF.
See, this example for data structure of input haplotype file
[sample input haplotype file01](https://github.com/everestial/pHASE-Stitcher/blob/master/example_01/haplotype_file01.txt)
- a tab separated text file with `PI` and `PG_al` value for each samples.



## Algorithm

For the **mcve** regarding the algorithm see this issue on [**stackoverflow**]() and/or [**my blog**](). 

## Tutorial

### Prerequisites

**phASE-Stitcher** is written in python3, so you need to have python3 installed on your system to run this code locally. If you don't have python installed then, you can install from [here](https://www.python.org/downloads/). For linux; you can get latest python3 by:

`sudo apt-get install python3`

### Installation  and setup

1. Clone this repo.

``` bash
git clone https://github.com/everestial/phASE-Stitcher
cd phASE-Stitcher
```

1. Make virtual env for python and install requirements.

``` bash
python3 -m venv .env
source .env/bin/activate   # for linux
.env\Scripts\activate      # for windows
pip install -r requirements.txt
```

OR, you can install latest versions individually by:

``` bash
pip install pandas numpy matplotlib

```


1. To run tests locally:

  ``` bash
    pip install pytest
    pytest .
   ```


## Usage

  Requires a readbackphased `haplotype file` as input and returns segregated and stitched haplotype file in both wide 
  and long format. Descriptive statistics of the final haplotype can also be produced if desired.


Check this detailed [step by step tutorial](https://github.com/everestial/pHASE-Stitcher/wiki) for preparation
of `input files` and know-how about running `phase-Stitcher`.
    
<br>

## Input data

***haplotype file (required):*** Input `haplotype` file. Should contain `PI` and `PG_al` values for each sample.

### Performance Related

* **--nt** _(1)_ - maximum number of processes to run at once. 
The maximum number of processes is limited to number of chromosomes (contigs) in the input haplotype file. 

## Arguments
### Required
* **--f1Sample** - name of the f1 hybrid sample of interest. 
It should refer to a single sample in the haplotype the file. 
* **--pat** - Paternal sample or comma separated sample names that
                        belong to "paternal" background. Sample group may also
                        be assigned using prefix. Options: 'paternal sample
                        name', 'comma separated samples', 'pre:...'. Unique
                        prefix (or comma separated prefixes) should begin with
                        'pre:'.
* **--mat** - Maternal sample or sample names (comma separated) that
                        belong to "maternal" background. Sample group can also
                        be assigned using unique prefix/es. Options: 'maternal
                        sample name', 'comma separated samples', 'pre:...'.
                        Unique prefix (or comma separated prefixes) should
                        begin with 'pre:'. 

### Optional
* **--python_string** _(python3)_ - Calls `python 3` interpreter to run the program.
* **--output** _(f1Sample_extended)_ - Name of the output directory.
* **----outPatMatID**_(pat,mat)_- Prefix of the **paternal (dad)** and **maternal (mom)** genotype in the output file.
                             This should be a maximum of three letter prefix separated by comma.
                             
* **--lods** _(5)_ - log(2) odds cutoff threshold required to assign maternal Vs. paternal 
                               haplotype segregation and stitching.  
               &emsp; &emsp; - Positive calculated log2Odd above lods cut off will assign the left haplotype to paternal allele.\
               &emsp; &emsp; - Negative calculated log2Odd above lods cut off will assign the left haplotype to maternal allele.\
               &emsp; &emsp; - Calculated abs |log2Odd| below lods cut off threshold will leave hapltoype unassigned.

* **--culLH** _(maxPd)_ - Method for Cumulative likelhood estimates -> The likelhoods for haplotype segregation can 
                            either be **max-sum** vs. **max-product**. 
                            ***Default*** is "max-product". ***Options:*** 'maxPd' or 'maxSum'.
* **--chr** - Restrict haplotype stitching to a specific chromosome.
* **--hapStats** _(no)_ - Computes the descriptive statistics of final haplotype. **Options:** 'yes', 'no'

<br>
<br>

## Output Files

### *f1Sample*_haplotype_long.txt
Final haplotype for **f1Sample** of interest after phase segregation in **long format**.
* 1 - **CHROM** - Contig name (or number).
* 2 - **POS** - Start position of haplotype (1 based).
* 3 - **REF** - Reference allele at that site.
* 4 - **all-alleles** - All the alleles represented by all the samples in the input file at that site.
* 5 - **_f1Sample_:PI** - Unique `PI` index of the haplotype blocks for sample of interest.
* 6 - **_f1Sample_:PG_al** - Phased GT (genotype) alleles at the genomic position that belong to unique `PI` indexes.
* 7 - **log2Odds** - log2 of Odds computed between the left vs. right haplotype against observed haplotype in 
paternal vs. maternal samples.
* 8 - **_pat_ _hap** - Haplotype that belongs to paternal background based on **_lods_** cutoff.
* 9 - **_mat_ _hap** - Haplotype that belongs to maternal background based on **_lods_** cutoff.

<br>

### *f1Sample*_haplotype_wide.txt
  - Final haplotype for **f1Sample** of interest after phase segregation in **wide format**.
  - All the headers are the same as file in **_long format_** except **_POS_Range_**

<br>

### *f1Sample*_haplotype_stats.txt
Descriptive haplotype statistics of the input haplotype file for the sample of interest. These statistics 
can be used to compute the distribution of several values (lods, number or variants etc.) between phased 
and unphased haplotype blocks and if they were assigned to final genome wide haplotype.

* 1 - **CHROM** - Contig name (or number).
* 2 - **phasedBlock** - Blocks that were phased to genome wide haplotype based on **_lods_** cutoff.
* 3 - **unphasedBlock** - Blocks that were not phased to genome wide haplotype based on **_lods_** cutoff.
* 4 - **numVarsInPhasedBlock** - Number of variants in each **_phasedBlock_**.
* 5 - **numVarsInUnPhasedBlock** - Number of variants in each **_unphasedBlock_**.
* 6 - **log2oddsInPhasedBlock** - Calculated **log2Odds** in each **_phasedBlock_**.
* 7 - **log2oddsInUnPhasedBlock** - Calculated **log2Odds** in each **_unphasedBlock_**.
* 8 - **totalNumOfBlock** - Total number of RBphased blocks in the given **_f1Sample_**.
* 9 - **totalNumOfVars** - Total number of variants in the given **_f1Sample_**.

**Note:** - The **block index i.e PI** in **_phasedBlock_** and in **_unphasedBlock_**, 
and it's associated statistics are in order.

## Some Q/A on phase-stitcher

The conjoined **Q/A** for **_phase stitcher_** is covered under **Q/A** for
[phase-extender](https://github.com/everestial/phase-Extender#some-qa-on-phase-extender)
