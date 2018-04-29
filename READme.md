# phASE-Stitcher

***A python program to segregate and stitch the ReadBackPhased genotypes in F1 hybrids using first order markov chain and transition probabilities.\
This tool can be used as a companion tool along with [`phase-Extender`](https://github.com/everestial/phase-Extender) or as a standalone tool.***

Developed by [Bishwa K. Giri](mailto:kirannbishwa01@gmail.com) in the [Remington Lab](https://biology.uncg.edu/people/david-remington/) at the University of North Carolina at Greensboro, Biology department.

## Citation
Giri, B. K., Remington D. L. Haplotype phase extension and preparation of diploid genome using phase-Extender and phase-Stitcher. biorxiv (2018) [not uploaded yet].

## AUTHOR/SUPPORT
Bishwa K. Giri (bkgiri@uncg.edu; kirannbishwa01@gmail.com) \
Support @ https://groups.google.com/d/forum/phase-extender

## Intro to ReadBackPhasing
**Check these links for details on readbackphasing**
- https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_phasing_ReadBackedPhasing.php
- https://github.com/secastel/phaser/tree/master/phaser

<br>
<br>

# BACKGROUND
Haplotype phasing is a second "go to" problem in bioinformatics after read alignment.
The importance of haplotype phasing applies directly to the analyses of ASE (allele specific expression),
preparation of extended haplotype for EHH (extended haplotype homozygosity) test,
and preparation of dipolid genome which will soon be a new standard in bioinformatics in coming years, etc.
The necessity for haplotype phasing (and eventually diploid genome) increases with the increase in heterozygosity
in the genome, because higher hetetogeneity leads to larger alignment bias and complicates the reliability of
the variants that are called using that alignment data (SAM, BAM files).

Gene expression quantification using sequence reads is occupying a dominant and standard sphere of functional
analyses of putative genetic loci.
Quantification of gene expression comes in two flavors:
  - DE (differential expression) > gene expression differences quantified between two individuals or groups categorized by population, treatments, space or time.
  - ASE (allele specific expression) > gene expression differences quantified between the two alleles of the same locus within a diploid individual, categorized mainly by haplotypes but may be further grouped by population, treatments, space or time.

Quantification of RNAseq reads using reference genome based approach mainly involves haploid genome.
ASE (allele specific expression) quantification using alignment of RNAseq reads from F1 hybrids or outbred individuals
on haploid reference genomes is more susceptible to biased ASE observation considering following factors:
  - alignment to reference genome will likely trigger higher mapping of the reads from the population closer to the reference genome.
  - using allele masking based approach on haploid reference to address alignment bias is again more likely to attract
  paralogous alignment creating futher biases.
  
Several new approaches have been therefore created for better estimation of ASE analyses. Most optimal approach to ASE however
involves **1)** preparation of phased genome of the hybrids or outbred individual, **2)** then preparation of diploid genome and/or
transcriptome and then **3)** competitive alignment of the reads on diploid genome using conserved strategy.
The first step involving haplotype phasing has mostly concentrated around fixing the phase state in humans who have highly
homogenous genome, in inbred lines of mice and in other model systems that have lots of haplotype reference panels available.
Existing phase methods of F1 involves `mom, dad, child` trio, which is not optimal when parental information
are missing, or in the natural hybrids where parental identification is not possible. Also, existing trio methods are mainly genotype
based which take allele at single position at once.

ASE (allele specific expression) analyses which aims to identify cis regulatory factors and mechanism underlying gene expression
differences are heading toward more genomic and rnaseq based approaches.
The full resolution of ASE therefore relies on the quality of the phased dipoid genome.


**`phase-Stitcher`** is designed to utilize RBphase data with population based haplotype phasing to solve phase states in F1s.
The approach is to take RBphased haplotype blocks of a single F1 hybrid individual and several haplotype from the two different
parental background of the F1, then segregate the haplotype of F1 by computing likelihood against two parental background.
The advantages of using `phase-Stitcher` is explained below:
  - With increase in the size of sequence reads we are able to generate larger RBphased fragments.
  These fragments are again considerably larger when a heterogenous population is sequenced.
  F1 hyrbids of these heterogenous population have even larger RBphased fragments.
  Thus, haplotype phasing using RBphased data with population based likelihood estimates provides more optimal approach
  to solving phase state.
  - This tool doesn't require exact `maternal, parental` genotype data to solve phase state in F1.
  Rather phasing can be casually approached by supplying genotype data from `maternal vs. parental` background.


## ![#f03c15](https://placehold.it/15/f03c15/000000?text=+)Data Requirements

**phASE-Stitcher** can be used with the multi-sample vcf files produced by GATK pipeline or other tools that generate
readbackphased haplotype blocks in the output VCF.
A haplotype file is created using the RBphased VCF and then piped into **phase-Stitcher**.
See, this example for data structure of input haplotype file
[sample input haplotype file01](https://github.com/everestial/pHASE-Stitcher/blob/master/example_data/test_file01.txt)
- a tab separated text file with `PI` and `PG_al` value for each samples.



## Algorithm
For the **mcve** regarding the algorithm see this issue on [**stackoverflow**]() and/or [**my blog**](). 

<br>
<br>

# Tutorial

## Setup
**phase-Stitcher** is written in python3 interpreter. So, it can be directly run using the **".py"** file,
given all the required modules (dependencies) are installed.

**Runs on `Python 3.x` and has the following dependencies:**

  - [argparse](https://docs.python.org/3/library/argparse.html)
  - [collections](https://docs.python.org/3/library/collections.html?highlight=collections#module-collections)
  - [itertools](https://docs.python.org/3/library/itertools.html?highlight=itertools)
  - [io](https://docs.python.org/3/library/io.html?highlight=io#module-io)
  - [math](https://docs.python.org/3/library/math.html)
  - [multiprocessing](https://docs.python.org/3/library/multiprocessing.html?highlight=multiprocessing#)
  - [pandas](http://pandas.pydata.org/)
  - [re](https://docs.python.org/3/library/re.html?highlight=re#module-re)
  - [resource](https://docs.python.org/3/library/resource.html?highlight=resource#module-resource)
  - [time](https://docs.python.org/3/library/time.html?highlight=time#module-time)
  - [shutil](https://docs.python.org/3/library/shutil.html?highlight=shutil#module-shutil)
  
<br>
  
## Installation
```
pip3 install -r requirements.txt

#Or, if sudo previlages required, use:
sudo pip3 install -r requirements.txt
```

**If there is issues with installation while using `requirements.txt` install each dependencies individually.**\
e.g: `sudo python3 -m pip install pandas`  
  
<br>

## Usage
  Requires a readbackphased `haplotype file` as input and returns segregated and stitched haplotype file in both wide and long format.


Check this detailed [step by step tutorial](https://github.com/everestial/pHASE-Stitcher/wiki) for preparation
of `input files` and know-how about running `phase-Stitcher`.
    
<br>

## Input data

***haplotype file (required):*** Input haplotype file. Should contain `PI` and `PG_al` values for each sample.
To convert the vcf file to haplotype file (from VCF to proper text format) use **Step 01 (a)** in
the [tutorial](https://github.com/everestial/pHASE-Stitcher/wiki).
The sample name should not contain "`_`" character.


## Output data
contd.... 





