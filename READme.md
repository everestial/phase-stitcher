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
Haplotype phasing is a second "go to" problem in bioinformatics after read alignment. The importance of haplotype phasing applies directly to the analyses of ASE (allele specific expression), preparation of extended haplotype for EHH (extended haplotype homozygosity) test, and preparation of dipolid genome which will soon be a new standard in bioinformatics in coming years, etc. The necessity for haplotype phasing (and eventually diploid genome) increases with the increase in heterozygosity in the genome, because higher hetetogeneity leads to larger alignment bias and complicates the reliability of the variants that are called using that alignment data (SAM, BAM files).

ASE (allele specific expression) analyses using haploid reference genomes in hybrid individuals or population are more susceptible to biased ASE observation considering following factors:
  - alignment to reference genome will likely trigger higher mapping of the reads from the population closer to the reference genome.
  - using allele masking to address alignment bias is again more likely to attract paralogous alignment creating futher biases.
  
Optimal approach to ASE therefore involves preparation of phased genome of the hybrids, then preparation of diploid genome and then competitive alignment of the reads using more conserved strategy. 

Haplotype phasing has mostly concerntrated around fixing the phase state in humans who have highly homogenous genome and/or in inbred lines of mice and largely studied model systems, which have lots of haplotype reference panels available.  Existing methods in haplotype phasing therefore have not addressed issues with genome that are highly heterogenous or are hybrids. ASE (allele specific expression) which aims to identify cis regulatory factors and mechanism underlying gene expression differeces at the same loci are heading toward more genomic and rnaseq approaches. The full resolution of ASE therefore relies on the quality of the phased dipoid genome, which requires genome wide phasing. The only possible method to phase F1 hybrids is by taking `mom, dad, child` trio, which is not optimal when parental information are missing, or in the natural hybrids where parental identification is not possible. Also, existing trio methods are genotype based which take allele at single position at once. 

**`phase-Stitcher`** is designed to overcome this limitations in haplotype phasing by utilizing RBphased haplotype blocks. The advantages of using `phase-Stiticher` is explained below:
  - With increase in the size of sequence reads we are able to generate larger RBphased fragments. These fragments are again considerably larger when a heterogenous population is sequenced. F1 hyrbids of these heterogenous population have even larger RBphased fragments. Thus, haplotype phasing using RBphase data is most optimal approach to solving haplotype phasing.
  - This tool doesn't require exact `maternal, parental` genotype data to solve phase state in F1. Rather phasing can be casually approached by supplying genotype data from `maternal vs. parental` background. 
  - This tools is also able to segregate haplotype into maternal vs. paternal background, thus solving haplotype phase state genome wide for each parents. 


## ![#f03c15](https://placehold.it/15/f03c15/000000?text=+)Data Requirements

**phASE-Stitcher** can be used with the multi-sample vcf files produced by GATK pipeline or other tools that generate readbackphased haplotype blocks in the output VCF. A haplotype file is created using the RBphased VCF and then piped into **phase-Stitcher**. See, this example for data structure of input haplotype file [sample input haplotype file01]() - a tab separated text file with `PI` and `PG_al` value for each samples.



## Algorithm
For the **mcve** regarding the algorithm see this issue on [**stackoverflow**]() and/or [**my blog**](). 

<br>
<br>

# Tutorial

## Setup
**phase-Extender** is written in python3 interpreter. So, it can be directly run using the **".py"** file, given all the required modules (dependencies) are installed.

**Runs on `Python 3.x` and has the following dependencies:**

  - [argparse](https://docs.python.org/3/library/argparse.html)
  - [collections](https://docs.python.org/3/library/collections.html?highlight=collections#module-collections)
  - [itertools](https://docs.python.org/3/library/itertools.html?highlight=itertools)
  - [io](https://docs.python.org/3/library/io.html?highlight=io#module-io)
  - math
  - [multiprocessing](https://docs.python.org/3/library/multiprocessing.html?highlight=multiprocessing#)
  - [pandas](http://pandas.pydata.org/)
  - re
  - [resource](https://docs.python.org/3/library/resource.html?highlight=resource#module-resource)
  - [time](https://docs.python.org/3/library/time.html?highlight=time#module-time)
  - shutil
  
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
  Requires a readbackphased `haplotype file` as input and returns segregated and stitched haplotype file and other results files containing statistics on the initial vs. extended haplotype. 


Check this detailed [step by step tutorial](https://github.com/everestial/pHASE-Stitcher/wiki) for preparation of `input files` and know-how about running `phase-Stitcher`.
    
<br>

## Input data

***haplotype file (required):*** Input haplotype file. Should contain `PI` and `PG_al` values for each sample.\
To convert the vcf file to haplotype file (from VCF to proper text format) use **Step 01 (a)** in the tutorial. \
The sample name should not contain "`_`" character.


## Output data
contd.... 





