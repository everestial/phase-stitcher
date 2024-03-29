# PhaseStitcher

***A python program to segregate and stitch the ReadBackPhased genotypes in F1 hybrids to prepare
a genome wide haplotype using first order markov chain and transition probabilities.\
This tool can be used as a companion tool along with
[`phase-Extender`](https://github.com/everestial/phase-Extender) or as a standalone tool.***

Developed by [Bishwa K. Giri](mailto:kirannbishwa01@gmail.com) in
the [Remington Lab](https://biology.uncg.edu/people/david-remington/) at the
University of North Carolina at Greensboro, Biology department.

- [PhaseStitcher](#phasestitcher)
  - [Citation](#citation)
  - [AUTHOR/SUPPORT](#authorsupport)
  - [Intro to ReadBackPhasing](#intro-to-readbackphasing)
  - [BACKGROUND](#background)
  - [Data Requirements](#data-requirements)
  - [Algorithm](#algorithm)
  - [Tutorial](#tutorial)
    - [Prerequisites](#prerequisites)
    - [Installation from pypi:](#installation-from-pypi)
    - [Installation  and setup from source (Optional)](#installation--and-setup-from-source-optional)
  - [Usage](#usage)
  - [Sample example](#sample-example)
  - [Output Files](#output-files)
    - [*f1Sample*_haplotype_long.txt](#f1sample_haplotype_longtxt)
    - [*f1Sample*_haplotype_wide.txt](#f1sample_haplotype_widetxt)
    - [*f1Sample*_haplotype_stats.txt](#f1sample_haplotype_statstxt)
  - [Some Q/A on phase-stitcher](#some-qa-on-phase-stitcher)

## Citation

Giri, B. K., Remington D. L. Haplotype phase extension and preparation of
diploid genome using phase-Extender and phase-Stitcher. biorxiv (2018) [not uploaded yet].

## AUTHOR/SUPPORT

Bishwa K. Giri (bkgiri@uncg.edu; kirannbishwa01@gmail.com) \
Support @ <https://groups.google.com/d/forum/phase-extender>

## Intro to ReadBackPhasing

**Check these links for details on readbackphasing*

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

### Installation from pypi:

PhaseStitcher is hosted on pypi. So, you can install it using pip as:
```bash
$ pip install phase-stitcher

# After installation is complete you can run help function to get its parameters.

$ phase-stitcher -h
```
Now you can jump to usage and replace `python3 phase_stitcher.py` with `phase-stitcher`.

### Installation  and setup from source (Optional)

1. Clone this repo.

``` bash
git clone https://github.com/everestial/phASE-Stitcher
cd phASE-Stitcher
```

2. Make virtual env for python and install requirements.

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

3. To run tests locally:

  ``` bash
    pip install pytest
    pytest .
   ```

## Usage

```
$ phase-stitcher --help
usage: phase-stitcher [-h] [--nt NT] --input INPUT --pat PAT --mat MAT --f1Sample F1SAMPLE [--outPatMatID OUTPATMATID] [--output OUTPUT]
                      [--lods LODS] [--culLH CULLH] [--chr CHR] [--hapStats HAPSTATS]

options:
  -h, --help            show this help message and exit
  --nt NT               number of process to run -> The maximum number of processes that can be run at once is the number of different chromosomes
                        (contigs) in the input haplotype file.
  --input INPUT         name of the input haplotype file -> This haplotype file should contain unique index represented by 'PI' and phased genotype
                        represented by 'PG_al' for all the samples.
  --pat PAT             Paternal sample or comma separated sample names that belong to Paternal background. Sample group may also be assigned using
                        prefix. Options: 'paternal sample name', 'comma separated samples', 'pre:...'. Unique prefix (or comma separated prefixes)
                        should begin with 'pre:'.
  --mat MAT             Maternal sample or sample names (comma separated) that belong to maternal background. Sample group can also be assigned
                        using unique prefix/es. Options: 'maternal sample name', 'comma separated samples', 'pre:...'. Unique prefix (or comma
                        separated prefixes) should begin with 'pre:'.
  --f1Sample F1SAMPLE   Name of the F1-hybrid sample. Please type the name of only one F1 sample.
  --outPatMatID OUTPATMATID
                        Prefix of the 'Paternal (dad)' and 'Maternal (mom)'genotype in the output file. This should be a maximum of three letter
                        prefix separated by comma. Default: 'pat,mat'.
  --output OUTPUT       Name of the output directory. Default: f1SampleName + '_stitched'
  --lods LODS           log(2) odds cutoff threshold required to assign maternal Vs. paternal haplotype segregation and stitching.
  --culLH CULLH         Cumulative likelhood estimates -> The likelhoods for haplotype segregation can either be max-sum vs. max-product. Default:
                        maxPd i.e max-product. Options: 'maxPd' or 'maxSum'.
  --chr CHR             Restrict haplotype stitching to a specific chromosome.
  --hapStats HAPSTATS   Computes the descriptive statistics of final haplotype. Default: 'no'.Option: 'yes', 'no' .

```

>NOTE Input haplotype file should contain `PI` and `PG_al` values for each sample.

  Requires a readbackphased `haplotype file` as input and returns segregated and stitched haplotype file in both wide
  and long format. Descriptive statistics of the final haplotype can also be produced if desired.

Check this detailed [step by step tutorial](https://github.com/everestial/pHASE-Stitcher/wiki) for preparation
of `input files` and know-how about running `phase-Stitcher`.

## Sample example 
```
$ phase-stitcher --nt 1 --input tests/inputs/haplotype_file01.txt --mat MA605 --pat Sp21 --f1Sample ms02g --culLH maxSum --lods 3 --hapStats yes

  - using haplotype file "tests/inputs/haplotype_file01.txt" 
  - F1-hybrid of interest: "ms02g" 
  - using "1" processes 
  - using log2 odds cut off of "3" 
  - using "max sum" to estimate the cumulative maximum likelyhood while segregating the diploid haplotype block into maternal vs. paternal haplotype 
  - statistics of the haplotype before and after extension will be prepared for the sample of interest i.e "ms02g" 
#######################################################################
        Welcome to phase-Stitcher version 1.2       
  Author: kiran N' bishwa (bkgiri@uncg.edu, kirannbishwa01@gmail.com) 
#######################################################################


##########################
 - Worker maximum memory usage: 85008.00 (mb)

Completed haplotype segregation and stitching for all the chromosomes.
Time elapsed: '0.341528' sec. 
 - Global maximum memory usage: 87792.00 (mb)
Completed writing the dataframes .....
The End :)
```

For exploratory analysis of statistics generated from file: visit this [notebook](EDA%20on%20hapstats%20data%20generated%20from%20phase%20stitcher.ipynb).

## Output Files

### *f1Sample*_haplotype_long.txt

Final haplotype for **f1Sample** of interest after phase segregation in **long format**.

- **CHROM** - Contig name (or number).
- **POS** - Start position of haplotype (1 based).
- **REF** - Reference allele at that site.
- **all-alleles** - All the alleles represented by all the samples in the input file at that site.
- **_f1Sample_:PI** - Unique `PI` index of the haplotype blocks for sample of interest.
- **_f1Sample_:PG_al** - Phased GT (genotype) alleles at the genomic position that belong to unique `PI` indexes.
- **log2Odds** - log2 of Odds computed between the left vs. right haplotype against observed haplotype in
paternal vs. maternal samples.
- **_pat_ _hap** - Haplotype that belongs to paternal background based on **_lods_** cutoff.
- **_mat_ _hap** - Haplotype that belongs to maternal background based on **_lods_** cutoff.

### *f1Sample*_haplotype_wide.txt

- Final haplotype for **f1Sample** of interest after phase segregation in **wide format**.
- All the headers are the same as file in **_long format_** except **_POS_Range_**

### *f1Sample*_haplotype_stats.txt

Descriptive haplotype statistics of the input haplotype file for the sample of interest. These statistics
can be used to compute the distribution of several values (lods, number or variants etc.) between phased
and unphased haplotype blocks and if they were assigned to final genome wide haplotype.

- **CHROM** - Contig name (or number).
- **phasedBlock** - Blocks that were phased to genome wide haplotype based on **_lods_** cutoff.
- **unphasedBlock** - Blocks that were not phased to genome wide haplotype based on **_lods_** cutoff.
- **numVarsInPhasedBlock** - Number of variants in each **_phasedBlock_**.
- **numVarsInUnPhasedBlock** - Number of variants in each **_unphasedBlock_**.
- **log2oddsInPhasedBlock** - Calculated **log2Odds** in each **_phasedBlock_**.
- **log2oddsInUnPhasedBlock** - Calculated **log2Odds** in each **_unphasedBlock_**.
- **totalNumOfBlock** - Total number of RBphased blocks in the given **_f1Sample_**.
- **totalNumOfVars** - Total number of variants in the given **_f1Sample_**.

**Note:** - The **block index i.e PI** in **_phasedBlock_** and in **_unphasedBlock_**,
and it's associated statistics are in order.

## Some Q/A on phase-stitcher

The conjoined **Q/A** for **_phase stitcher_** is covered under **Q/A** for
