# pHASE-Stitcher
a python program to stitch the ReadBack phased haplotypes in F1 hybrids.

This program is designed to extend the phase states of paritally phased F1 hybrids using 1) allele frequency of the population the F1 hybrid came from 2) paritally phased haplotype blocks of the population resequence data.
This program exclusively uses the partial phase states generated using the tool "phaser" https://github.com/secastel/phaser
and is designed to work exclusively for F1 hybrid for now.

# Note: This is written in python3 and requires the following modules
pandas
io
pyvcf
itertools
collections
functools


#Usage: 

python Stitcher_using_1stOrderMarkov_InteractiveMode.py --vcf1 MY_subSample.vcf --vcf2 SP_subSample.vcf --pop1 My --pop2 Sp --output test --het_vcf F1_subSample.vcf --f1_sample 2ms04h


