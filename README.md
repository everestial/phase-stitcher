# pHASE-Stitcher
a python program to stitch the ReadBack phased haplotypes in F1 hybrids.

This program is designed to extend the phase states of paritally phased F1 hybrids using 1) allele frequency of the population the F1 hybrid came from 2) paritally phased haplotype blocks of the population resequence data.

Usage: python Stitcher_using_1stOrderMarkov_InteractiveMode.py --vcf1 MY_subSample.vcf --vcf2 SP_subSample.vcf --pop1 My --pop2 Sp --output test --het_vcf F1_subSample.vcf --f1_sample 2ms04h


