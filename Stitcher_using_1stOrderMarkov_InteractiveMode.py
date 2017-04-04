
#!/home/everestial007/anaconda3/bin/python
#/usr/local/bin/python3

# Note: this program uses python 3
# preparing a interactive version of the pHASE-Stitcher (using 1st order Markov chain)

import sys, os, time
from io import StringIO

try: import argparse
except: sys.exit('ArgumentParser module not found.\nPlease install it before.')

# for mining SNPs and phase states from vcf files
try: import vcf
except: sys.exit('pyvcf module not found.\nPlease install it before.')

from functools import  reduce

# for running markov chain process and testing of odds ratio
try: import pandas as pd
except: sys.exit('pandas module not found.\nPlease install it before.')
import itertools as it
from itertools import product
from collections import defaultdict


### To Do - Problem Stage #1:  !!! needs change
# 1) Add different arguments (pop1 and pop2 - these are name of the populations and will be added as header in the output file)
# 2) Add arguments (vcf1 and vc2) - these are vcf files that you should proide
# 3) Add argument for the name of output file

### !!! needs change
## 4a) Create an output file which should contain following headers
# contig	pos	id	ref	alt	ref-freq-My	ref-freq-Sp	alt-A-freq-My	alt-A-freq-Sp
# the headers are separted by tab
# 4b) write the following values from vcf files: contig (from CHR), pos, id, ref - without any changes
# When merging information from two vcfs the the intersecting values can be kept and non intersecting values can be added
# 4c) Read the AF values in the vcf1 and vcf2 files and write the corresponding values



#Step01: first define the main arguments.
# Description - parser.add_argument adds the argument to the program.
# This step is added just to make the program interactive so user can supply their data

def main():
    # Arguments that are created as input
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf1",
                        help="sorted vcf file from population 1", required=True)
    parser.add_argument("--vcf2",
                        help="sorted vcf file from population 2", required=True)
    parser.add_argument("--pop1", help="Name or prefix for population 1."
                                       " Explanation: Please type the name of the population this vcf was derived from.", required=True)
    parser.add_argument("--pop2", help="Name or prefix for population 2."
                                       " Explanation: Please type the name of the population this vcf was derived from.", required=True)
    parser.add_argument("--output", help="Prefix for the name of the output file", default="output", required=True)

    # agruments for vcf file from heterozygous individual
    parser.add_argument("--het_vcf",
                        help="sorted vcf file from heterozygote individual or f1 hybrid",
                        required=True)
    parser.add_argument("--f1_sample", help="the name of the sample in the vcf file for heterozygote individuals or f1"
                                            " Explanation: Please type one F1 sample name from vcf of F1 hybrids.", required=True)

    # debug / development / reporting - this function not complete yet
    parser.add_argument("--chr", default="", help="Restrict haplotype stitching to a specific chromosome.")

    ## other arguments to add - to do !!!!
    # number of minimum SNPs in a haplotype
    # cutoff for the odds ratios
    # phase Indels or not

    global args; #creating a global argument variable
    args = parser.parse_args()

    # setup
    version = "0.6";  ## The below message is printed to the screen when the program is run with provided files/variables.
    fun_flush_print("");
    fun_flush_print("##################################################")
    fun_flush_print(" Welcome to pHASE_Stitcher Markov Order 01");
    fun_flush_print("  Author: kiranNbishwa (bkgiri@uncg.edu) ")
    fun_flush_print("##################################################");
    fun_flush_print("");

    global devnull;
    devnull = open(os.devnull, 'w')  ##starts devnull and provides a writing permission

    fun_flush_print("#Step 01:  Checking the required files......\n");

    #Step 01-B checks if the required files (vcfs) are provided or not
    check_files = [args.vcf1, args.vcf2, args.het_vcf];
    for xfile in check_files:
        if xfile != "":
            if os.path.isfile(xfile) == False:
                fatal_error("File: %s not found." % (xfile));
                # reports error message if the checked files isn't found

    ## To Do -  check if the names of the population is given
    check_files = [args.pop1, args.pop2];
    for xfile in check_files:
        if xfile == "":
            fatal_error("File: %s name should be provided for each population." % (xfile));

    ## To Do - Do this for each stages..
    #print(time)

    #Step 01-C Now, load the VCF file cotaining allele frequency for populations
    if args.vcf1 != "":
        if os.path.isfile(args.vcf1) == True:
            #vcf1_file = vcf.Reader(filename=args.vcf1);  #here we created a new object called vcf1_data
            fun_flush_print("#1: Checking if vcf file for %s population exists ......" % (args.pop1));
        else:
            fatal_error("Allele frequency VCF (--vcf1) specified for population %s does not exist.\n" % (args.pop1));



    if args.vcf2 != "":
        if os.path.isfile(args.vcf2) == True:
            #vcf2_file = vcf.Reader(filename=args.vcf2);  #here we created a new object called vcf2_data
            fun_flush_print("#1: Checking if vcf file for %s population exists ......" % (args.pop2));
        else:
            fatal_error("Allele frequency VCF (--vcf2) specified for population %s does not exist.\n" % (args.pop2));



    if args.het_vcf != "":
        if os.path.isfile(args.het_vcf) == True:
            #het_vcf_file = vcf.Reader(filename=args.het_vcf);  #here we created a new object called het_vcf
            fun_flush_print("#1:  Checking the presence of vcf file %s exists for heterozygous sample %s ......\n" % (args.f1_sample, args.het_vcf));
        else:
            fatal_error("Phased VCF (--het_vcf) %s for f1 hybrid %s does not exist.\n" % (args.het_vcf, args.f1_sample));


    ## Defining the parameters for new functions   ## !!!! - needs verification
    ## These functions will be used for streaming the data from one function to another


    ## Now, read the following functions  - !!! find a way or making this more clean
    #stream_pop1_allele_table, stream_pop2_allele_table, stream_f1_allele_table = read_vcf()
    #merge_tables(stream_pop1_allele_table, stream_pop2_allele_table, stream_f1_allele_table)
    create_table_headers()
    #read_vcf(x,y,z)


## Step 02 - Now, we read the vcf files.
## Read the vcf files from two different populations and F1 hybrid
# And, extract data out of it.
# it can be done by creating a new independent function or by using the function builtin within the main function (Step01)

def create_table_headers():
    global args;

    fun_flush_print("#Step 02:  Reading the vcf file to prepare vcf-allele table......\n");
    fun_flush_print("#Step 02:  Creating the header lines for each vcf-allele table......\n");

    # Step 02-A **start**
    ## Before reading a vcf file, we create a file to write the output of the vcf
    ## Step 02...A: Create a table to write the information from the vcf file of population #01
    output1 = open(args.output + '_' +  args.pop1 + "_allele_table.txt", "w"); ## may not be necessary

    ## contd..... to do: Read the vcf file and prepare the output table at the same time - do this !!!

    pop1_data = open(args.vcf1, 'r')
    for lines in pop1_data.read().split('\n'):
        if '##' in lines:
            continue
        if '#CHROM' in lines:
            header = lines.split('\t')
            front_header = ('contig\tpos\tref\t' + 'all_alleles_' + args.pop1 +'\tall_freq_' + args.pop1 +'\t')
            sample_genotype = header[9:]  # starting from 10th column position of vcf file
            rear_header = list()

            for name in sample_genotype:
                rear_header.append(name + "_GT")
                rear_header.append(name + "_allele")
                rear_header.append(name + "_PG")
                rear_header.append(name + "_PG_allele")
                rear_header.append(name + "_PI")

            ## common way of writing output # deactivated though
            #output = open(args.output + '_' +  args.pop1 + "_allele_table.txt", "w")
            #output.write(front_header + "\t".join(rear_header) + '\n')
            #output.close()

            # this method write the data to file and also stores it in the terminal/memory/console
            header_pop1_allele_table = (front_header + "\t".join(rear_header))

            with open(args.output + '_' +  args.pop1 + "_allele_table.txt", "w") as f:
                print(header_pop1_allele_table, file=f)
                f.close()

            break

    fun_flush_print("#2-A: Creating a table to write data from the file %s ......" %(args.vcf1));


    ## Data from population #02 (file 02):    Read the data from vcf #02
    ## Create header for vcf 02

    pop2_data = open(args.vcf2, 'r')
    for lines in pop2_data.read().split('\n'):
        if '##' in lines:
            continue
        if '#CHROM' in lines:
            header = lines.split('\t')
            #front_header = "contig\tpos\tref\tall_alleles_Sp\tall_freq_Sp\t"
            front_header = ('contig\tpos\tref\t' + 'all_alleles_' + args.pop2 + '\tall_freq_' + args.pop2 + '\t')
            sample_genotype = header[9:]
            rear_header = list()

            for name in sample_genotype:
                rear_header.append(name + "_GT")
                rear_header.append(name + "_allele")
                rear_header.append(name + "_PG")
                rear_header.append(name + "_PG_allele")
                rear_header.append(name + "_PI")

            # this method write the data to file and also stores it in the terminal/memory/console
            header_pop2_allele_table = (front_header + "\t".join(rear_header))

            with open(args.output + '_' + args.pop2 + "_allele_table.txt", "w") as f:
                print(header_pop2_allele_table, file=f)
                f.close()

            break

    fun_flush_print("#2-A: Creating a table to write data from the file %s ......" % (args.vcf2));


    ## Step 02...A...contd....: Create a table to write the information from the vcf file of F1 hybrid
    ## Data from F1 hybrid vcf:    Read the data from vcf #03
    ## Create header for vcf 03

    f1_data = open(args.het_vcf, 'r')
    for lines in f1_data.read().split('\n'):
        if '##' in lines:
            continue
        if '#CHROM' in lines:
            header = lines.split('\t')
            front_header = ('contig\tpos\tref\tall_alleles_hybrid\tall_freq_hybrid\t')
            sample_genotype = header[9:]
            rear_header = list()

            for name in sample_genotype:
                rear_header.append(name + "_GT")
                rear_header.append(name + "_allele")
                rear_header.append(name + "_PG")
                rear_header.append(name + "_PG_allele")
                rear_header.append(name + "_PI")

            # also create a variable to store the data from f1 hybrid **
            header_f1_allele_table = front_header + "\t".join(rear_header)

            with open(args.output + 'f1_allele_table.txt', "w") as f:
                print(header_f1_allele_table, file=f)
                f.close()

            break

    fun_flush_print("#2-A: Creating a table to write data from the file %s ......\n" % (args.het_vcf));

    # clear memory
    del pop1_data, pop2_data, f1_data

    return read_vcf(header_pop1_allele_table, header_pop2_allele_table, header_f1_allele_table)


    ## ** Step 02-A end ** ## ************


    ## Step 02-B **start**
    ## Now, read the data from the vcf files to mine the data
def read_vcf(header_pop1_allele_table, header_pop2_allele_table, header_f1_allele_table):

    stream_pop1_allele_table = header_pop1_allele_table
    stream_pop2_allele_table = header_pop2_allele_table
    stream_f1_allele_table = header_f1_allele_table

    del header_pop1_allele_table, header_pop2_allele_table, header_f1_allele_table



    fun_flush_print("#2-B: Mining the GT, allele, PG (phased genotype), PG_allele and PI (phase Index) "
                    "values for each sample from each vcf ......\n")
    ## Reading vcf-01
    fun_flush_print("#2-B: Mining data from %s ......" % (args.vcf1))

    stream_pop1_vcf = vcf.Reader(open(args.vcf1, 'r'))
    sample_name = stream_pop1_vcf.samples
    sample_size = len(sample_name)

    for record in stream_pop1_vcf:
        contig1 = record.CHROM
        pos1 = record.POS
        ref_allele1 = record.REF

        # since alt alleles are represented as a list
        # alt_alleles1 = str(record.ALT[::]) would report this list as csv values and brackets
        # use the below code to report list as strings and without brackets []
        alt_alleles1 = ",".join(map(str, (record.ALT[::])))

        alt_freq1 = ",".join(map(str, (record.INFO['AF'])))
        # this values is returned as a list, so need to be separated or converted using index values..

        alt_freq_sum1 = sum(record.INFO['AF'])

        ref_freq1 = round(1 - alt_freq_sum1, 3)  # calculates ref frequency and rounds to 3.
        # Note: we can also write,   ref_freq1 = format((1-alt_freq_sum1), '.3f')

        ## Creating a string of all alleles (which are separated by commas).
        # This actually preserves the list and index of all the alleles in the vcf
        # The reference allele is placed first (so its index is 0)
        # the rest of the alleles (alt allele) are joined together (with their conserved)
        all_alleles1 = ref_allele1 + "," + alt_alleles1

        ## create the frequency values for all_alleles.
        # The index of the allele and it's corresponding frequency are conserved
        all_freq1 = str(ref_freq1) + "," + str(alt_freq1)

        ## Creating a list of all_alleles from the string
        # The separating commas are stripped
        # this is useful while mining the correpsonding allele and its phased state...
        # ... in the for-loop below
        allele_list = all_alleles1.split(',')

        # Now, write some of the values from this line of the vcf file
        output = open(args.output + '_' +  args.pop1 + "_allele_table.txt", "a")
        output.write("{}\t{}\t{}\t{}\t{}"
                     .format(contig1, pos1, ref_allele1, all_alleles1, all_freq1))

        stream_pop1_allele_table += ("\n{}\t{}\t{}\t{}\t{}"
                     .format(contig1, pos1, ref_allele1, all_alleles1, all_freq1))
        ## Note: adding extra '\n' in the stream_pop1_allele_table data since the new lines are added ...
        # ... right after the last lines

        # the rest of the record (information) are called ...
        # ... using for-loop for each samples in that vcf line.

        # for records in record:
        # print('\n printing records\n')
        # print(record) # prints the main records: chr#, pos, Ref allele and alt allele
        # print(record.samples) # prints the remaining record for each sample in the vcf file

        # the rest of the record (information) are called using for-loop for samples in that vcf line.
        for sample in record.samples:

            # extracting the gt_bases requires the sample name from vcf file
            # it was difficult to extract sample name in for-loop
            # so using below method to extract sample name

            # print(sample) gives me:
            # Call(sample=MA605, CallData(GT=. /., AD = [0, 0], DP = 0, PG =./., PL = [0, 0, 0], PW =./.))
            # print(type(sample))
            # which is following type of data: <class 'vcf.model._Call'>
            # print((str(sample).split())[0].split("=")[1].strip(','))
            # gives me the sample name as: MA605
            #print('\n printing samples \n')
            #print(sample)

            ### Now, call the the sample named 'sample_id' for further genotype information

            # old method of calling a sample is ... (below)

            # extracting the gt_bases requires the sample name from vcf file
            # it was difficult to extract sample name in for-loop
            # so using below method to extract sample name
            # print(sample) gives me:
            # Call(sample=MA605, CallData(GT=. /., AD = [0, 0], DP = 0, PG =./., PL = [0, 0, 0], PW =./.))
            # print(type(sample))
            # which is following type of data: <class 'vcf.model._Call'>
            # print((str(sample).split())[0].split("=")[1].strip(','))
            # gives me the sample name as: MA605
            # sample_id = (str(sample).split())[0].split("=")[1].strip(',')
            # call_id = record.genotype(sample_id)
            # genotype_base = call_id.gt_bases

            ### or new method is simply done
            # this method is explained as 'Call' objects in pyVCF document
            sample_id = record.genotype(sample.sample)
            genotype_base = sample_id.gt_bases

            ## Now, extract the phased allele using the 'PG' tag, not the phased 'GT' tag
            if sample['PG'] != './.':  # which could be anyone of '0/0', 1/1, 1/2, 1|2, 2|2, etc. etc.
                phased_PG = list(sample['PG'])
                first_allele_index = phased_PG[0]
                separator = phased_PG[1]
                second_allele_index = phased_PG[2]

                # Therefore...
                PG_allele = allele_list[int(first_allele_index)] + \
                            separator + allele_list[int(second_allele_index)]

            else:
                PG_allele = sample['PG']

            ## Now, write the remaining values to the file
            try:
                output.write("\t{}\t{}\t{}\t{}\t{}".format(sample['GT'], genotype_base, sample['PG'], PG_allele,
                                                           sample['PI']).replace('./.', '.').replace('None', '.'))

                stream_pop1_allele_table += ("\t{}\t{}\t{}\t{}\t{}".format(sample['GT'], genotype_base, sample['PG'], PG_allele,
                                                           sample['PI']).replace('./.', '.').replace('None', '.'))

            ## If the record.samples doesn't contain a particular tag (like PI in this case)..
            # ... it will throw an AttributeError (attribute not present) which can be simply handle using None..
            # .. replace by period(.).
            except AttributeError:
                output.write("\t{}\t{}\t{}\t{}\t{}".format(sample['GT'], genotype_base, sample['PG'],
                                                           PG_allele, '.').replace('./.', '.').replace('None', '.'))
                stream_pop1_allele_table += ("\t{}\t{}\t{}\t{}\t{}".format(sample['GT'], genotype_base, sample['PG'],
                                                           PG_allele, '.').replace('./.', '.').replace('None', '.'))

        output.write('\n')
        output.close()


    #### Data from population #02 (file 02):    Read the data from vcf #02
    fun_flush_print("#2-B: Mining data from %s ......" % (args.vcf2))
    stream_pop2_vcf = vcf.Reader(open(args.vcf2, 'r'))
    sample_name = stream_pop2_vcf.samples
    sample_size = len(sample_name)

    for record in stream_pop2_vcf:
        contig2 = record.CHROM
        pos2 = record.POS
        ref_allele2 = record.REF
        alt_alleles2 = ",".join(map(str, (record.ALT[::])))

        alt_freq2 = ",".join(map(str, (record.INFO['AF'])))
        alt_freq_sum2 = sum(record.INFO['AF'])

        ref_freq2 = round(1 - alt_freq_sum2, 3)
        all_alleles2 = ref_allele2 + "," + alt_alleles2
        all_freq2 = str(ref_freq2) + "," + str(alt_freq2)

        allele_list = all_alleles2.split(',')

        # Now, write some of the values from this line of the vcf file
        output = open(args.output + '_' + args.pop2 + "_allele_table.txt", "a")
        output.write("{}\t{}\t{}\t{}\t{}"
                     .format(contig2, pos2, ref_allele2, all_alleles2, all_freq2))

        stream_pop2_allele_table += ("\n{}\t{}\t{}\t{}\t{}"
                     .format(contig2, pos2, ref_allele2, all_alleles2, all_freq2))

        for sample in record.samples:
            sample_id = record.genotype(sample.sample)
            genotype_base = sample_id.gt_bases

            ## Now, extract the phased allele using the 'PG' tag, not the phased 'GT' tag
            if sample['PG'] != './.':  # which could be anyone of '0/0', 1/1, 1/2, 1|2, 2|2, etc. etc.
                phased_PG = list(sample['PG'])
                first_allele_index = phased_PG[0]
                separator = phased_PG[1]
                second_allele_index = phased_PG[2]

                # Therefore...
                PG_allele = allele_list[int(first_allele_index)] + \
                            separator + allele_list[int(second_allele_index)]

            else:
                PG_allele = sample['PG']

            ## Now, write the remaining values to the file
            try:
                output.write("\t{}\t{}\t{}\t{}\t{}".format(sample['GT'], genotype_base, sample['PG'], PG_allele,
                                                           sample['PI']).replace('./.', '.').replace('None', '.'))
                stream_pop2_allele_table += ("\t{}\t{}\t{}\t{}\t{}".format(sample['GT'], genotype_base, sample['PG'], PG_allele,
                                                           sample['PI']).replace('./.', '.').replace('None', '.'))

            except AttributeError:
                output.write("\t{}\t{}\t{}\t{}\t{}".format(sample['GT'], genotype_base, sample['PG'],
                                                           PG_allele, '.').replace('./.', '.').replace('None', '.'))
                stream_pop2_allele_table += ("\t{}\t{}\t{}\t{}\t{}".format(sample['GT'], genotype_base, sample['PG'],
                                                           PG_allele, '.').replace('./.', '.').replace('None', '.'))

        output.write('\n')
        output.close()


    ## Read the phased vcf file for F1 hybrid (file vcf 03)
    fun_flush_print("#2-B: Mining data from %s ......" % (args.het_vcf))
    stream_f1_vcf = vcf.Reader(open(args.het_vcf, 'r'))
    sample_name = stream_f1_vcf.samples
    sample_size = len(sample_name)

    for record in stream_f1_vcf:
        contig3 = record.CHROM
        pos3 = record.POS
        ref_allele3 = record.REF
        alt_alleles3 = ",".join(map(str, (record.ALT[::])))

        alt_freq3 = ",".join(map(str, (record.INFO['AF'])))
        alt_freq_sum3 = sum(record.INFO['AF'])

        ref_freq3 = round(1 - alt_freq_sum3, 3)
        all_alleles3 = ref_allele3 + "," + alt_alleles3
        all_freq3 = str(ref_freq3) + "," + str(alt_freq3)

        allele_list = all_alleles3.split(',')

        # Now, write some of the values from this line of the vcf file
        output = open(args.output + 'f1_allele_table.txt', "a")
        output.write("{}\t{}\t{}\t{}\t{}"
                     .format(contig3, pos3, ref_allele3, all_alleles3, all_freq3))

        stream_f1_allele_table += ("\n{}\t{}\t{}\t{}\t{}"
                     .format(contig3, pos3, ref_allele3, all_alleles3, all_freq3))

        for sample in record.samples:
            sample_id = record.genotype(sample.sample)
            genotype_base = sample_id.gt_bases

            ## Now, extract the phased allele using the 'PG' tag, not the phased 'GT' tag
            if sample['PG'] != './.':  # which could be anyone of '0/0', 1/1, 1/2, 1|2, 2|2, etc. etc.
                phased_PG = list(sample['PG'])
                first_allele_index = phased_PG[0]
                separator = phased_PG[1]
                second_allele_index = phased_PG[2]

                # Therefore...
                PG_allele = allele_list[int(first_allele_index)] + \
                            separator + allele_list[int(second_allele_index)]

            else:
                PG_allele = sample['PG']

            ## Now, write the remaining values to the file
            try:
                output.write("\t{}\t{}\t{}\t{}\t{}".format(sample['GT'], genotype_base, sample['PG'], PG_allele,
                                                           sample['PI']).replace('./.', '.').replace('None', '.'))

                stream_f1_allele_table += ("\t{}\t{}\t{}\t{}\t{}".format(sample['GT'], genotype_base, sample['PG'], PG_allele,
                                                           sample['PI']).replace('./.', '.').replace('None', '.'))

            except AttributeError:
                output.write("\t{}\t{}\t{}\t{}\t{}".format(sample['GT'], genotype_base, sample['PG'],
                                                           PG_allele, '.').replace('./.', '.').replace('None', '.'))

                stream_f1_allele_table += ("\t{}\t{}\t{}\t{}\t{}".format(sample['GT'], genotype_base, sample['PG'],
                                                           PG_allele, '.').replace('./.', '.').replace('None', '.'))


        output.write('\n')
        output.close()


    fun_flush_print("#2-B: Completed allele_table table for each vcf and samples ......\n\n");

    # Pass the data to another function and then delete it
    return merge_tables(stream_pop1_allele_table, stream_pop2_allele_table, stream_f1_allele_table)


    ## ** Step 02-B complete ** ## ************



## Step 03-A: Read the data frame as table using tab as delimiter
# Note: read_csv and read_table is going to take the first row info as a header by default
def merge_tables(stream_pop1_allele_table, stream_pop2_allele_table, stream_f1_allele_table):
    global args

    fun_flush_print("Stage #3: merging the allele tables from %s, %s and f1-sample %s into one dataframe..."
                    "...using 'contig', 'pos' and 'ref_allele' as keys" % (args.pop1, args.pop2, args.f1_sample));

    ## Below codes are deactivated for now
    #df1 = pd.read_table(args.output + '_' + args.pop1 + "_allele_table.txt", sep='\t')
    #df2 = pd.read_table(args.output + '_' + args.pop2 + "_allele_table.txt", sep='\t')
    #df_f1 = pd.read_table(args.output + 'f1_allele_table.txt', sep='\t')

    fun_flush_print("#3-A: reading allele tables into pandas dataframe.......");

    ## Note: if we desire to read the files directly from console use the below code
    df1 = pd.read_table(StringIO(stream_pop1_allele_table), sep='\t')
    df2 = pd.read_table(StringIO(stream_pop2_allele_table), sep='\t')
    df_f1 = pd.read_table(StringIO(stream_f1_allele_table), sep='\t')

    # Clear memory
    del (stream_pop1_allele_table, stream_pop2_allele_table, stream_f1_allele_table)

    ## Step 03-B: Select only required columns from the allele-table file.

####################  continue .... tomorrow ... introduce indel removal method .......
    #### also introduce read only from a specific chromosome method .......    ###########
    ### write data in landscape vs. potrait mode ....
    ### write variants back to vcf file
    ### see if the variants-table can be in GATK format

    ##### Note: after reading the above files    #### do this later  - to do !!
    ### To do !: Introduce a method to remove the site and it's index if the variant is InDel
    ## only from population vcf not from F1 hybrid vcf   #### do later  !!

##################################################################################################33

    ### only selecting the required columns i.e: sample + _PG_allele and sample + _PI ..
    ### ... from each data frame
    fun_flush_print("#3-B: Selecting only the required column:\n"
                    "1) PG_allele (bases) 2) PI (phase block index of that PG allele) ....... \n\n");

    fun_flush_print("#3-B: from %s population vcf......." % (args.pop1))
    ## for data from population #1
    # set the `contig pos` to index
    df1.set_index(['contig', 'pos'], append=True, inplace=True)

    # and prepare a filter labels to select only required columns
    filter_col = [col for col in df1.columns.values if col.endswith('_PG_allele')
                  or col.endswith('_PI')]


    # filter the data and reset index and add prefix to the sample names
    df1 = df1[filter_col].add_prefix(args.pop1 + '_').reset_index(level=['contig', 'pos'])


    fun_flush_print("#3-B: from %s population vcf......." % (args.pop2))
    ## for data from population #2
    # apply the above method to the data from population #2
    df2.set_index(['contig', 'pos'], append=True, inplace=True)

    filter_col = [col for col in df2.columns.values if col.endswith('_PG_allele')
                  or col.endswith('_PI')]

    df2 = df2[filter_col].add_prefix(args.pop2 + '_').reset_index(level=['contig', 'pos'])

    fun_flush_print("#3-B: from %s population %s sample......." % (args.pop1, args.f1_sample))
    ## for data from F1 hybrids
    # apply the above method to the data from F1 hybrid - a little different though
    df_f1.set_index(['contig', 'pos', 'ref'], append=True, inplace=True)

    # difference is: selecting columns as well as the specific sample
    filter_col = [col for col in df_f1.columns.values if col.startswith(args.f1_sample)
                  and (col.endswith('_PG_allele') or col.endswith('_PI'))]

    # filter the columns, reset the index and rename the columns names
    # this three methods can be done separately if desired
    df_f1 = df_f1[filter_col].reset_index(level=['contig', 'pos', 'ref']) \
        .rename(columns={args.f1_sample + '_PG_allele': 'F1_' + args.f1_sample,
                         args.f1_sample+'_PI': 'F1_' + args.f1_sample + '_PI'})



    fun_flush_print("stage (#3 C): now merging alleles tables for %s, %s & %s\n"
                    % (args.pop1, args.pop2, args.f1_sample));
    ## Step 03-C: Now, merge the data frame using 'contig', 'pos' as the keys...
    ## merge the data frame using 'contig', 'pos' as the keys...
    # ...to output a union of the data frames

    ## Make the list of 3 data frames
    # The below method can be used for multiple dataframes
    data_frames = [df1, df2, df_f1]
    df_merged = reduce(lambda left, right: pd.merge(left, right, on=['contig', 'pos'],
                                                    how='outer', sort=True), data_frames)

    ## Clear memory
    del df1, df2, df_f1, filter_col


    ## select only the data (rows)..
    # .. where the given F1_sample name (e.g 2ms04h) has a meaningful data
    df_merged = df_merged[df_merged['F1_' + args.f1_sample + '_PI'] != '.'].fillna('.')


    fun_flush_print(" Replacing the sample + PG_allele with only sample name\n");
    # replace the string '_PG_allele' from all the column names if present
    # the new column names actually have values from the phased genotypes
    df_merged.rename(columns=lambda x: x.replace('_PG_allele', ''), inplace=True)

    # reset the index after the rows are dropped
    df_merged = df_merged.reset_index(drop=True)


    fun_flush_print(" Writing the merged dataframe to file ..... \n");
    # write the updated dataframe to file
    pd.DataFrame.to_csv(df_merged,
                        args.output + '_' + args.pop1 + '_' + args.pop2 + '_' + args.f1_sample +
                        "_allele_table.txt", sep='\t', index=False)
    ## Note: this output allele_table contains the contig, pos and then genotypes and indexes of the phased genotypes
    # this file is used below for markov chain process.

    fun_flush_print("#3-B: merge complete...\n")

    ## Step 04:
    ## Starting the markov chain process in the merged dataset

    ## Also, we want to convert this dataframe to csv like object to be used downstream
    stream_merged_data = pd.DataFrame.to_csv(df_merged, sep='\t', index=False)
    #fun_flush_print("Next stage (#4): Calculate Odds Ratio to segregate the haplotypes...\n")
    fun_flush_print("Next stage (#4): Starting the 1st order Markov Chain process ...\n")

    ## clear memory
    del data_frames, df_merged

    return markov_chain_process (stream_merged_data)
    #return calc_OddsRatio(stream_merged_data)

    ## end of Step 04


def markov_chain_process (stream_merged_data):
    global args;
    stream_csv = pd.read_table(StringIO(stream_merged_data), sep='\t')

    del stream_merged_data

    fun_flush_print("stage (#5): Count the number of haplotype in each 1st order markov chain ....\n")
    # 05-B: Make a duplicate copy of the `F1 index'
    stream_csv['F1_index'] = stream_csv['F1_' + args.f1_sample + '_PI']

    ## 05-C:
    ## Create a dictionary
    # .. this dictionary will be used downstream to do calculation of odds ratio
    ## convert the counts to dictionary to do Markov model odds ratio calculation (Step # )

    ## select the required columns
    cols = ['contig', 'pos', 'ref', 'F1_index', 'F1_x', 'F1_y']
    # Note: if all columns are desired do this: cols = mcve_data.columns

    # Split the phased genotype/haplotype column into two column
    # one haplotype (on the left is named 'F1_x') and another 'F1_y'
    stream_csv['F1_x'], stream_csv['F1_y'] = zip(*stream_csv['F1_' + args.f1_sample].str.split('|'))

    haplotype_dict = defaultdict(dict)
    haplotype_dict = stream_csv[cols].groupby('F1_index'). \
        apply(lambda x: x.set_index('F1_index').to_dict('list')).to_dict()
    ## Note: this dictionary ('haplotype_dict') will be merged with another dict downstream in the pipeline

    ## eliminating the repeatitive values in the `contig` indentifier
    #  may need to bring it back if needed (when writing the haplotypes back to vcf file)
    for k in haplotype_dict:
        haplotype_dict[k]['contig'] = haplotype_dict[k]['contig'][0]


    # remove the 'F1_x' and 'F1_y' columns and 'ref'
    del stream_csv['F1_x']
    del stream_csv['F1_y']
    del stream_csv['ref']



    # 05-D:
    # -> now group the data using the F1 index values of the F1 genotype/sample ..
    # -> repeat the first and the last line of each block to assist in the begin and end states of a markov chain
    # -> Also, make sure the group.by function doesn't change the order of the lines.
    # The order of the lines are preserved using the 'contig' and 'pos' as sort-keys
    stream_csv = stream_csv.groupby('F1_index', group_keys=False). \
        apply(lambda g: pd.concat([g.head(1), g, g.tail(1)])). \
        sort_values(by=['contig', 'pos'])

    # print('\n mcve after repeat \n')
    # print(mcve_data)

    ## drop the repeat index labels
    # this drops the repeated index label of pandas dataframe
    # but, doesn't drop the repeated line (or row values) of the dataframe
    stream_csv = stream_csv.reset_index(drop=True)

    fun_flush_print("stage (#4): Preparing Markov chain for each haplotype block .....\n")
    ### Step 06: Write a function to identify the possible observed haplotypes states..

    # 06-A:
    # defining a function to do product vs. zip based on condition of the data - important function ..
    # ... based on the index values of the phased genotypes with in population samples

    def find_phase_states(frame):
        # get the labels which do not end with _index
        labels = [(l, l + '_PI')
                  for l in frame.columns.values if not l.endswith('_PI')]

        def separate_rows(row):
            # get row n and row n-1
            front = row[:len(row) >> 1]
            back = row[len(row) >> 1:]

            # loop through the labels
            results = []
            for l, i in labels:
                if '|' or '/' in front[l]:
                    x = list(front[l])[::2]
                else:
                    x = front[l]

                if '|' or '/' in back[l]:
                    y = list(back[l])[::2]
                else:
                    y = back[l]

                if front[i] is not '.' and back[i] is not '.' \
                        and front[i] == back[i]:
                    results.append(x[0] + 'g' + y[0] + ',' + x[1] + 'g' + y[1])

                elif '.' in x or '.' in y:
                    results.append('.')

                else:
                    results.append(
                        ','.join([x + 'g' + y for x, y in it.product(x, y)]))

            return pd.Series(results)

        # take the above function and apply it to pandas dataframe:
        # This runs a first order markov chain to process the data
        df = pd.concat([frame, frame.shift(1)], axis=1)[1:].apply(
            separate_rows, axis=1)

        # put the header information back on the dataframe
        df.rename(columns={i: x[0] for i, x in enumerate(labels)},
                  inplace=True)
        return pd.concat([df], axis=1)

    # 06-B: set the 'contig', 'pos' and 'F1_index' column values to index
    # ... so the find_phase_states(frame) can run properly
    stream_csv_reset = stream_csv.set_index(['contig', 'pos', 'F1_index'], append=True)  # remove F1_index later


    ## 06-C: Take the dataframe to run the Markov-Chain process
    markov_chain_df = (find_phase_states(stream_csv_reset))

    del stream_csv_reset, stream_csv

    # The Markov chain process also creates a row with values that are 4 states for F1 phased genotype
    # this needs to be removed
    markov_chain_df = markov_chain_df[markov_chain_df['F1_' + args.f1_sample].apply(lambda x: len(x.split(',')) < 3)]

    fun_flush_print("stage (#4): Writing markov chain to a file ...\n")
    # write the updated data to a file  # put appropriate names later
    pd.DataFrame.to_csv(markov_chain_df, args.output + '_' + 'markov_chain.txt', sep='\t', index=True)   #- !! keep or remove

    ## 06-D: Now, split the F1 genotypes Markov Chain...
    # ... to separate it into haplotype X vs. Y
    markov_chain_df['F1_x'], markov_chain_df['F1_y'] = zip(*markov_chain_df['F1_' + args.f1_sample].str.split(','))

    # then delete the redundant column
    del markov_chain_df['F1_' + args.f1_sample]  # rename the columns later


    return count_haplotypes(markov_chain_df, haplotype_dict)


def count_haplotypes(df_chain, haplotype_dict):
    global args

    ### Step 07: Now, count the number of phased states from F1 in My vs. Sp samples

    hapX_count = pd.DataFrame(index=df_chain.index[0:0])
    hapY_count = pd.DataFrame(index=df_chain.index[0:0])
    for index, lines in df_chain.iterrows():
        hap_x = lines['F1_x']
        hap_y = lines['F1_y']
        x_count = lines.apply(lambda x: x.split(',').count(hap_x) / 2 if len(x.split(',')) > 2 else x.count(hap_x))
        y_count = lines.apply(lambda x: x.split(',').count(hap_y) / 2 if len(x.split(',')) > 2 else x.count(hap_y))

        hapX_count = hapX_count.append(x_count)
        hapY_count = hapY_count.append(y_count)

    fun_flush_print("stage (#4): Count the number of haplotype in each 1st order markov chain ....\n")
    ## keep this ....
    df_counts = hapX_count.astype(str).add('|').add(hapY_count.astype(str))

    fun_flush_print("stage (#4): Writing the number of haplotype counts to a file ....\n")
    pd.DataFrame.to_csv(df_counts, args.output + '_' + 'phase_counts.txt', sep='\t')


    ## Delete the column that have 'F1' in the hapX_count and hapY_count Dataframe.
    # or do a direct deleting of the column using 'del' function
    for col in hapX_count:
        if 'F1' in col:
            del hapX_count[col]

    for col in hapY_count:
        if 'F1' in col:
            del hapY_count[col]


    fun_flush_print("stage (#4): Sum the counts of haplotype in each 1st order markov chain ....\n")
    ##### summing the count for two population for each haplotype

    ### Note: to keep scripting easy, pop1 is named as 'My' and pop2 as 'Sp'
    ## now separate the number of hapX vs hapY for each population
    # these data can be directly summed if need be using 'pandas.sum(axis=1)`
    hapX_My_count = hapX_count.filter(like=args.pop1 + '_')
    hapX_Sp_count = hapX_count.filter(like=args.pop2 + '_')

    hapY_My_count = hapY_count.filter(like=args.pop1 + '_')
    hapY_Sp_count = hapY_count.filter(like=args.pop2 + '_')

    ## print this to file if need be



    ## now sum the numbers of hapX vs hapY for each population
    # first create a null dataframe
    ## assuming pop1 is My and pop2 is Sp. This doesn't effect the algorithm though
    hapX_My_Sum = pd.DataFrame();
    hapX_Sp_Sum = pd.DataFrame()
    hapY_My_Sum = pd.DataFrame();
    hapY_Sp_Sum = pd.DataFrame()

    # and then sum
    hapX_My_Sum['hapX_My_Sum'] = hapX_My_count.sum(axis=1)
    hapX_Sp_Sum['hapX_Sp_Sum'] = hapX_Sp_count.sum(axis=1)

    hapY_My_Sum['hapY_My_Sum'] = hapY_My_count.sum(axis=1)
    hapY_Sp_Sum['hapY_Sp_Sum'] = hapY_Sp_count.sum(axis=1)

    ## Now, merge the hapX and hapY value for each population ...
    # ... using concatenation
    frames = [hapX_My_Sum, hapY_My_Sum, hapX_Sp_Sum, hapY_Sp_Sum]

    merged_hap_counts = pd.concat(frames, axis=1) \
        .reset_index(level=['contig', 'pos', 'F1_index'])  # .sort_values(by=['contig','pos'])


    ## Clear memory
    del hapX_count, hapY_count, hapX_My_count, hapY_My_count, hapX_Sp_count, hapY_Sp_count


    ### 04 : Merge two dictionaries to start Odds ratio test
    ## convert the counts to dictionary to do Markov model odds ratio calculation

    # select on the index and summed columns for each haplotype-population combination
    # other columns were selected in the earlier dictionary (Step # ... ?? )
    # The resulting dictionary will be merged with the dictionary created earlier in the pipeline
    cols = ['F1_index', 'hapX_My_Sum', 'hapY_My_Sum', 'hapX_Sp_Sum', 'hapY_Sp_Sum']

    haplotype_summed_dict = defaultdict(dict)
    haplotype_summed_dict = merged_hap_counts[cols].groupby('F1_index'). \
        apply(lambda x: x.set_index('F1_index').to_dict('list')).to_dict()
    ## to do: merge this dict with the dict we make early in the pipeline


    ## Now, merge the two haplotype dictionaries: 'haplotype_summed_dict' with 'haplotype_dict'
    # 'F1_index' key will be used to merge the two dictionaries
    # this will update the information in haplotype_dict
    for k, v in haplotype_summed_dict.items():
        for subk, subv in v.items():
            haplotype_dict[k][subk] = subv



    ## Create a file to write the output of the odds ratio
    haplo_output = open(args.f1_sample + "_haplotype_landscape.txt", "w")
    haplo_output.write('contig\thaplotype_block\thap_X\thap_Y\todds_ratio\t' + args.pop1+'_hap\t' + args.pop2+ '_hap\n')
    # Comment: haplotype_block is actually phase index
    haplo_output.close()


    # Creating null haplotype string  ## rename !!!
    My_haplotype = ""
    Sp_haplotype = ""

    # and a variable to stream the landscape mode of the data to potriat mode
    stream_segregated_haplotype = 'contig\tpos\tref\thaplotype_block\thap_X\thap_Y\todds_ratio\t'\
                                  + args.pop1 + '_hap\t' + args.pop2 + '_hap'

    ### to do ! - in the below code we need to account for the number of sample, add pseudo count, etc.
    ## also print the haplotype output in portrait - so it can be imported easily into vcf file
    ## rescue the phase state of the indels using the phase state of the indel
    # fix the zip vs. product between adjacent phases of matching vs. non-matching index label.


    for k in haplotype_dict:
        contig = haplotype_dict[k]['contig']
        posi = (haplotype_dict[k]['pos'])
        print('\n posi \n')
        print(posi)
        print(type(posi))
        print('-'.join(str(x) for x in posi))
        pos = "-".join(str(x) for x in haplotype_dict[k]['pos'])


        ref = "-".join(haplotype_dict[k]['ref'])

        hap_X = "-".join(haplotype_dict[k]['F1_x'])
        hap_Y = "-".join(haplotype_dict[k]['F1_y'])

        # sum the number of haplotypes for alleles within ...
        # ... each haplotype (using the Max-Sum method) for each population
        ## Note. Max-Sum vs. Max-Product give same results, when taking ratios of the ratio

        hapX_prob_My = sum(float(x) for x in (haplotype_dict[k]['hapX_My_Sum'])) + 1
        hapY_prob_My = sum(float(x) for x in (haplotype_dict[k]['hapY_My_Sum'])) + 1

        hapX_prob_Sp = sum(float(x) for x in (haplotype_dict[k]['hapX_Sp_Sum'])) + 1
        hapY_prob_Sp = sum(float(x) for x in (haplotype_dict[k]['hapY_Sp_Sum'])) + 1

        # calculate the odds-ratio
        OddsRatio_hX_vs_hY_is_MyVsSp = Odds_Ratio = ((hapX_prob_My / hapX_prob_Sp) / (hapY_prob_My / hapY_prob_Sp))

        del OddsRatio_hX_vs_hY_is_MyVsSp

        if Odds_Ratio > 1.0:
            My_haplotype = hap_X
            Sp_haplotype = hap_Y
            print("the first haplotype %s is Mayodan" % str(hap_X))
        elif Odds_Ratio < 1.0:
            My_haplotype = hap_Y
            Sp_haplotype = hap_X
            # print("the first haplotype %s is Spiterstulen" % str(hap_X))

        elif Odds_Ratio == 1.0:
            My_haplotype = ""  # use a replace method to sub haplotype allele with Ns
            Sp_haplotype = ""

        Odds_Ratio = round(Odds_Ratio, 3)

        ## Now, write the index and phase type into a file
        output_hap = open(args.f1_sample + "_haplotype_landscape.txt", "a")
        output_hap.write(str(contig) + '\t' + str(k) + '\t' + hap_X + '\t' + hap_Y + '\t' + str(
            Odds_Ratio) + '\t' + My_haplotype + '\t' + Sp_haplotype + '\n')


        # writing haplotype in potrait mode:
        stream_segregated_haplotype += '\n'+ str(contig) + '\t' + str(pos) +  '\t' + str(ref) + '\t' + str(k) + '\t' + hap_X + '\t' + hap_Y + '\t' + str(
            Odds_Ratio) + '\t' + My_haplotype + '\t' + Sp_haplotype

    return haplotype_in_potrait(stream_segregated_haplotype)

    output_hap.close()


from itertools import chain
import numpy as np

def haplotype_in_potrait(stream_segregated_haplotype):
    global args

    landscape_df = df = pd.read_table(StringIO(stream_segregated_haplotype), sep='\t')

    del landscape_df # was just added for naming purpose

    ## create list by split
    cols = ['pos', 'ref', 'hap_X', 'hap_Y', args.pop1 +'_hap', args.pop2 +'_hap']
    df[cols] = df[cols].apply(lambda x: x.str.split('-'))


    ## flatten the lists and repeat the values where the column values are same ...
    ## for e.g columns like 'contig', 'haplotype_block', odds-ratio

    length = df.ref.str.len() # can use any filled columns (ref, hapX, hapY)
    ## for calculating the length

    potrait_data = pd.DataFrame({"contig": np.repeat(df.contig.values, length),
        "pos": list(chain.from_iterable(df.pos.values)),
        "ref": list(chain.from_iterable(df.ref.values)),
        "haplotype_block": np.repeat(df.haplotype_block.values, length),
        "hap_X": list(chain.from_iterable(df.hap_X)),
        "hap_Y": list(chain.from_iterable(df.hap_Y)),
        "odds_ratio": np.repeat(df.odds_ratio.values, length),
        args.pop1 + "_hap": list(chain.from_iterable(df[args.pop1 + '_hap'])),
        args.pop2 + "_hap": list(chain.from_iterable(df[args.pop1 + '_hap']))
                                 }).reindex_axis(df.columns, axis=1)

    ## correct the repeatition of the duplicates in the column (odds ratio)
    s = pd.Series(np.repeat(df.index.values, length))
    potrait_data.loc[s.duplicated(), 'odds_ratio'] = '-'

    ## clear memory
    del cols, length, s

    pd.DataFrame.to_csv(potrait_data, args.f1_sample + "_haplotype_potrait.txt", sep='\t', index=False)


    print("**end**")



# we are creating a new function called "fun_flush_print" which is/will-be used for working with several other functions within the program
def fun_flush_print(text):
    print(text);
    sys.stdout.flush();


def fatal_error(text):
    fun_flush_print("     FATAL ERROR: " + text);
    quit();


if __name__ == "__main__":
    main();


