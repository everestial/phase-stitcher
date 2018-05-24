#!/home/bin/python3

# Note: this program uses python 3
# Function of the program: preparing a interactive version of the pHASE-Stitcher (using 1st order Markov chain)

## Checking required modules
print("Checking and importing required modules... ")


import argparse
from collections import OrderedDict
from decimal import Decimal
import os, time
from io import StringIO
import itertools
from itertools import product
from multiprocessing import Pool
import pandas as pd
import resource
import re
import shutil
import sys


'''Overall structure of the program workflow:
    01 > Load haplotype file with maternal, paternal and F1 hybrid data.
    02 > Assign sample names to appropriate genotype background.
    03 > Split data by chromosome and write into multiple files for multiprocessing.
    04 > Now, use imap (lazy mapping) of each split of data.
    05 > Group the data by unique "PI" index and then pass it to marko process.
    06 > Compute the likelihood of each haplotype (left and right) in F1 RBphased blocks
         beloging to maternal vs. paternal background.
    07 > Based on the lods2CutOff score segregate the haplotype into maternal paternal background.
'''


def main():

    print()
    ''' printing authorship '''

    print("#######################################################################")
    print("        Welcome to phase-Stitcher version %s       " % 1.1)
    print("  Author: kiran N' bishwa (bkgiri@uncg.edu, kirannbishwa01@gmail.com) ")
    print("#######################################################################")
    print()


    print("Loading the argument variables and assigning values ....")

    ''' Define required argument for interactive mode program. '''
    parser = argparse.ArgumentParser()

    parser.add_argument("--nt",
                        help="number of process to run -> "
                             "The maximum number of processes that can be run at once is the number "
                             "of different chromosomes (contigs) in the input haplotype file.",
                        default=1, required=False)

    parser.add_argument("--input",
                        help="name of the input haplotype file -> "
                             "This haplotype file should contain unique index represented by 'PI' and "
                             "phased genotype represented by 'PG_al' for all the samples. ",
                        required=True)

    parser.add_argument("--pat",
                        help="Paternal sample or comma separated sample names that "
                             "belong to Paternal background. "
                             "Sample group may also be assigned using prefix. "
                             "Options: 'paternal sample name', 'comma separated samples', 'pre:...'. "
                             "Unique prefix (or comma separated prefixes) should begin with 'pre:'. ",
                        required=True)

    parser.add_argument("--mat",
                        help="Maternal sample or sample names (comma separated) that "
                             "belong to maternal background. Sample group can also be "
                             "assigned using unique prefix/es. "
                             "Options: 'maternal sample name', 'comma separated samples', 'pre:...'. "
                             "Unique prefix (or comma separated prefixes) should begin with 'pre:'. ",
                        required=True)

    parser.add_argument("--f1Sample",
                        help="Name of the F1-hybrid sample. Please type the name of only one F1 sample.",
                        required=True)

    parser.add_argument("--outPatMatID",
                        help="Prefix of the 'Paternal (dad)' and 'Maternal (mom)'genotype in the output file. "
                             "This should be a maximum of three letter prefix separated by comma. "
                             "Default: 'pat,mat'.",
                        default='pat,mat', required=False)

    parser.add_argument("--output", default='', type=str,
                        help="Name of the output directory. "
                             "Default: f1SampleName + '_stitched' ", required=False)

    parser.add_argument("--lods",
                        help="log(2) odds cutoff threshold required "
                             "to assign maternal Vs. paternal haplotype segregation and stitching. ",
                        default=5)
    # note: log_odds_ratio 5 is 2^5 = 32 times likely

    parser.add_argument("--culLH",
                        help="Cumulative likelhood estimates -> "
                             "The likelhoods for haplotype segregation can either be max-sum vs. max-product. "
                             "Default: maxPd i.e max-product. "
                             "Options: 'maxPd' or 'maxSum'. ",
                        default='maxPd', required=False)

    parser.add_argument("--chr", required = False, default="",
                        help="Restrict haplotype stitching to a specific chromosome.")

    parser.add_argument("--hapStats",
                        help="Computes the descriptive statistics of final haplotype. "
                             "Default: 'no'."
                             "Option: 'yes', 'no' .",
                        default='no', required=False)



    ''' create a global argument variables and declare values of the global variables.
            The values for the global variable are assigned using arguments (args.()) passed by user.'''
    # this is necessary if arguments are declared at several functions using "args.variable" parameter.
    global args;
    args = parser.parse_args()  # .. but keep it as it is.

    global output  # ** now used with "outputdir" - may change in future to control filename in output


    print("#Step 01:  Checking the required files......\n");

    #Step 01-B checks if the required files (vcfs) are provided or not
    check_files = [args.mat, args.pat, args.f1Sample, args.input];
    for xfile in check_files:
        if xfile == "":
            print('fatal error !!! ')
            print("Please assign the required parameter : '%s' ." % (xfile));
                # reports error message if the checked files isn't found
            sys.exit()



    ################## begins at ............
    global output
    global outputdir
    global use_bed  # ** prolly deprecate
    global soif1
    global soimom
    global mom_id
    global soidad
    global dad_id
    global snp_threshold
    global num_of_hets
    global lods_cut_off
    global maxed_as
    global hapstats
    global time01  # set a start time (used to time computational run for each chromosome vs. all processes)


    print("Assigning values to the global variables ....")
    # advantage of this method (assignment of global variable values at the beginning) ...
    # ... is that we can change the program to non-interactive mode easily and
    # point to specific values/files on this top part. This becomes helpful while debugging.

    input_file = args.input  # the input haplotype file
    # input_file = 'allele_table_for_phase_extender.txt'
    print('  - using haplotype file "%s" ' % (input_file))

    time01 = time.time()
    soif1 = args.f1Sample
    print('  - F1-hybrid of interest: "%s" ' % (soif1))


    # only read the header of the input haplotype file and ...
    # .. find the required samples (for maternal and paternal background)
    header_ = open(input_file, 'r').readline()
    print(header_)

    soimom = args.mat
    soimom = find_samples(soimom, header_)

    soidad = args.pat
    soidad = find_samples(soidad, header_)



    # prefix for mat,pat haplotype id
    # ** by default treats the samples in "--mat" as maternal and "--pat" as paternal
    pat_mat_id = args.outPatMatID
    if pat_mat_id == 'pat,mat':
        dad_id = 'pat_hap'
        mom_id = 'mat_hap'
    else:
        dad_id = pat_mat_id.split(',')[0] + '_hap'
        mom_id = pat_mat_id.split(',')[1] + '_hap'



    # Assign the output directory
    if args.output == '':
        outputdir = soif1 + '_stitched'
    elif args.output != '':
        outputdir = args.output
    if os.path.exists(outputdir):
        shutil.rmtree(outputdir, ignore_errors=False, onerror=None)
    os.makedirs(outputdir, exist_ok=True)


    if args.chr != "":
        chr_list = (args.chr).split(',')
    else:
        chr_list = ""


    # assign number of process to be used
    nt = int(args.nt)  # default, nt = 1
    print('  - using "%s" processes ' % (nt))


    lods_cut_off = int(args.lods)  # log_of_odds_cut_off, default = 5
    print('  - using log2 odds cut off of "%s" ' % (lods_cut_off))


    # add argument for max sum vs. max product of likelyhood estimates before calculating the LOD-score
    maxed_as = args.culLH  # default, maxed_as = "*"
    if maxed_as == 'maxSum':
        max_is = 'max sum'
        maxed_as = '+'
    elif maxed_as == 'maxPd':
        max_is = 'max product'
        maxed_as = '*'
    print('  - using "%s" to estimate the cumulative maximum likelyhood while segregating '
          'the diploid haplotype block into maternal vs. paternal haplotype ' % (max_is))


    # Add argument to compute descriptive statistics of the final
    # print the hapstats to file and also plot histogram
    if args.hapStats == 'yes':  # default, hapstats = 'no' ** ??
        hapstats = 'yes'
        print('  - statistics of the haplotype before and after extension will '
              'be prepared for the sample of interest i.e "%s" ' %(soif1))
    else:
        hapstats = 'no'
        print('  - statistics of the haplotype before and after extension will not '
          'be prepared for the sample of interest i.e "%s". '
              '    Only extendent haplotype block will be prepared.' % (soif1))

    print('hapstats :', hapstats)



    '''Assign the number of process.
       **note: number of process should be declared after all the global variables are declared,
       because each pool will need to copy the variable/value of global variables. '''
    pool = Pool(processes=nt)  # number of pool to run at once; default at 1

    #### completed assignment of values to argument variables


    ''' Step 01: Read the input file and and prepare two files as output '''
    # a) One output file contains extended phase-states for the sample of interest (soi)
    # b) another output file contains the lines that have missing data for sample of interest

    ''' Step 01 - A: read the input haplotype file and prepare output files '''
    with open(input_file) as input_data,\
            open(outputdir + '/' + soif1 + '_missing_data.txt', 'w+') as missing_data, \
            open(outputdir + '/' + soif1 + '_haplotype_long.txt', 'w+') as output_long, \
            open(outputdir + '/' + soif1 + '_haplotype_wide.txt', 'w+') as output_wide:

        print()
        print('Reading the input haplotype file. "%s" ' % input_data.name)


        # load data to pandas dataframe
        my_df = pd.read_csv(input_data, sep='\t')


        ''' write header for the output data (data with missing PI for soif1, stitched haplotype (long and wide)).
            - rest of the data will be appended downstream in the pipeline. '''

        missing_data.write('\t'.join(list(my_df)) + '\n')
        output_long.write('\t'.join(['CHROM', 'POS', 'REF', 'all-alleles',
                                     soif1 + ':PI', soif1 + ':PG_al', 'log2odds',
                                     dad_id, mom_id]) + '\n')
        output_wide.write('\t'.join(['CHROM', 'POS_Range', soif1 + ':PI', 'hap_left',
                                     'hap_right', 'log2odds', dad_id, mom_id]) + '\n')


        ''' using pandas dataframe to load the data as df and dictionary and group by contig'''
        # ** Note: Possible future optimization using "dask" module
        my_df_by_contig = my_df.groupby('CHROM', sort=False)


        '''Step 01 - B: Starting multiprocessing after splitting the dataframe by chromosome.'''
        # multiprocessing of data (by lazy "imap" process)

        # folder for storing the splitted dataframes
        if os.path.exists('chunked_Data_' + soif1):
            shutil.rmtree('chunked_Data_' + soif1, ignore_errors=False, onerror=None)
        os.makedirs('chunked_Data_' + soif1, exist_ok=True)


        # split the data by grouping and write it to the disk.
        for chr_, data_by_chr in my_df_by_contig:
            if chr_list == '':
                pd.DataFrame.to_csv(data_by_chr, 'chunked_Data_' + soif1 + '/mydata_' + str(chr_),
                                    sep='\t', index=False, header=True)
            elif chr_list != '':
                if str(chr_) in chr_list:
                    pd.DataFrame.to_csv(data_by_chr, 'chunked_Data_' + soif1 + '/mydata_' + str(chr_),
                                        sep='\t', index=False, header=True)



    '''Step 02 - A: Now, pipe the procedure to next function for multiprocessing (i.e Step 02 - B) '''
    multiproc(pool)


    '''** workflow, between step (02-A and 06) is passed onto other functions. '''


    '''Step 06: Clean the directory that is storing splitted dataframe. '''

    # remove the chunked data folder
    #shutil.rmtree('chunked_Data_' + soif1, ignore_errors=False, onerror=None)

    print('The End :)')




'''Step 02 - B : Pipe the data to multiprocessing environment. '''
def multiproc(pool):

    print()
    time_all = time.time()


    ''' Step 02 - B: Start, multiprocessing/threading - process each contig separately. '''
    path = 'chunked_Data_' + soif1 + '/'  # opens the a temp folder
    file_path = sorted([path + item for item in os.listdir(path)], key=numericalSort)


    '''This imap multiprocess take the pipeline to Step 03. '''
    results = pool.imap(process_by_contig, file_path)

    pool.close()
    pool.join()
    pool.terminate()


    # ** deprecated: if no returns desired when appending dataframe directly on fly
    # pool.imap(process_by_contig, file_path)


    ###################################################################################
    ##############  codes below are reactivated again.  ******************
    ###### In future we can decided to append the data directly to a file ..
    ## .. rather than returning it to a imap (pool) - it might be better ??
    ## this was avoided because there was extra step, time and memory consumption while taking the list of dataframes ..
      # .. splitting it and then writing the long vs. wide format separately.
    # but, will that lead to a problem when multiple data/process are being written or are in imap queue??
    # ** rethink in the future if this is a problem with large data, and while writing to the disk.


    results = list(results)

    result_long = [result[0] for result in results]
    result_wide = [result[1] for result in results]


    ## merge all the data both length and widthwise
    result_merged_length = pd.concat(result_long)
    result_merged_width = pd.concat(result_wide)


    ''' Step : merge the returned result (extended haplotype block by each contigs)
            #together and write it to an output file. '''

    # write this updated dataframe to file.
    # ** for future: or write is directly after "imap'ing". The main problem is that they complete at different times
    # and complete result is not obtained. ** problem not known.
    pd.DataFrame.to_csv(result_merged_length, outputdir + '/' + soif1 + '_haplotype_long.txt',
                        header=True, index=False, sep='\t')

    pd.DataFrame.to_csv(result_merged_width, outputdir + '/' + soif1 + '_haplotype_wide.txt',
                        header=True, index=False, sep='\t')


    print()
    print("Completed haplotype segregation and stitching for all the chromosomes.\n"
          "Time elapsed: '%f' sec. " % (time.time() - time_all))
    print(' - Global maximum memory usage: %.2f (mb)' % current_mem_usage())
    print()

    print("Completed writing the dataframes .....")
    print()


    # ''' Prepare the descriptive statistics of the phased (stitched) haplotypes if desired: '''
    if hapstats == 'yes':
        with open(outputdir + '/' + soif1 + '_haplotype_stats.txt', 'w+') as out_stats:
            stats = result_merged_width.groupby(['CHROM'], sort = False)
            out_stats.write('\t'.join(['CHROM', 'phasedBlock', 'unphasedBlock', 'numVarsInPhasedBlock',
                            'numVarsInUnPhasedBlock', 'log2oddsInPhasedBlock',
                            'log2oddsInUnPhasedBlock',	'totalNumOfBlock',	'totalNumOfVars']) + '\n')

            for xath, valth in stats:
                print('Computing descriptive statistics of the final haplotype for contig "%s". ' %str(xath))

                #distribution of the log2odds values in that CHROM
                log2odds_dist = valth['log2odds'].tolist()

                # name of all haplotype blocks (i.e PI) in that CHROM
                hapblocks = valth[soif1 + ':PI'].tolist()

                # among all the haplotype blocks find the the ..
                # .. indexes of the ones that are phased vs. unphased in final haplotype
                idx_ofphasedblocks = [ith for ith in range(len(log2odds_dist))
                                      if abs(log2odds_dist[ith]) >= lods_cut_off]
                idx_ofunphasedblocks = [ith for ith in range(len(log2odds_dist))
                                        if abs(log2odds_dist[ith]) < lods_cut_off]

                # based on that index find the phased blocks vs. unphased blocks
                phasedblocks = [hapblocks[zth] for zth in idx_ofphasedblocks]
                unphasedblocks = [hapblocks[zth] for zth in idx_ofunphasedblocks]

                # compute number of variants in all vs. phased vs. unphased blocks
                num_vars_in_allblocks = [len(ith.split('-')) for ith in valth['hap_left'].tolist()]
                num_vars_phased_blocks = [num_vars_in_allblocks[zth] for zth in idx_ofphasedblocks]
                num_vars_unphased_blocks = [num_vars_in_allblocks[zth] for zth in idx_ofunphasedblocks]

                # convert values to appropriate structure and string type before writing
                phasedblocks = ','.join(str(x) for x in phasedblocks)
                unphasedblocks = ','.join(str(x) for x in unphasedblocks)
                num_vars_phased_blocks = ','.join([str(x) for x in num_vars_phased_blocks])
                num_vars_unphased_blocks = ','.join([str(x) for x in num_vars_unphased_blocks])

                # find the log2odds of the phased vs. unphased blocks
                log2_phased_blocks = [log2odds_dist[zth] for zth in idx_ofphasedblocks]
                log2_unphased_blocks = [log2odds_dist[zth] for zth in idx_ofunphasedblocks]

                log2_phased_blocks = ','.join([str(x) for x in log2_phased_blocks])
                log2_unphased_blocks = ','.join([str(x) for x in log2_unphased_blocks])

                # compute total number of blocks and total number of variants
                count_blocks = str(len(log2odds_dist))
                count_vars = str(sum(num_vars_in_allblocks))


                data_to_write = [str(xath), phasedblocks, unphasedblocks, num_vars_phased_blocks,
                                 num_vars_unphased_blocks, log2_phased_blocks, log2_unphased_blocks,
                                 count_blocks, count_vars]
                data_to_write = ['.' if x == '' else x for x in data_to_write]

                out_stats.write('\t'.join(data_to_write))

                out_stats.write('\n')
            print()

    elif hapstats == 'no':
        print('Proceeding without preparing descriptive statistics for final haplotype.')



'''Step 03: Process each contig separately with lazy mapping using "imap" multiprocessing. '''
def process_by_contig(file_path):

    #print()
    time_chr = time.time()

    ''' Step 03 - A: Process data for each chromosome (contig) separately. '''
    good_data_by_contig = open(file_path, 'r')
    chr_ = good_data_by_contig.name.split('_')[-1]
    # this name identification might cause problem if "chr" has "_" in it's name i.e scaff_18

    print()
    print('Starting markov transition for chromosome "%s" ' % (chr_))
    print('##########################')


    ## Load the data into pandas dataframe and again groupby unique "PI" values
    contigs_group = pd.read_csv(StringIO(good_data_by_contig.read()), sep='\t')

    # groupby "PI index" and feed it to transition computation
    my_df_grouped = contigs_group.groupby(soif1 + ':PI', sort=False)

    # empty list to store the pandas dataframe
    haplotype_result_long = []
    haplotype_result_wide = []


    '''Step 03 - B: process the data within each unique "PI" separately '''
    for pi_index, data_by_pi in my_df_grouped:
        if pi_index == '.':
            pd.DataFrame.to_csv(data_by_pi, outputdir + '/' + soif1 + '_missing_data.txt',
                                header=False, index=False, mode='a', sep='\t')

        else:
            ## pipe the data for computing transition and segregating the haplotypes
            my_df_dict = data_by_pi.to_dict(orient='list')

            ### Split the RBphase block into the left and right haplotype type
            haplotype_left = [x.split('|')[0] for x in my_df_dict[soif1 + ':PG_al']]
            haplotype_right = [x.split('|')[1] for x in my_df_dict[soif1 + ':PG_al']]


            ''' Enter Step 04: '''
            ## Pass the data to a function to compute likelihood of each configuration
              # - estimate likelihood that each haplotype (left vs. right) belong to mom vs. dad.
              # - estimate likelihood using both forward and reverse markov chain.
              # - ** for future (likelihood estimates can be multiprocessed)


            '''Step 04 - A: Estimate likelihood on both forward and reverse chain.'''

            ''' likelihood of a haplotype (left vs. right) belonging to a parent (dad). '''
            # on forward chain
            likelihood_hap_left_dad_fwd, likelihood_hap_right_dad_fwd = \
                compute_transition(my_df_dict, haplotype_left, haplotype_right, parent='dad', orientation=lambda x: x)

            # on reverse chain
            likelihood_hap_left_dad_rev, likelihood_hap_right_dad_rev = \
                compute_transition(my_df_dict, haplotype_left, haplotype_right, parent='dad', orientation=reversed)



            ''' likelihood of a haplotype (left vs. right) belonging to (mom) or maternal background. '''
            # on forward chain
            likelihood_hap_left_mom_fwd, likelihood_hap_right_mom_fwd = \
                compute_transition(my_df_dict, haplotype_left, haplotype_right, parent='mom', orientation=lambda x: x)

            # on reverse chain
            likelihood_hap_left_mom_rev, likelihood_hap_right_mom_rev = \
                compute_transition(my_df_dict, haplotype_left, haplotype_right, parent='mom', orientation=reversed)



            '''Step 04 - B: Maximize the likelihood estimate (using maxSum or maxProduct). '''
            likelihood_hap_left_dad = cumulate_likelihoods([
                likelihood_hap_left_dad_fwd, likelihood_hap_left_dad_rev])

            likelihood_hap_right_dad = cumulate_likelihoods([
                likelihood_hap_right_dad_fwd, likelihood_hap_right_dad_rev])


            likelihood_hap_left_mom = cumulate_likelihoods([
                likelihood_hap_left_mom_fwd, likelihood_hap_left_mom_rev])

            likelihood_hap_right_mom = cumulate_likelihoods([
                likelihood_hap_right_mom_fwd, likelihood_hap_right_mom_rev])



            '''Step 04 - B: Segregate the haplotypes: Pass the likelihood estimates to another function
                            to compute log2Odds and segregate the haplotypes into mom vs. dad background. '''
            updated_df_long, updated_df_wide = compute_lods(likelihood_hap_left_dad, likelihood_hap_right_dad,
                         likelihood_hap_left_mom, likelihood_hap_right_mom, haplotype_left, haplotype_right,
                                                            data_by_pi)


            ### Store the pandas dataframe as list
            haplotype_result_long.append(updated_df_long)
            haplotype_result_wide.append(updated_df_wide)


    if len(haplotype_result_long) == 0 or len(haplotype_result_wide) == 0:
        result_merged_long = pd.DataFrame()
        result_merged_wide = pd.DataFrame()
        result_merged_stats = pd.DataFrame()


    '''Returns to Step 06: '''
    ## merge the list of pandas dataframe into one dataframe
        # and then return it to the imap pool
        # ** - we can directly write the data to file using pd.csv (append) method
    try:
        result_merged_long = pd.concat(haplotype_result_long, axis=0)
        result_merged_wide = pd.concat(haplotype_result_wide, axis=0)

    except ValueError:  # Sometimes empty list is returned so, raising this exception
        print('No RBphased haplotypes available in this contig.')
        result_merged_long = pd.DataFrame()
        result_merged_wide = pd.DataFrame()



    print(' - Phase-extension completed for contig "%s" in %.6f seconds' % (chr_, time.time() - time_chr))
    print(' - Worker maximum memory usage: %.2f (mb)' % (current_mem_usage()))
    print()


    return (result_merged_long, result_merged_wide)


    ### Deprecated for now !!!
    ''' Write the result: now, write the segregated haplotype to a output file.
        In both long and wide format.
        **Note: I decided to append each returned dataframe directly to output file, rather
        than passing it back to "imap" pool. I think this is more efficient (both time and memory). '''


    ## write results for long format
    #pd.DataFrame.to_csv(result_merged_long, outputdir + '/' + soif1 + '_haplotype_long.txt',
    #                    header=False, index=False, mode='a', sep='\t')

    ## write results for wide format
    #pd.DataFrame.to_csv(result_merged_wide, outputdir + '/' + soif1 + '_haplotype_wide.txt',
    #                    header=False, index=False, mode='a' ,sep='\t')




'''Step 04 - A : Estimate the likelihood of a haplotype being mom vs. dad. This function needs inputs:
    - RBphased data (as my df dict), left vs. right haplotype for soif1, parents, orientation (fwd vs. rev). '''
def compute_transition(my_df_dict, haplotype_left, haplotype_right, parent, orientation):

    #print()
    #print('parent: ', parent)
    #print('orientation: ', orientation)
    ## pick parental type for likelihood estimation
    if parent == 'dad':
        parent  = soidad

    elif parent == 'mom':
        parent = soimom


    '''Create a variable that stores likelihood values, that ..
       left vs. right haplotype belongs to a particular parent.
       This depends upon which parent was passed in this function'''

    ''' We assign the prior likelihood that left and right haplotype transition equally belong to that parent.
        This likelihood estimates then cumulates (maxSum or maxPd) as we get more evidence from emissions and
        transitions frequency. '''
    ## likelihood = prob(for any random emission) x prob(of any random transition) = 0.5 * 0.0625
    ## **Note: different level of value (prior) doesn't affect maxPd but only maxSum.
    likelihood_hap_left = Decimal(0.25 * 0.0625)  # lambda x: x
    likelihood_hap_right = Decimal(0.25 * 0.0625)  # lambda x: x


    ## Deprecated !!!
    # if maxed_as == '+':
    # likelihood_hap_left = 0.25*0.0625
    # likelihood_hap_right = 0.25*0.0625

    # if maxed_as == '*':
    # likelihood_hap_left = 1
    # likelihood_hap_right = 1



    '''Step 04 - A (i): Run 1st order markov chain on one block of unique "PI" data. '''
    for v1, v2 in zip(orientation(list(enumerate(my_df_dict[soif1 + ':PG_al']))),
                      itertools.islice(orientation(list(enumerate(my_df_dict[soif1 + ':PG_al']))), 1, None)):

        # if reverse chain isn't considered, we can use this:
        #for v1, v2 in zip(enumerate(my_df_dict[soif1 + '_PG_al']),
        #                  itertools.islice(enumerate(my_df_dict[soif1 + '_PG_al']), 1, None)):

        n1 = v1[0]
        n2 = v2[0]



        '''Skip the markov chain if indels is present in either "n1" or "n2" level. '''
        if len(v1[1]) > 3 or len(v2[1]) > 3:
            continue

        '''Skip the markov chain if symbolic allele (e.g "*") is present in either "n1" or "n2" level. '''
        if "*" in v1[1] or "*" in v2[1]:
            continue



        ''' Step: Estimate emission counts and emission probabilities at "n1" level. '''

        # empty dict to store emission data
          # using 0.25 as prior probability
        #nucleotide_count_dict = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        nucleotide_count_dict = {'A': 0.25, 'T': 0.25, 'G': 0.25, 'C': 0.25}

        nucleotide_prob_dict = nucleotide_count_dict.copy()


        # using for loop over the given sample to compute emission counts
        for (x, y) in parent:
            for nucleotide in 'ATGC':
                # count only after splitting. else 'AT|C' will have one count of 'A'
                nucleotide_count_dict[nucleotide] += my_df_dict[y][n1].split('|').count(nucleotide)


        # now, compute emission probabilities
        total_nucl = sum(list(nucleotide_count_dict.values()))

        for element in nucleotide_prob_dict:
            nucleotide_prob_dict[element] = \
                compute_probs(nucleotide_count_dict[element], total_nucl, modeis='emission')

        #print('nucleotide probs')
        #print(nucleotide_prob_dict)
        #print()



        ''' Step : Estimate transition counts and probabilities from "n1" to "n2". '''

        # empty dict to store transition data
          # here we use 1/16 as a prior probability distribution of each transition under random assumption
          # this is supplied as pseudo estimate and serves as minimal observation (probability) when there is no data ..
          # .. for the given transition
        transition_count_dict = {
            ('A', 'A'): 1/16, ('A', 'T'): 1/16, ('A', 'G'): 1/16, ('A', 'C'): 1/16,
            ('T', 'A'): 1/16, ('T', 'T'): 1/16, ('T', 'G'): 1/16, ('T', 'C'): 1/16,
            ('G', 'A'): 1/16, ('G', 'T'): 1/16, ('G', 'G'): 1/16, ('G', 'C'): 1/16,
            ('C', 'A'): 1/16, ('C', 'T'): 1/16, ('C', 'G'): 1/16, ('C', 'C'): 1/16 }

        # shallow copy should be fine  # **for future - do deep copy if need be
        transition_prob_dict = transition_count_dict.copy()


        ## Estimate transition counts
        for x, y in parent:
            nucl_n1 = (my_df_dict[y][n1]).split('|')  # nucleotides at n1 level
            nucl_n2 = (my_df_dict[y][n2]).split('|')  # nucleotides at n2 level

            # create zip (n1 to n2) haplotype if they have same "PI" index
            if my_df_dict[x][n1] == my_df_dict[x][n2]:
                ziped_nuclb1b2 = list(zip(nucl_n1, nucl_n2))

                # then count the number of transitions
                for (from_, to) in transition_count_dict:
                    transition_count_dict[(from_, to)] += ziped_nuclb1b2.count((from_, to))

            # create prod(n1 to n2) haplotype if they have non matching "PI" index
            # ** for future - prod can also be used if the PI index is "." i.e has non phased SNPs.
            elif my_df_dict[x][n1] != my_df_dict[x][n2]:
                prod_nuclb1b2 = list(product(nucl_n1, nucl_n2))

                # then count the number of transitions
                for (from_, to) in transition_count_dict:
                    transition_count_dict[(from_, to)] += prod_nuclb1b2.count((from_, to))/2


        ## convert transition counts to probabilities
        for (from_, to) in transition_prob_dict:
            transition_prob_dict[(from_, to)] = \
                compute_probs(transition_count_dict[(from_, to)], nucleotide_count_dict[from_], modeis='transition')

        #print('transition counts')
        #print(transition_count_dict)
        #print()

        #print('transition probs dict')
        #print(transition_prob_dict)
        #print()



        ''' find observed transition configuration for left and right
        haplotype from "n-1" to "n" position'''

        hap_transition_left = (haplotype_left[n1], haplotype_left[n2])
        hap_transition_right = (haplotype_right[n1], haplotype_right[n2])

        #print('transitions are: ')
        #print('left: ', haplotype_left[n1], 'to', haplotype_left[n2])
        #print('right: ', haplotype_right[n1], 'to', haplotype_right[n2])
        #print()


        ''' Now, cumulate the likelihood (as cumulation of the (emission prob * transition prob)). '''
        ## cumulation of the probability can be either maxSum or maxProduct.

        likelihood_hap_left = cumulate_likelihoods(
            [likelihood_hap_left,
             nucleotide_prob_dict[haplotype_left[n1]] * transition_prob_dict[hap_transition_left]])

        likelihood_hap_right = cumulate_likelihoods(
            [likelihood_hap_right,
             nucleotide_prob_dict[haplotype_right[n1]] * transition_prob_dict[hap_transition_right]])

        #print('likelihoods left and right')
        #print(likelihood_hap_left, likelihood_hap_right)


    return likelihood_hap_left, likelihood_hap_right



    ##########################################################
    ###################### deprecated method !!
        #if maxed_as == '+':
    #    likelihood_hap_left += \
    #           nucleotide_prob_dict[haplotype_left[n1]] * transition_prob_dict[hap_transition_left]

    #       likelihood_hap_right += \
    #           nucleotide_prob_dict[haplotype_right[n1]] * transition_prob_dict[hap_transition_right]


    #        elif maxed_as == '*':
    #       likelihood_hap_left *= \
    #           nucleotide_prob_dict[haplotype_left[n1]] * transition_prob_dict[hap_transition_left]
    #
    #       likelihood_hap_right *= \
    #           nucleotide_prob_dict[haplotype_right[n1]] * transition_prob_dict[hap_transition_right]
    ##############################################################



'''Step 04 - C: Estimate likelihood ratio and then segregate haplotypes. '''
def compute_lods(likelihood_hap_left_dad, likelihood_hap_right_dad,
                 likelihood_hap_left_mom, likelihood_hap_right_mom,
                 haplotype_left, haplotype_right, data_by_pi):


    kite = 'hello'  # just a mock variable to highlight the doc string below

    '''Segregate the haplotype and write it to a file. '''
    lh_hapL_vs_hapR_is_pat_vs_mat = Decimal(
        (likelihood_hap_left_dad / likelihood_hap_left_mom) /
        (likelihood_hap_right_dad / likelihood_hap_right_mom))

    log_Odds_Ratio = Decimal(lh_hapL_vs_hapR_is_pat_vs_mat).ln() / (Decimal('2').ln())

    log_Odds_Ratio = round(float(log_Odds_Ratio), 3)

    del lh_hapL_vs_hapR_is_pat_vs_mat


    if log_Odds_Ratio >= float(lods_cut_off):
        haplotype_dad = haplotype_left
        haplotype_mom = haplotype_right

    elif log_Odds_Ratio <= -float(lods_cut_off):
        haplotype_dad = haplotype_right
        haplotype_mom = haplotype_left

    else:  # Cannot assign the haplotype, so sub the allele with Ns
        haplotype_mom = ['N']*len(haplotype_left)
        haplotype_dad = ['N']*len(haplotype_left)


    #### Write the final output for the given PI
    # pull the required part of original dataframe to write it to final output
    my_soif1_df = data_by_pi[[
        'CHROM', 'POS', 'REF', 'all-alleles', soif1 + ':PI', soif1 + ':PG_al']].reset_index(drop=True)


    ## Output in long format
    hap_segregated_df_long = pd.DataFrame(
        OrderedDict((('log2odds', log_Odds_Ratio), (dad_id, haplotype_dad), (mom_id, haplotype_mom))))

    # merge this segregated df to main df for f1soi
    updated_df_long = pd.concat([my_soif1_df, hap_segregated_df_long], axis=1)


    ## Output in wide format
       # transform the data structure to appropriate wide format
    contig_w = my_soif1_df['CHROM'].tolist()[0]
    pos_range = '-'.join([str(my_soif1_df['POS'].min()), str(my_soif1_df['POS'].max())])
    soif1_PI = my_soif1_df[soif1 + ':PI'].tolist()[0]

    updated_df_wide = pd.DataFrame(
        [[contig_w, pos_range, soif1_PI , '-'.join(haplotype_left), '-'.join(haplotype_right), log_Odds_Ratio,
          '-'.join(haplotype_dad), '-'.join(haplotype_mom)]],
        columns=['CHROM', 'POS_Range', soif1 +':PI', 'hap_left', 'hap_right', 'log2odds', dad_id, mom_id])


    ## ** For future: update the haplotype statistics (at this position) if desired.
    ## The stats can have following headers
    # contig	ms02g_PI	phased_PI	total_haplotypes	phased_haplotypes	\
    # total_Vars	phased_Vars	Ref_in_Mat	Ref_in_Pat

    return updated_df_long, updated_df_wide




'''function for name sorting while reading file.
   - This function helps to read the file in alpha-numerical order when multiprocessing. '''

numbers = re.compile(r'(\d+)')

def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])

    return parts



''' function to return the transitions probabilities from transition counts'''
def compute_probs(pX_Y, pX, modeis):
    # returns emission or transition probs from "X" to "Y"
    #** variable "modeis" not used now. Use it in future.

    return Decimal(pX_Y / pX)
    #return round(pX_Y / pX, 7)


''' function to find appropriate sample names in each parental background '''
def find_samples(samples, header_):

    header_ = header_.rstrip('\n').split('\t')

    if 'pre:' in samples:
        samples = samples.lstrip('pre:').split(',')
        sample_list = []
        for pref in samples:
            for names in header_:
                if names.startswith(pref):
                    sample_list.append(names.split(':')[0])
        sample_list = list(set(sample_list))
        sample_list = [((x + ':PI'), (x + ':PG_al')) for x in sample_list]

    else:
        sample_list = samples.split(',')
        sample_list = [((x + ':PI'), (x + ':PG_al')) for x in sample_list]

    return sample_list



''' function to control the cumulation of the likelihoods.
    - We can either maxSum or maxProduct
    - cumulation is done two times:
        - once within markov chain
        - and futher cumulation of estimates from forward and reverse chain.
'''
def prod(iterable):
    value = 1
    for nth in iterable:
        value *= nth
    return value


def cumulate_likelihoods(items):
    if maxed_as == '+':
        cuml_is = sum(items)

    elif maxed_as == '*':
        cuml_is = prod(items)

    return cuml_is



''' to monitor memory usage. '''
def current_mem_usage():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.



if __name__ == "__main__":
    main();


