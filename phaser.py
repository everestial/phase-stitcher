import itertools
import os
import shutil
import time

from pathlib import Path
from functools import partial
from multiprocessing import Pool

import pandas as pd
import dask.dataframe as dd
from compute_stats import compute_stats
from compute_transition import compute_lods, compute_transition, cumulate_likelihoods
from utils import current_mem_usage

# from dowser.utils import launch_memory_usage_server
# launch_memory_usage_server()


def phase_stich(
    input_file, soimeta, outputdir, chr_list, nt, lods_cut_off, maxed_as, hapstats
):

    # number of pool to run at once; default at 1

    """Step 01: Read the input file and and prepare two files as output"""

    # load data as dask dataframe
    data = dd.read_csv(input_file, sep="\t")

    # write column header in missing datafile
    with open(
        outputdir + "/" + soimeta.soif1 + "_missing_data.txt", "w+"
    ) as missing_data:
        missing_data.write("\t".join(list(data)) + "\n")

    my_df_by_contig = data.groupby("CHROM")
    # df_list = (my_df_by_contig.get_group(x).compute() for x in good_data['CHROM'].unique())

    # get list of chrom values and sort them
    sorted_ch_list = sorted([x for x in data["CHROM"].unique()])
    # missing_df= my_df_by_contig
    df_list = (my_df_by_contig.get_group(chr) for chr in sorted_ch_list)

    """This imap multiprocess take the pipeline to Step 03. """

    time_all = time.time()

    df_list_c = (df.compute() for df in df_list)

    partial_func = partial(
        process_by_contig,
        outputdir=outputdir,
        soimeta=soimeta,
        maxed_as=maxed_as,
        lods_cut_off=lods_cut_off,
    )
    with Pool(processes=nt) as pool:
        results = pool.imap(partial_func, df_list_c)
        pool.close()
        pool.join()
        # pool.terminate()

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
    # long_exp = [result[0] for result in results]
    # wide_exp = [result[1] for result in results]
    # result_long= pd.concat(long_exp, axis = 0)
    # result_wide = pd.concat(wide_exp, axis = 0)
    # # results = list(results)
    # for result in results:
    #     result_long.concat(result[0], axis = 1)
    #     result_wide.concat(result[1], axis = 1)

    path_long = outputdir + "/" + soimeta.soif1 + "_haplotype_long.txt"
    path_wide = outputdir + "/" + soimeta.soif1 + "_haplotype_wide.txt"
    for long_df, wide_df in results:
        if not os.path.isfile(path_long):
            long_df.to_csv(path_long, header=True, index=False, sep="\t")
        else:
            long_df.to_csv(path_long, mode="a", header=False, index=False, sep="\t")

        if not os.path.isfile(path_wide):
            wide_df.to_csv(path_wide, header=True, index=False, sep="\t")
        else:
            wide_df.to_csv(path_wide, mode="a", header=False, index=False, sep="\t")

    # with open(path_long, 'a') as long_f, open(path_wide, 'a') as wide_f:
    # for long_df, wide_df in results:
    #     long_df.to_csv(path_long, mode='a', header=not long_f.tell(), index=False, sep='\t')

    #     wide_df.to_csv(path_wide, mode='a', header=not wide_f.tell(), index=False, sep='\t')

    # result_long = (result[0] for result in results)
    # result_wide = (result[1] for result in results)

    # ## merge all the data both length and widthwise
    # result_merged_length = pd.concat(result_long)
    # result_merged_width = pd.concat(result_wide)

    """ Step : merge the returned result (extended haplotype block by each contigs)
            #together and write it to an output file. """

    # write this updated dataframe to file.
    # ** for future: or write is directly after "imap'ing". The main problem is that they complete at different times
    # and complete result is not obtained. ** problem not known.
    # result_merged_length.to_csv(outputdir + '/' + soimeta.soif1 + '_haplotype_long.txt',
    #                     header=True, index=False, sep='\t')

    # result_merged_width.to_csv(outputdir + '/' + soimeta.soif1 + '_haplotype_wide.txt',
    #                     header=True, index=False, sep='\t')

    print(
        "Completed haplotype segregation and stitching for all the chromosomes.\n"
        "Time elapsed: '%f' sec. " % (time.time() - time_all)
    )
    print(" - Global maximum memory usage: %.2f (mb)" % current_mem_usage())
    print("Completed writing the dataframes .....")

    # ''' Prepare the descriptive statistics of the phased (stitched) haplotypes if desired: '''
    # ** future - move this code to a separate function.
    # result_merged_width = pd.read_csv(path_wide, sep = '\t')
    if hapstats == "yes":
        compute_stats(outputdir, soimeta.soif1, path_wide, lods_cut_off)

    elif hapstats == "no":
        print(
            "Proceeding without preparing descriptive statistics for final haplotype."
        )

    """Step 02 - A: Now, pipe the procedure to next function for multiprocessing (i.e Step 02 - B) """
    # multiproc(pool,  outputdir,soimeta, maxed_as, lods_cut_off, hapstats, df_list)

    """** workflow, between step (02-A and 06) is passed onto other functions. """

    """Step 06: Clean the directory that is storing splitted dataframe. """

    # remove the chunked data folder
    # shutil.rmtree('chunked_Data_' + soif1, ignore_errors=False, onerror=None)

    print("The End :)")


def dask_multiproc(
    contig_df, outputdir, soimeta, maxed_as, lods_cut_off, hapstats, my_df_by_contig
):
    pass


"""Step 02 - B : Pipe the data to multiprocessing environment. """


def multiproc(pool, outputdir, soimeta, maxed_as, lods_cut_off, hapstats, df_list):

    print()
    time_all = time.time()

    """ Step 02 - B: Start, multiprocessing/threading - process each contig separately. """
    # path = 'chunked_Data_' + soif1 + '/'  # opens the a temp folder
    # file_path = sorted([path + item for item in os.listdir(path)], key=numericalSort)

    """This imap multiprocess take the pipeline to Step 03. """
    # df_list = (my_df_by_contig.get_group(x) for x in my_df_by_contig.groups)
    # df_list = my_df_by_contig
    df_list_c = (df.compute() for df in df_list)

    partial_func = partial(
        process_by_contig,
        outputdir=outputdir,
        soimeta=soimeta,
        maxed_as=maxed_as,
        lods_cut_off=lods_cut_off,
    )
    results = pool.imap(partial_func, df_list_c)

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

    result_long = (result[0] for result in results)
    result_wide = (result[1] for result in results)

    ## merge all the data both length and widthwise
    result_merged_length = pd.concat(result_long)
    result_merged_width = pd.concat(result_wide)

    """ Step : merge the returned result (extended haplotype block by each contigs)
            #together and write it to an output file. """

    # write this updated dataframe to file.
    # ** for future: or write is directly after "imap'ing". The main problem is that they complete at different times
    # and complete result is not obtained. ** problem not known.
    result_merged_length.to_csv(
        outputdir + "/" + soimeta.soif1 + "_haplotype_long.txt",
        header=True,
        index=False,
        sep="\t",
    )

    result_merged_width.to_csv(
        outputdir + "/" + soimeta.soif1 + "_haplotype_wide.txt",
        header=True,
        index=False,
        sep="\t",
    )

    print(
        "Completed haplotype segregation and stitching for all the chromosomes.\n"
        "Time elapsed: '%f' sec. " % (time.time() - time_all)
    )
    print(" - Global maximum memory usage: %.2f (mb)" % current_mem_usage())
    print("Completed writing the dataframes .....")

    # ''' Prepare the descriptive statistics of the phased (stitched) haplotypes if desired: '''
    # ** future - move this code to a separate function.
    if hapstats == "yes":
        compute_stats(outputdir, soimeta.soif1, result_merged_width, lods_cut_off)

    elif hapstats == "no":
        print(
            "Proceeding without preparing descriptive statistics for final haplotype."
        )


"""Step 03: Process each contig separately with lazy mapping using "imap" multiprocessing. """


def process_by_contig(df_chr, outputdir, soimeta, maxed_as, lods_cut_off):

    # print()
    time_chr = time.time()

    # ''' Step 03 - A: Process data for each chromosome (contig) separately. '''
    # good_data_by_contig = open(file_path, 'r')
    # chr_ = good_data_by_contig.name.split(':')[-1]
    # this name identification might cause problem if "chr" has "_" in it's name i.e scaff_18

    print()
    # print('Starting markov transition for chromosome "%s" ' % (chr_))
    print("##########################")

    ## Load the data into pandas dataframe and again groupby unique "PI" values
    # contigs_group = pd.read_csv(StringIO(good_data_by_contig.read()), sep='\t')

    # groupby "PI index" and feed it to transition computation
    my_df_grouped = df_chr.groupby(soimeta.soif1 + ":PI", sort=False)

    # empty list to store the pandas dataframe
    haplotype_result_long = []
    haplotype_result_wide = []

    """Step 03 - B: process the data within each unique "PI" separately """
    for pi_index, data_by_pi in my_df_grouped:
        if pi_index == ".":
            pd.DataFrame.to_csv(
                data_by_pi,
                outputdir + "/" + soimeta.soif1 + "_missing_data.txt",
                header=False,
                index=False,
                mode="a",
                sep="\t",
            )

        else:
            ## pipe the data for computing transition and segregating the haplotypes
            my_df_dict = data_by_pi.to_dict(orient="list")

            ### Split the RBphase block into the left and right haplotype type
            haplotype_left = [
                x.split("|")[0] for x in my_df_dict[soimeta.soif1 + ":PG_al"]
            ]
            haplotype_right = [
                x.split("|")[1] for x in my_df_dict[soimeta.soif1 + ":PG_al"]
            ]

            """ Enter Step 04: """
            ## Pass the data to a function to compute likelihood of each configuration
            # - estimate likelihood that each haplotype (left vs. right) belong to mom vs. dad.
            # - estimate likelihood using both forward and reverse markov chain.
            # - ** for future (likelihood estimates can be multiprocessed)

            """Step 04 - A: Estimate likelihood on both forward and reverse chain."""

            """ likelihood of a haplotype (left vs. right) belonging to a parent (dad). """
            # on forward chain
            (
                likelihood_hap_left_dad_fwd,
                likelihood_hap_right_dad_fwd,
            ) = compute_transition(
                my_df_dict,
                haplotype_left,
                haplotype_right,
                soimeta,
                parent="dad",
                orientation=lambda x: x,
                maxed_as=maxed_as,
            )

            # on reverse chain
            (
                likelihood_hap_left_dad_rev,
                likelihood_hap_right_dad_rev,
            ) = compute_transition(
                my_df_dict,
                haplotype_left,
                haplotype_right,
                soimeta,
                parent="dad",
                orientation=reversed,
                maxed_as=maxed_as,
            )

            """ likelihood of a haplotype (left vs. right) belonging to (mom) or maternal background. """
            # on forward chain
            (
                likelihood_hap_left_mom_fwd,
                likelihood_hap_right_mom_fwd,
            ) = compute_transition(
                my_df_dict,
                haplotype_left,
                haplotype_right,
                soimeta,
                parent="mom",
                orientation=lambda x: x,
                maxed_as=maxed_as,
            )

            # on reverse chain
            (
                likelihood_hap_left_mom_rev,
                likelihood_hap_right_mom_rev,
            ) = compute_transition(
                my_df_dict,
                haplotype_left,
                haplotype_right,
                soimeta,
                parent="mom",
                orientation=reversed,
                maxed_as=maxed_as,
            )

            """Step 04 - B: Maximize the likelihood estimate (using maxSum or maxProduct). """
            likelihood_hap_left_dad = cumulate_likelihoods(
                [likelihood_hap_left_dad_fwd, likelihood_hap_left_dad_rev], maxed_as
            )

            likelihood_hap_right_dad = cumulate_likelihoods(
                [likelihood_hap_right_dad_fwd, likelihood_hap_right_dad_rev], maxed_as
            )

            likelihood_hap_left_mom = cumulate_likelihoods(
                [likelihood_hap_left_mom_fwd, likelihood_hap_left_mom_rev], maxed_as
            )

            likelihood_hap_right_mom = cumulate_likelihoods(
                [likelihood_hap_right_mom_fwd, likelihood_hap_right_mom_rev], maxed_as
            )

            """Step 04 - B: Segregate the haplotypes: Pass the likelihood estimates to another function
                            to compute log2Odds and segregate the haplotypes into mom vs. dad background. """
            updated_df_long, updated_df_wide = compute_lods(
                likelihood_hap_left_dad,
                likelihood_hap_right_dad,
                likelihood_hap_left_mom,
                likelihood_hap_right_mom,
                haplotype_left,
                haplotype_right,
                data_by_pi,
                soimeta,
                lods_cut_off,
            )

            ### Store the pandas dataframe as list
            haplotype_result_long.append(updated_df_long)
            haplotype_result_wide.append(updated_df_wide)

        if len(haplotype_result_long) == 0 or len(haplotype_result_wide) == 0:
            result_merged_long = pd.DataFrame()
            result_merged_wide = pd.DataFrame()
            result_merged_stats = pd.DataFrame()

        """Returns to Step 06: """
        ## merge the list of pandas dataframe into one dataframe
        # and then return it to the imap pool
        # ** - we can directly write the data to file using pd.csv (append) method
        try:
            result_merged_long = pd.concat(haplotype_result_long, axis=0)
            result_merged_wide = pd.concat(haplotype_result_wide, axis=0)

        except ValueError:  # Sometimes empty list is returned so, raising this exception
            print("No RBphased haplotypes available in this contig.")
            result_merged_long = pd.DataFrame()
            result_merged_wide = pd.DataFrame()

        # print(' - Phase-extension completed for contig "%s" in %.6f seconds' % (chr_, time.time() - time_chr))
        print(" - Worker maximum memory usage: %.2f (mb)" % (current_mem_usage()))
        print()

    return (result_merged_long, result_merged_wide)

    # ### Deprecated for now !!!
    # ''' Write the result: now, write the segregated haplotype to a output file.
    #     In both long and wide format.
    #     **Note: I decided to append each returned dataframe directly to output file, rather
    #     than passing it back to "imap" pool. I think this is more efficient (both time and memory). '''

    ## write results for long format
    # pd.DataFrame.to_csv(result_merged_long, outputdir + '/' + soif1 + '_haplotype_long.txt',
    #                    header=False, index=False, mode='a', sep='\t')

    ## write results for wide format
    # pd.DataFrame.to_csv(result_merged_wide, outputdir + '/' + soif1 + '_haplotype_wide.txt',
    #                    header=False, index=False, mode='a' ,sep='\t')
