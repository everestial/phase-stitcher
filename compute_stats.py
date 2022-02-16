import numpy as np
import pandas as pd


def compute_stats(outputdir, soif1, path_wide, lods_cut_off):
    result_merged_width = pd.read_csv(path_wide, sep="\t")
    # result_merged_width['log2odds'].round(3)

    grouped_df = result_merged_width.groupby("CHROM")
    output_df = grouped_df.apply(compute_func, soif1=soif1, lods_cut_off=lods_cut_off)
    result = output_df.replace("", ".").reset_index()
    result.columns = [
        "CHROM",
        "totalPhasedBlocks",
        "totalUnPhasedBlocks",
        "totalVarsPhasedBlocks",
        "totalVarsUnPhasedBlocks",
        "phasedBlocksPI",
        "unphasedBlocksPI",
        "numVarsInPhasedBlock",
        "numVarsInUnPhasedBlock",
        "log2oddsInPhasedBlock",
        "log2oddsInUnPhasedBlock",
        "totalNumOfBlock",
        "totalNumOfVars",
    ]

    file_path = outputdir + "/" + soif1 + "_haplotype_stats.txt"
    result.to_csv(file_path, sep="\t", index=False, float_format="%g")


def compute_func(df, soif1, lods_cut_off):
    log2odds_dist = np.round(df["log2odds"].values, 3)

    # name of all haplotype blocks (i.e PI) in that CHROM
    hapblocks = df[soif1 + ":PI"].values

    # among all the haplotype blocks find the the ..
    # .. indexes of the ones that are phased vs. unphased in final haplotype
    idx_ofphasedblocks = [
        ith
        for ith in range(len(log2odds_dist))
        if abs(log2odds_dist[ith]) >= lods_cut_off
    ]
    idx_ofunphasedblocks = [
        ith
        for ith in range(len(log2odds_dist))
        if abs(log2odds_dist[ith]) < lods_cut_off
    ]

    # based on that index find the phased blocks vs. unphased blocks
    phasedblocks = [hapblocks[zth] for zth in idx_ofphasedblocks]
    unphasedblocks = [hapblocks[zth] for zth in idx_ofunphasedblocks]

    # compute number of variants in all vs. phased vs. unphased blocks
    num_vars_in_allblocks = [len(ith.split("-")) for ith in df["hap_left"].tolist()]
    num_vars_phased_blocks = [num_vars_in_allblocks[zth] for zth in idx_ofphasedblocks]
    num_vars_unphased_blocks = [
        num_vars_in_allblocks[zth] for zth in idx_ofunphasedblocks
    ]

    # count total number of phased vs. unphased blocks (by PI)
    # and the total number of variants in each category
    total_phased_blocks = len(phasedblocks)
    total_unphased_blocks = len(unphasedblocks)

    total_vars_in_phasedblocks = sum(num_vars_phased_blocks)
    total_vars_in_unphasedblocks = sum(num_vars_unphased_blocks)

    # convert values to appropriate structure and string type before writing
    phasedblocks = ",".join(str(x) for x in phasedblocks)
    unphasedblocks = ",".join(str(x) for x in unphasedblocks)
    num_vars_phased_blocks = ",".join([str(x) for x in num_vars_phased_blocks])
    num_vars_unphased_blocks = ",".join([str(x) for x in num_vars_unphased_blocks])

    # find the log2odds of the phased vs. unphased blocks
    log2_phased_blocks = [log2odds_dist[zth] for zth in idx_ofphasedblocks]
    log2_unphased_blocks = [log2odds_dist[zth] for zth in idx_ofunphasedblocks]

    log2_phased_blocks = ",".join([str(x) for x in log2_phased_blocks])
    log2_unphased_blocks = ",".join([str(x) for x in log2_unphased_blocks])

    # log2_phased_blocks = ','.join([format(x, "10.3f") for x in log2_phased_blocks])
    # log2_unphased_blocks = ','.join([format(x, "10.3f") for x in log2_unphased_blocks])

    # compute total number of blocks and total number of variants
    count_blocks = len(log2odds_dist)
    count_vars = sum(num_vars_in_allblocks)
    return pd.Series(
        [
            total_phased_blocks,
            total_unphased_blocks,
            total_vars_in_phasedblocks,
            total_vars_in_unphasedblocks,
            phasedblocks,
            unphasedblocks,
            num_vars_phased_blocks,
            num_vars_unphased_blocks,
            log2_phased_blocks,
            log2_unphased_blocks,
            count_blocks,
            count_vars,
        ]
    )


# def compute_stats(outputdir, soif1, result_merged_width,lods_cut_off):
#         with open(outputdir + '/' + soif1 + '_haplotype_stats.txt', 'w+') as out_stats:
#             stats = result_merged_width.groupby(['CHROM'], sort = False)
#             out_stats.write('\t'.join(['CHROM', 'totalPhasedBlocks', 'totalUnPhasedBlocks', 'totalVarsPhasedBlocks',
#                                        'totalVarsUnPhasedBlocks', 'phasedBlocksPI', 'unphasedBlocksPI', 'numVarsInPhasedBlock',
#                             'numVarsInUnPhasedBlock', 'log2oddsInPhasedBlock',
#                             'log2oddsInUnPhasedBlock',	'totalNumOfBlock',	'totalNumOfVars']) + '\n')

#             # ** future: additional column to add: totalPhasedBlocks and totalunphasedBlocks,
#             # totalPhasedVars, totalunPhasedVars,
#             # Rename: phaseBlock to phasedBlocksPI

#             for xath, valth in stats:
#                 print('Computing descriptive statistics of the final haplotype for contig "%s". ' %str(xath))

#                 #distribution of the log2odds values in that CHROM
#                 log2odds_dist = valth['log2odds'].tolist()

#                 # name of all haplotype blocks (i.e PI) in that CHROM
#                 hapblocks = valth[soif1 + ':PI'].tolist()

#                 # among all the haplotype blocks find the the ..
#                 # .. indexes of the ones that are phased vs. unphased in final haplotype
#                 idx_ofphasedblocks = [ith for ith in range(len(log2odds_dist))
#                                       if abs(log2odds_dist[ith]) >= lods_cut_off]
#                 idx_ofunphasedblocks = [ith for ith in range(len(log2odds_dist))
#                                         if abs(log2odds_dist[ith]) < lods_cut_off]

#                 # based on that index find the phased blocks vs. unphased blocks
#                 phasedblocks = [hapblocks[zth] for zth in idx_ofphasedblocks]
#                 unphasedblocks = [hapblocks[zth] for zth in idx_ofunphasedblocks]

#                 # compute number of variants in all vs. phased vs. unphased blocks
#                 num_vars_in_allblocks = [len(ith.split('-')) for ith in valth['hap_left'].tolist()]
#                 num_vars_phased_blocks = [num_vars_in_allblocks[zth] for zth in idx_ofphasedblocks]
#                 num_vars_unphased_blocks = [num_vars_in_allblocks[zth] for zth in idx_ofunphasedblocks]


#                 # count total number of phased vs. unphased blocks (by PI)
#                 # and the total number of variants in each category
#                 total_phased_blocks = str(len(phasedblocks))
#                 total_unphased_blocks = str(len(unphasedblocks))

#                 total_vars_in_phasedblocks = str(sum(num_vars_phased_blocks))
#                 total_vars_in_unphasedblocks = str(sum(num_vars_unphased_blocks))


#                 # convert values to appropriate structure and string type before writing
#                 phasedblocks = ','.join(str(x) for x in phasedblocks)
#                 unphasedblocks = ','.join(str(x) for x in unphasedblocks)
#                 num_vars_phased_blocks = ','.join([str(x) for x in num_vars_phased_blocks])
#                 num_vars_unphased_blocks = ','.join([str(x) for x in num_vars_unphased_blocks])

#                 # find the log2odds of the phased vs. unphased blocks
#                 log2_phased_blocks = [log2odds_dist[zth] for zth in idx_ofphasedblocks]
#                 log2_unphased_blocks = [log2odds_dist[zth] for zth in idx_ofunphasedblocks]

#                 log2_phased_blocks = ','.join([str(x) for x in log2_phased_blocks])
#                 log2_unphased_blocks = ','.join([str(x) for x in log2_unphased_blocks])

#                 # compute total number of blocks and total number of variants
#                 count_blocks = str(len(log2odds_dist))
#                 count_vars = str(sum(num_vars_in_allblocks))


#                 data_to_write = [str(xath), total_phased_blocks, total_unphased_blocks, total_vars_in_phasedblocks,
#                                  total_vars_in_unphasedblocks, phasedblocks, unphasedblocks, num_vars_phased_blocks,
#                                  num_vars_unphased_blocks, log2_phased_blocks, log2_unphased_blocks,
#                                  count_blocks, count_vars]
#                 data_to_write = ['.' if x == '' else x for x in data_to_write]

#                 out_stats.write('\t'.join(data_to_write))
#                 out_stats.write('\n')
