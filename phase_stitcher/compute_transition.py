from collections import OrderedDict
from decimal import Decimal
from itertools import product, islice

import pandas as pd


"""Step 04 - A : Estimate the likelihood of a haplotype being mom vs. dad. This function needs inputs:
    - RBphased data (as my df dict), left vs. right haplotype for soif1, parents, orientation (fwd vs. rev). """


def compute_transition(
    my_df_dict, haplotype_left, haplotype_right, soimeta, parent, orientation, maxed_as
):

    soif1 = soimeta.soif1
    parent = soimeta.soidad if parent == "dad" else soimeta.soimom

    """Create a variable that stores likelihood values, that ..
       left vs. right haplotype belongs to a particular parent.
       This depends upon which parent was passed in this function"""

    """ We assign the prior likelihood that left and right haplotype transition equally belong to that parent.
        This likelihood estimates then cumulates (maxSum or maxPd) as we get more evidence from emissions and
        transitions frequency. """
    ## likelihood = prob(for any random emission) x prob(of any random transition) = 0.5 * 0.0625
    ## **Note: different level of value (prior) doesn't affect maxPd but only maxSum.
    likelihood_hap_left = Decimal(0.25 * 0.0625)  # lambda x: x
    likelihood_hap_right = Decimal(0.25 * 0.0625)  # lambda x: x

    """Step 04 - A (i): Run 1st order markov chain on one block of unique "PI" data. """
    for v1, v2 in zip(
        orientation(list(enumerate(my_df_dict[soif1 + ":PG_al"]))),
        islice(orientation(list(enumerate(my_df_dict[soif1 + ":PG_al"]))), 1, None),
    ):

        n1 = v1[0]
        n2 = v2[0]

        """Skip the markov chain if indels is present in either "n1" or "n2" level. """
        if len(v1[1]) > 3 or len(v2[1]) > 3:
            continue

        """Skip the markov chain if symbolic allele (e.g "*") is present in either "n1" or "n2" level. """
        if "*" in v1[1] or "*" in v2[1]:
            continue

        """ Step: Estimate emission counts and emission probabilities at "n1" level. """
        nucleotide_prob_dict, nucleotide_count_dict = compute_emission_prob(
            parent, my_df_dict, n1
        )

        """ Step : Estimate transition counts and probabilities from "n1" to "n2". """

        transition_prob_dict = compute_tr_prob(
            parent, my_df_dict, n1, n2, nucleotide_count_dict
        )

        """ find observed transition configuration for left and right
        haplotype from "n-1" to "n" position"""

        hap_transition_left = (haplotype_left[n1], haplotype_left[n2])
        hap_transition_right = (haplotype_right[n1], haplotype_right[n2])

        """ Now, cumulate the likelihood (as cumulation of the (emission prob * transition prob)). """
        ## cumulation of the probability can be either maxSum or maxProduct.

        likelihood_hap_left = cumulate_likelihoods(
            [
                likelihood_hap_left,
                nucleotide_prob_dict[haplotype_left[n1]]
                * transition_prob_dict[hap_transition_left],
            ],
            maxed_as,
        )

        likelihood_hap_right = cumulate_likelihoods(
            [
                likelihood_hap_right,
                nucleotide_prob_dict[haplotype_right[n1]]
                * transition_prob_dict[hap_transition_right],
            ],
            maxed_as,
        )

        # print('likelihoods left and right')
        # print(likelihood_hap_left, likelihood_hap_right)

    return likelihood_hap_left, likelihood_hap_right


def compute_lods(
    likelihood_hap_left_dad,
    likelihood_hap_right_dad,
    likelihood_hap_left_mom,
    likelihood_hap_right_mom,
    haplotype_left,
    haplotype_right,
    data_by_pi,
    soimeta,
    lods_cut_off,
):

    """Segregate the haplotype and write it to a file."""
    lh_hapL_vs_hapR_is_pat_vs_mat = Decimal(
        (likelihood_hap_left_dad / likelihood_hap_left_mom)
        / (likelihood_hap_right_dad / likelihood_hap_right_mom)
    )

    log_Odds_Ratio = Decimal(lh_hapL_vs_hapR_is_pat_vs_mat).ln() / (Decimal("2").ln())

    log_Odds_Ratio = round(float(log_Odds_Ratio), 3)

    del lh_hapL_vs_hapR_is_pat_vs_mat

    if log_Odds_Ratio >= float(lods_cut_off):
        haplotype_dad = haplotype_left
        haplotype_mom = haplotype_right

    elif log_Odds_Ratio <= -float(lods_cut_off):
        haplotype_dad = haplotype_right
        haplotype_mom = haplotype_left

    else:  # Cannot assign the haplotype, so sub the allele with Ns
        haplotype_mom = ["N"] * len(haplotype_left)
        haplotype_dad = ["N"] * len(haplotype_left)

    #### Write the final output for the given PI
    # pull the required part of original dataframe to write it to final output
    my_soif1_df = data_by_pi[
        [
            "CHROM",
            "POS",
            "REF",
            "all-alleles",
            soimeta.soif1 + ":PI",
            soimeta.soif1 + ":PG_al",
        ]
    ].reset_index(drop=True)

    ## Output in long format
    hap_segregated_df_long = pd.DataFrame(
        OrderedDict(
            (
                ("log2odds", log_Odds_Ratio),
                (soimeta.dad_id, haplotype_dad),
                (soimeta.mom_id, haplotype_mom),
            )
        )
    )

    # merge this segregated df to main df for f1soi
    updated_df_long = pd.concat([my_soif1_df, hap_segregated_df_long], axis=1)

    ## Output in wide format
    # transform the data structure to appropriate wide format
    contig_w = my_soif1_df["CHROM"].tolist()[0]
    pos_range = "-".join([str(my_soif1_df["POS"].min()), str(my_soif1_df["POS"].max())])
    soif1_PI = my_soif1_df[soimeta.soif1 + ":PI"].tolist()[0]

    updated_df_wide = pd.DataFrame(
        [
            [
                contig_w,
                pos_range,
                soif1_PI,
                "-".join(haplotype_left),
                "-".join(haplotype_right),
                log_Odds_Ratio,
                "-".join(haplotype_dad),
                "-".join(haplotype_mom),
            ]
        ],
        columns=[
            "CHROM",
            "POS_Range",
            soimeta.soif1 + ":PI",
            "hap_left",
            "hap_right",
            "log2odds",
            soimeta.dad_id,
            soimeta.mom_id,
        ],
    )

    ## ** For future: update the haplotype statistics (at this position) if desired.
    ## The stats can have following headers
    # contig	ms02g_PI	phased_PI	total_haplotypes	phased_haplotypes	\
    # total_Vars	phased_Vars	Ref_in_Mat	Ref_in_Pat

    return updated_df_long, updated_df_wide


""" function to return the transitions probabilities from transition counts"""


def compute_probs(pX_Y, pX, modeis):
    # returns emission or transition probs from "X" to "Y"
    # ** variable "modeis" not used now. Use it in future.

    return Decimal(pX_Y / pX)
    # return round(pX_Y / pX, 7)


""" function to control the cumulation of the likelihoods.
    - We can either maxSum or maxProduct
    - cumulation is done two times:
        - once within markov chain
        - and futher cumulation of estimates from forward and reverse chain.
"""


def prod(iterable):
    value = 1
    for nth in iterable:
        value *= nth
    return value


def cumulate_likelihoods(items, maxed_as):
    if maxed_as == "+":
        cuml_is = sum(items)

    elif maxed_as == "*":
        cuml_is = prod(items)

    return cuml_is


def compute_emission_prob(parent, my_df_dict, n1):
    nucleotide_count_dict = {"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}

    nucleotide_prob_dict = nucleotide_count_dict.copy()

    # using for loop over the given sample to compute emission counts
    for (x, y) in parent:
        for nucleotide in "ATGC":
            # count only after splitting. else 'AT|C' will have one count of 'A'
            nucleotide_count_dict[nucleotide] += (
                my_df_dict[y][n1].split("|").count(nucleotide)
            )

    # now, compute emission probabilities
    total_nucl = sum(list(nucleotide_count_dict.values()))

    for element in nucleotide_prob_dict:
        nucleotide_prob_dict[element] = compute_probs(
            nucleotide_count_dict[element], total_nucl, modeis="emission"
        )

    return nucleotide_prob_dict, nucleotide_count_dict


def compute_tr_prob(parent, my_df_dict, n1, n2, nucleotide_count_dict):
    # empty dict to store transition data
    # here we use 1/16 as a prior probability distribution of each transition under random assumption
    # this is supplied as pseudo estimate and serves as minimal observation (probability) when there is no data ..
    # .. for the given transition
    transition_count_dict = {
        ("A", "A"): 1 / 16,
        ("A", "T"): 1 / 16,
        ("A", "G"): 1 / 16,
        ("A", "C"): 1 / 16,
        ("T", "A"): 1 / 16,
        ("T", "T"): 1 / 16,
        ("T", "G"): 1 / 16,
        ("T", "C"): 1 / 16,
        ("G", "A"): 1 / 16,
        ("G", "T"): 1 / 16,
        ("G", "G"): 1 / 16,
        ("G", "C"): 1 / 16,
        ("C", "A"): 1 / 16,
        ("C", "T"): 1 / 16,
        ("C", "G"): 1 / 16,
        ("C", "C"): 1 / 16,
    }

    # shallow copy should be fine  # **for future - do deep copy if need be
    transition_prob_dict = transition_count_dict.copy()

    ## Estimate transition counts
    for x, y in parent:
        nucl_n1 = (my_df_dict[y][n1]).split("|")  # nucleotides at n1 level
        nucl_n2 = (my_df_dict[y][n2]).split("|")  # nucleotides at n2 level

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
                transition_count_dict[(from_, to)] += (
                    prod_nuclb1b2.count((from_, to)) / 2
                )

    ## convert transition counts to probabilities
    for (from_, to) in transition_prob_dict:
        transition_prob_dict[(from_, to)] = compute_probs(
            transition_count_dict[(from_, to)],
            nucleotide_count_dict[from_],
            modeis="transition",
        )

    return transition_prob_dict
