"""
This function add the arguments and returns parsed values values.
"""


def get_args(parser):
    parser.add_argument(
        "--nt",
        help="number of process to run -> "
        "The maximum number of processes that can be run at once is the number "
        "of different chromosomes (contigs) in the input haplotype file.",
        default=1,
        required=False,
    )

    parser.add_argument(
        "--input",
        help="name of the input haplotype file -> "
        "This haplotype file should contain unique index represented by 'PI' and "
        "phased genotype represented by 'PG_al' for all the samples. ",
        required=True,
    )

    parser.add_argument(
        "--pat",
        help="Paternal sample or comma separated sample names that "
        "belong to Paternal background. "
        "Sample group may also be assigned using prefix. "
        "Options: 'paternal sample name', 'comma separated samples', 'pre:...'. "
        "Unique prefix (or comma separated prefixes) should begin with 'pre:'. ",
        required=True,
    )

    parser.add_argument(
        "--mat",
        help="Maternal sample or sample names (comma separated) that "
        "belong to maternal background. Sample group can also be "
        "assigned using unique prefix/es. "
        "Options: 'maternal sample name', 'comma separated samples', 'pre:...'. "
        "Unique prefix (or comma separated prefixes) should begin with 'pre:'. ",
        required=True,
    )

    parser.add_argument(
        "--f1Sample",
        help="Name of the F1-hybrid sample. Please type the name of only one F1 sample.",
        required=True,
    )

    parser.add_argument(
        "--outPatMatID",
        help="Prefix of the 'Paternal (dad)' and 'Maternal (mom)'genotype in the output file. "
        "This should be a maximum of three letter prefix separated by comma. "
        "Default: 'pat,mat'.",
        default="pat,mat",
        required=False,
    )

    parser.add_argument(
        "--output",
        default="",
        type=str,
        help="Name of the output directory. " "Default: f1SampleName + '_stitched' ",
        required=False,
    )

    parser.add_argument(
        "--lods",
        help="log(2) odds cutoff threshold required "
        "to assign maternal Vs. paternal haplotype segregation and stitching. ",
        default=5,
    )
    # note: log_odds_ratio 5 is 2^5 = 32 times likely

    parser.add_argument(
        "--culLH",
        help="Cumulative likelhood estimates -> "
        "The likelhoods for haplotype segregation can either be max-sum vs. max-product. "
        "Default: maxPd i.e max-product. "
        "Options: 'maxPd' or 'maxSum'. ",
        default="maxPd",
        required=False,
    )

    parser.add_argument(
        "--chr",
        required=False,
        default="",
        help="Restrict haplotype stitching to a specific chromosome.",
    )

    parser.add_argument(
        "--hapStats",
        help="Computes the descriptive statistics of final haplotype. "
        "Default: 'no'."
        "Option: 'yes', 'no' .",
        default="no",
        required=False,
    )

    args = parser.parse_args()
    return args
