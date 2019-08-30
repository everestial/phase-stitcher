import os
import shutil

def args_to_values(args):
    input_file = args.input  # the input haplotype file
    # input_file = 'allele_table_for_phase_extender.txt'
    print('  - using haplotype file "%s" ' % (input_file))

    soif1 = args.f1Sample
    print('  - F1-hybrid of interest: "%s" ' % (soif1))


    # only read the header of the input haplotype file and ...
    # .. find the required samples (for maternal and paternal background)
    header_ = open(input_file, 'r').readline()
    #print(header_)

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
    
    return input_file, soif1, soimom, soidad, dad_id, mom_id, outputdir, chr_list, nt, lods_cut_off, max_is, maxed_as, hapstats

    


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