__author__ = 'Lindo Nkambule'

import hail as hl


def pca_filter_mt(
        in_mt: hl.MatrixTable,
        maf: float = 0.05,
        hwe: float = 1e-3,
        call_rate: float = 0.98):

    print("\nInitial number of SNPs before filtering: {}".format(in_mt.count_rows()))
    mt = hl.variant_qc(in_mt)
    mt_filt = mt.annotate_rows(maf=hl.min(mt.variant_qc.AF))
    mt_filt = mt_filt.filter_rows(mt_filt.maf > maf)
    # mt_filt = mt.filter_rows((mt.variant_qc.AF[0] > 0.001) & (mt.variant_qc.AF[0] < 0.999))
    print("\nNumber of SNPs after MAF filtering: {}".format(mt_filt.count_rows()))

    mt_filt = mt_filt.filter_rows(mt_filt.variant_qc.p_value_hwe > hwe)
    print("\nNumber of SNPs after HWE filtering: {}".format(mt_filt.count_rows()))

    mt_filt = mt_filt.filter_rows(mt_filt.variant_qc.call_rate >= call_rate)
    print("\nNumber of SNPs after Call Rate filtering: {}".format(mt_filt.count_rows()))

    # no strand ambiguity
    mt_filt = mt_filt.filter_rows(~hl.is_strand_ambiguous(mt_filt.alleles[0], mt_filt.alleles[1]))
    print("\nNumber of SNPs after strand ambiguity filtering: {}".format(mt_filt.count_rows()))

    # MHC chr6:25-35Mb
    # chr8.inversion chr8:7-13Mb
    intervals = ['chr6:25M-35M', 'chr8:7M-13M']
    mt_filt = hl.filter_intervals(mt_filt, [hl.parse_locus_interval(x, reference_genome='GRCh38') for x in intervals],
                                  keep=False)
    print("\nNumber of SNPs after MHC and chr8 inversions filtering: {}".format(mt_filt.count_rows()))

    # This step is expensive (on local machine)
    mt_ld_prune = hl.ld_prune(mt_filt.GT, r2=0.2, bp_window_size=250000)
    mt_ld_pruned = mt.filter_rows(hl.is_defined(mt_ld_prune[mt.row_key]))
    print("\nNumber of SNPs after LD pruning: {}".format(mt_ld_pruned.count_rows()))

    return mt_ld_pruned


def relatedness_check(
        in_mt: hl.MatrixTable = None,
        method: str = 'ibd',
        outdir: str = None,
        ibd: float = 0.98):

    global mt

    in_mt = hl.variant_qc(in_mt)
    in_mt = hl.sample_qc(in_mt)

    if method == 'pc_relate':
        print("\nUsing PC-Relate for relatedness checks")

    elif method == 'ibd':
        print("\nUsing PLINK-style identity by descent for relatedness checks")
        mt = in_mt.annotate_rows(maf=hl.min(in_mt.variant_qc.AF))
        relatedness_ht = hl.identity_by_descent(mt, maf=mt['maf'])  # this returns a Hail Table with the sample pairs
        samples_to_remove_ht = relatedness_ht.filter(relatedness_ht.ibd.PI_HAT > ibd)

        # _localize=False means don't put this in Python, keep it as a Hail expr
        call_rate_dict = mt.aggregate_cols(hl.dict(hl.agg.collect((mt.s, mt.sample_qc.call_rate))), _localize=False)
        # get call rates for both samples so we remove the one with lower call rate between the two
        samples_to_remove = samples_to_remove_ht.annotate(cr_s1=call_rate_dict[samples_to_remove_ht.i],
                                                          cr_s2=call_rate_dict[samples_to_remove_ht.j])
        samples_list = samples_to_remove.annotate(sample_to_remove=hl.cond(
            samples_to_remove.cr_s1 >= samples_to_remove.cr_s2, samples_to_remove.i, samples_to_remove.j))
        samples = samples_list.sample_to_remove.collect()

        if len(samples) > 0:
            mt = mt.filter_cols(hl.literal(samples).contains(mt['s']), keep=False)
            print("Number of samples removed after relatedness checks: {}".format(len(samples)))
            with open(outdir + 'relatedness_removed_samples.tsv', 'w') as f:
                for sample in samples:
                    f.write(sample + "\n")

        else:
            print("No samples failed the relatedness check")

    else:
        print("\nUsing KING for relatedness checks")

    return mt
