__author__ = 'Lindo Nkambule'

import hail as hl


def pca_filter_mt(
        in_mt: hl.MatrixTable,
        maf: float = 0.05,
        hwe: float = 1e-3,
        call_rate: float = 0.98,
        ld_cor: float = 0.2,
        ld_window: int = 250000):

    print("\nInitial number of SNPs before filtering: {}".format(in_mt.count_rows()))
    mt = hl.variant_qc(in_mt)
    print(f'\nFiltering out variants with MAF < {maf}')
    mt_filt = mt.annotate_rows(maf=hl.min(mt.variant_qc.AF))
    mt_filt = mt_filt.filter_rows(mt_filt.maf > maf)

    print(f'\nFiltering out variants with HWE < {hwe:1e}')
    mt_filt = mt_filt.filter_rows(mt_filt.variant_qc.p_value_hwe > hwe)

    print(f'\nFiltering out variants with Call Rate < {call_rate}')
    mt_filt = mt_filt.filter_rows(mt_filt.variant_qc.call_rate >= call_rate)

    # no strand ambiguity
    print('\nFiltering out strand ambigous variants')
    mt_filt = mt_filt.filter_rows(~hl.is_strand_ambiguous(mt_filt.alleles[0], mt_filt.alleles[1]))

    # MHC chr6:25-35Mb
    # chr8.inversion chr8:7-13Mb
    print('\nFiltering out variants in MHC [chr6:25M-35M] and chromosome 8 inversions [chr8:7M-13M]')
    intervals = ['chr6:25M-35M', 'chr8:7M-13M']
    mt_filt = hl.filter_intervals(mt_filt, [hl.parse_locus_interval(x, reference_genome='GRCh38') for x in intervals],
                                  keep=False)

    # This step is expensive (on local machine)
    print(f'\nLD pruning using correlation threshold of {ld_cor} and window size of {ld_window}')
    mt_ld_prune = hl.ld_prune(mt_filt.GT, r2=ld_cor, bp_window_size=ld_window)
    mt_ld_pruned = mt_filt.filter_rows(hl.is_defined(mt_ld_prune[mt_filt.row_key]))
    print("\nNumber of SNPs after filtering: {}".format(mt_ld_pruned.count_rows()))

    return mt_ld_pruned


def relatedness_check(
        in_mt: hl.MatrixTable = None,
        method: str = 'pc_relate',
        outdir: str = None,
        kin_estimate: float = 0.98):

    global mt, samples_to_remove

    in_mt = hl.variant_qc(in_mt)
    in_mt = hl.sample_qc(in_mt)

    # _localize=False means don't put this in Python, keep it as a Hail expr
    call_rate_dict = in_mt.aggregate_cols(hl.dict(hl.agg.collect((in_mt.s, in_mt.sample_qc.call_rate))), _localize=False)

    if method == 'pc_relate':
        print("\nUsing PC-Relate for relatedness checks")
        relatedness_ht = hl.pc_relate(in_mt.GT, 0.01, k=10, min_kinship=0.1, statistics='kin')
        samples_to_remove_ht = relatedness_ht.filter(relatedness_ht.kin > kin_estimate)

        # get call rates for both samples so we remove the one with lower call rate between the two
        samples_to_remove = samples_to_remove_ht.annotate(cr_s1=call_rate_dict[samples_to_remove_ht.i.s],
                                                          cr_s2=call_rate_dict[samples_to_remove_ht.j.s])

        samples_list = samples_to_remove.annotate(sample_to_remove=hl.cond(
            samples_to_remove.cr_s1 <= samples_to_remove.cr_s2, samples_to_remove.i, samples_to_remove.j))

    elif method == 'ibd':
        print("\nUsing PLINK-style identity by descent for relatedness checks")
        in_mt = in_mt.annotate_rows(maf=hl.min(in_mt.variant_qc.AF))
        relatedness_ht = hl.identity_by_descent(in_mt, maf=in_mt['maf'])  # this returns a Hail Table with the sample pairs
        samples_to_remove_ht = relatedness_ht.filter(relatedness_ht.ibd.PI_HAT > kin_estimate)

        # get call rates for both samples so we remove the one with lower call rate between the two
        samples_to_remove = samples_to_remove_ht.annotate(cr_s1=call_rate_dict[samples_to_remove_ht.i],
                                                          cr_s2=call_rate_dict[samples_to_remove_ht.j])

        samples_list = samples_to_remove.annotate(sample_to_remove=hl.cond(
            samples_to_remove.cr_s1 <= samples_to_remove.cr_s2, samples_to_remove.i, samples_to_remove.j))

    else:
        print("\nUsing KING for relatedness checks")
        if kin_estimate > 0.5:
            raise Exception("\nThe maximum kinship coefficient is for KING 0.5")
        relatedness_mt = hl.king(in_mt.GT)
        filtered_relatedness_mt = relatedness_mt.filter_entries((relatedness_mt.s_1 != relatedness_mt.s) &
                                                                (relatedness_mt.phi >= kin_estimate), keep=True)
        samples_to_remove_ht = filtered_relatedness_mt.entries()
        samples_to_remove = samples_to_remove_ht.annotate(cr_s1=call_rate_dict[samples_to_remove_ht.s_1],
                                                          cr_s2=call_rate_dict[samples_to_remove_ht.s])

        samples_list = samples_to_remove.annotate(sample_to_remove=hl.cond(
            samples_to_remove.cr_s1 <= samples_to_remove.cr_s2, samples_to_remove.s_1, samples_to_remove.s))

    samples = samples_list.sample_to_remove.collect()

    if len(samples) > 0:
        in_mt = in_mt.filter_cols(hl.literal(samples).contains(in_mt['s']), keep=False)
        print("\nNumber of samples that fail relatedness checks: {}".format(len(samples)))
        with open(outdir + 'relatedness_removed_samples.tsv', 'w') as f:
            for sample in samples:
                f.write(sample + "\n")

    else:
        print("\nNo samples failed the relatedness check")

    return in_mt
