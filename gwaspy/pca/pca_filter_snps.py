__author__ = 'Lindo Nkambule'

import hail as hl
import pandas as pd


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
        method: str = 'king',
        outdir: str = None,
        kin_estimate: float = 0.1,
        include_kinself: bool = False):

    if method == 'pc_relate':
        print("\nUsing PC-Relate for relatedness checks")
        # compute kinship statistic for every sample-pair
        if include_kinself:
            print("\nkinself will be included in exported tsv file")
        relatedness_ht = hl.pc_relate(in_mt.GT, 0.01, k=10, statistics='kin', include_self_kinship=include_kinself)
        
        print('exporting relatedness statistics to a tsv file')
        ht_export = relatedness_ht.key_by()
        ht_export = ht_export.select(ht_export.kin, i=ht_export.i.s, j=ht_export.j.s)
        ht_export.export(f'{outdir}relatedness_checks_pc_relate.tsv.bgz')

        print('getting related samples to be removed using maximal independent set')
        # only run maximal independent set step on sample-pairs with kinship above specified threshold

        # when include_kinself is True, not removing kinself will result in all samples failing relatedness because for
        # every kin between a sample with itself, the kin estimate will be ~0.5 in most cases (excluding inbreeding)
        if include_kinself:
            relatedness_ht = relatedness_ht.filter(relatedness_ht.i == relatedness_ht.j, keep=False)
        else:
            relatedness_ht = relatedness_ht
        pairs = relatedness_ht.filter(relatedness_ht['kin'] > kin_estimate)
        samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, False)
        samples = samples_to_remove.node.s.collect()

    elif method == 'ibd':
        print("\nUsing PLINK-style IBD for relatedness checks")
        in_mt = hl.variant_qc(in_mt)
        in_mt = in_mt.annotate_rows(maf=hl.min(in_mt.variant_qc.AF))
        relatedness_ht = hl.identity_by_descent(in_mt, maf=in_mt['maf'])
        
        print('exporting relatedness statistics to a tsv file')
        relatedness_ht.export(f'{outdir}relatedness_checks_ibd.tsv.bgz')

        print('getting related samples to be removed using maximal independent set')
        # only run maximal independent set step on sample-pairs with kinship above specified threshold
        pairs = relatedness_ht.filter(relatedness_ht['ibd.PI_HAT'] > kin_estimate)
        samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, False)
        samples = samples_to_remove.node.collect()

    else:
        print("\nUsing KING for relatedness checks")
        if kin_estimate > 0.5:
            raise Exception("\nThe maximum kinship coefficient in KING is 0.5")
        relatedness_mt = hl.king(in_mt.GT)
        relatedness_ht = relatedness_mt.filter_entries((relatedness_mt.s_1 != relatedness_mt.s) &
                                                       (relatedness_mt.phi >= kin_estimate)).entries()

        print('exporting relatedness statistics to a tsv file')
        relatedness_ht.export(f'{outdir}relatedness_checks_king.tsv.bgz')

        print('getting related samples to be removed using maximal independent set')
        samples_to_remove = hl.maximal_independent_set(relatedness_ht.s_1, relatedness_ht.s, False)
        samples = samples_to_remove.node.collect()

    if len(samples) > 0:
        # Do not remove samples that fail relatedness check
        # in_mt = in_mt.filter_cols(hl.literal(samples).contains(in_mt['s']), keep=False)
        print(f"\nNumber of samples that fail relatedness checks: {len(samples)}")
        
        df = pd.DataFrame(samples, columns=['Sample'])
        ht = hl.Table.from_pandas(df)
        ht.export(f'{outdir}samples_failing_relatedness_checks.tsv')

    else:
        print("\nNo samples failed the relatedness check")

    return in_mt, samples
