__author__ = 'Dan King'

import hail as hl


def variant_qc_aggregator(mt) -> hl.MatrixTable:
    """:func:`.variant_qc` as an aggregator."""
    bound_exprs = {}
    gq_dp_exprs = {}

    def has_field_of_type(name, dtype):
        return name in mt.entry and mt[name].dtype == dtype

    if has_field_of_type('DP', hl.tint32):
        gq_dp_exprs['dp_stats'] = hl.agg.stats(mt.DP).select('mean', 'stdev', 'min', 'max')
    if has_field_of_type('GQ', hl.tint32):
        gq_dp_exprs['gq_stats'] = hl.agg.stats(mt.GQ).select('mean', 'stdev', 'min', 'max')
    if not has_field_of_type('GT', hl.tcall):
        raise ValueError("'variant_qc': expect an entry field 'GT' of type 'call'")
    bound_exprs['n_called'] = hl.agg.count_where(hl.is_defined(mt['GT']))
    bound_exprs['n_not_called'] = hl.agg.count_where(hl.is_missing(mt['GT']))
    n_cols = hl.agg.count()
    bound_exprs['n_filtered'] = hl.int64(n_cols) - hl.agg.count()
    bound_exprs['call_stats'] = hl.agg.call_stats(mt.GT, mt.alleles)
    return hl.rbind(hl.struct(**bound_exprs),
                    lambda e1: hl.rbind(
                        hl.case().when(hl.len(mt.alleles) == 2,
                                       hl.hardy_weinberg_test(e1.call_stats.homozygote_count[0],
                                                              e1.call_stats.AC[1] - 2
                                                              * e1.call_stats.homozygote_count[1],
                                                              e1.call_stats.homozygote_count[1])
                                       ).or_missing(),
                        lambda hwe: hl.struct(**{
                            **gq_dp_exprs,
                            **e1.call_stats,
                            'call_rate': hl.float(e1.n_called) / (e1.n_called + e1.n_not_called + e1.n_filtered),
                            'n_called': e1.n_called,
                            'n_not_called': e1.n_not_called,
                            'n_filtered': e1.n_filtered,
                            'n_het': e1.n_called - hl.sum(e1.call_stats.homozygote_count),
                            'n_non_ref': e1.n_called - e1.call_stats.homozygote_count[0],
                            'het_freq_hwe': hwe.het_freq_hwe,
                            'p_value_hwe': hwe.p_value})))


def agg_call_rate(mt: hl.MatrixTable):
    # DOES NOT HANDLE filter_entries CORRECTLY!
    n_called = hl.agg.count_where(hl.is_defined(mt['GT']))

    return hl.agg.filter(
        ~(mt.exclude_row | mt.exclude_col),
        n_called / hl.agg.count())


def impute_sex_aggregator(call,
                          aaf,
                          aaf_threshold=0.0,
                          include_par=False,
                          female_threshold=0.4,
                          male_threshold=0.8) -> hl.Table:
    """:func:`.impute_sex` as an aggregator."""
    mt = call._indices.source
    rg = mt.locus.dtype.reference_genome
    x_contigs = hl.literal(
        hl.eval(
            hl.map(lambda x_contig: hl.parse_locus_interval(x_contig, rg), rg.x_contigs)))
    inbreeding = hl.agg.inbreeding(call, aaf)
    is_female = hl.if_else(inbreeding.f_stat < female_threshold,
                           True,
                           hl.if_else(inbreeding.f_stat > male_threshold,
                                      False,
                                      hl.is_missing('tbool')))
    expression = hl.struct(is_female=is_female, **inbreeding)
    if not include_par:
        interval_type = hl.tarray(hl.tinterval(hl.tlocus(rg)))
        par_intervals = hl.literal(rg.par, interval_type)
        expression = hl.agg.filter(
            ~par_intervals.any(lambda par_interval: par_interval.contains(mt.locus)),
            expression)
    expression = hl.agg.filter((aaf > aaf_threshold) & (aaf < (1 - aaf_threshold)), expression)
    expression = hl.agg.filter(
        x_contigs.any(lambda contig: contig.contains(mt.locus)),
        expression)

    return expression


def allele_types(mt):
    from hail.expr.functions import _num_allele_type, _allele_types
    allele_types = _allele_types[:]
    allele_types.extend(['Transition', 'Transversion'])
    allele_enum = {i: v for i, v in enumerate(allele_types)}
    allele_ints = {v: k for k, v in allele_enum.items()}

    def allele_type(ref, alt):
        return hl.bind(lambda at: hl.if_else(at == allele_ints['SNP'],
                                             hl.if_else(hl.is_transition(ref, alt),
                                                        allele_ints['Transition'],
                                                        allele_ints['Transversion']),
                                             at),
                       _num_allele_type(ref, alt))

    return mt.alleles[1:].map(lambda alt: allele_type(mt.alleles[0], alt))
