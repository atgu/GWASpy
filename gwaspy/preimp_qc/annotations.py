__author__ = 'Lindo Nkambule & Zan Koenig'

import hail as hl
from gwaspy.preimp_qc.aggregators import agg_call_rate, variant_qc_aggregator, impute_sex_aggregator
from gwaspy.preimp_qc.plots import plt_hist, fstat_plot, qqplot, manhattan_plot
import pandas as pd


class BaseFilter:
    def __init__(self):
        pass

    def filter(self, mt):
        pass

    def plot(self, mt):
        pass


class pre_geno(BaseFilter):
    def __init__(self, pre_geno_cr: float = 0.95, pre_row_filter: str = None, pre_col_filter: str = None):
        super().__init__()
        self._pre_geno_cr = pre_geno_cr
        self._row_filter = pre_row_filter
        self._col_filter = pre_col_filter

    def filter(self, mt):
        row_filter = mt[self._row_filter].filters if self._row_filter else mt.exclude_row
        col_filter = mt[self._col_filter].filters if self._col_filter else mt.exclude_col

        pre_filter = row_filter | col_filter

        mt = mt.annotate_rows(**{
            'pre_geno_noy': hl.struct(
                filters=hl.agg.filter(((pre_filter == False) & (mt.geno_y_excluded == True)),
                                      variant_qc_aggregator(mt).call_rate) < self._pre_geno_cr),
            'pre_geno_y': hl.struct(
                filters=hl.agg.filter(((pre_filter == False) & (mt.geno_y_only == True)),
                                      variant_qc_aggregator(mt).call_rate) < self._pre_geno_cr)})

        mt = mt.annotate_rows(**{
            'pre_geno': hl.struct(
                filters=((hl.agg.any(mt['pre_geno_noy'].filters) == True) |
                         (hl.agg.any(mt['pre_geno_y'].filters) == True))
            )})

        return mt

    def plot(self, mt):
        pass


class id_call_rate(BaseFilter):
    def __init__(self, mind: float = 0.98, pre_row_filter: str = None, pre_col_filter: str = None, data_type: str = None):
        super().__init__()
        self._mind = mind
        self._row_filter = pre_row_filter
        self._col_filter = pre_col_filter
        self._data_type = data_type

    def filter(self, mt):
        row_filter = mt[self._row_filter].filters if self._row_filter else mt.exclude_row
        col_filter = mt[self._col_filter].filters if self._col_filter else mt.exclude_col

        pre_filter = row_filter | col_filter

        mt = mt.annotate_cols(**{
            'mind': hl.struct(
                filters=hl.agg.filter(pre_filter == False, agg_call_rate(mt)) < self._mind)})

        return mt

    def plot(self, mt):
        global id_call_rate_plts
        mt = mt.annotate_cols(
            mind_cr_pre=hl.agg.filter(mt.pre_geno.filters == False, agg_call_rate(mt)))

        if 'is_case' in mt.col:
            if self._data_type == "Case-only":
                mt_cases: hl.MatrixTable = mt.filter_cols(mt.is_case == True)
                cas_pre = plt_hist(mt_cases.mind_cr_pre, title="Cases", threshold=self._mind, x_label='Call Rate')
                id_call_rate_plts = [cas_pre]

            if self._data_type == "Control-only":
                mt_controls: hl.MatrixTable = mt.filter_cols(mt.is_case == False)
                con_pre = plt_hist(mt_controls.mind_cr_pre, title="Controls", threshold=self._mind, x_label='Call Rate')
                id_call_rate_plts = [con_pre]

            if self._data_type == "Case-Control":
                mt_controls: hl.MatrixTable = mt.filter_cols(mt.is_case == False)
                mt_cases: hl.MatrixTable = mt.filter_cols(mt.is_case == True)

                con_pre = plt_hist(mt_controls.mind_cr_pre, title="Controls", threshold=self._mind, x_label='Call Rate')
                cas_pre = plt_hist(mt_cases.mind_cr_pre, title="Cases", threshold=self._mind, x_label='Call Rate')

                id_call_rate_plts = [con_pre, cas_pre]

        else:
            all_cas_con = plt_hist(mt.mind_cr_pre, title="Cases+Controls", threshold=self._mind, x_label='Call Rate')
            id_call_rate_plts = [all_cas_con]

        return id_call_rate_plts


class fhet_autosomes(BaseFilter):
    def __init__(self, fhet_thresh: float = 0.2, pre_row_filter: str = None,
                 pre_col_filter: str = None):
        super().__init__()
        self._fhet_th = fhet_thresh
        self._row_filter = pre_row_filter
        self._col_filter = pre_col_filter

    def filter(self, mt):
        mt = mt.annotate_rows(variant_qc=variant_qc_aggregator(mt))

        row_filter = mt[self._row_filter].filters if self._row_filter else mt.exclude_row
        col_filter = mt[self._col_filter].filters if self._col_filter else mt.exclude_col

        pre_filter = row_filter | col_filter

        mt = mt.annotate_cols(**{
            'fstat': hl.struct(
                filters=hl.agg.filter(pre_filter == False & mt.locus.in_autosome(),
                                      (hl.agg.inbreeding(mt.GT, hl.min(mt.variant_qc.AF)).f_stat < -self._fhet_th) |
                                      (hl.agg.inbreeding(mt.GT, hl.min(mt.variant_qc.AF)).f_stat > self._fhet_th))
            )})

        return mt

    def plot(self, mt):
        pass


class fhet_sex(BaseFilter):
    def __init__(self, fstat_x: float = 0.5, fstat_y: float = 0.5, pre_row_filter: str = None,
                 pre_col_filter: str = None, figsize: tuple = (12, 8)):
        super().__init__()
        self._fstat_x = fstat_x
        self._fstat_y = fstat_y
        self._row_filter = pre_row_filter
        self._col_filter = pre_col_filter
        self._figsize = figsize

    def filter(self, mt):
        mt = mt.annotate_rows(aaf=hl.agg.call_stats(mt.GT, mt.alleles).AF[1])
        row_filter = mt[self._row_filter].filters if self._row_filter else mt.exclude_row
        col_filter = mt[self._col_filter].filters if self._col_filter else mt.exclude_col

        pre_filter = row_filter | col_filter

        mt = mt.annotate_cols(**{
            'sex_violations': hl.struct(
                filters=hl.agg.filter(pre_filter == False,
                                      ((hl.agg.filter(mt.is_female == True,
                                                      impute_sex_aggregator(mt.GT, mt.aaf).f_stat)) > self._fstat_y) |
                                      ((hl.agg.filter(mt.is_female == False,
                                                      impute_sex_aggregator(mt.GT, mt.aaf).f_stat)) < self._fstat_x))
            )})

        return mt

    def plot(self, mt):
        mt = mt.annotate_cols(**{
            'fstat_sex': hl.struct(
                filters=hl.agg.filter(mt.pre_geno.filters == False, impute_sex_aggregator(mt.GT, mt.aaf).f_stat)
            )})

        mt_xy: hl.MatrixTable = mt.filter_cols(mt.is_female == True)
        mt_xx: hl.MatrixTable = mt.filter_cols(mt.is_female == False)

        exprs_xy = mt_xy.fstat_sex.collect()
        exprs_xx = mt_xx.fstat_sex.collect()

        df_xy = pd.DataFrame(exprs_xy)
        df_xx = pd.DataFrame(exprs_xx)

        fig = fstat_plot(df_female=df_xy, df_male=df_xx, f_stat_y=self._fstat_y, f_stat_x=self._fstat_x,
                         figsize=self._figsize)

        return fig


class fhet_sex_warnings(BaseFilter):
    def __init__(self, fstat_x: float = 0.8, fstat_y: float = 0.2, pre_row_filter: str = None,
                 pre_col_filter: str = None):
        super().__init__()
        self._fstat_x = fstat_x
        self._fstat_y = fstat_y
        self._row_filter = pre_row_filter
        self._col_filter = pre_col_filter

    def filter(self, mt):
        mt = mt.annotate_rows(aaf=hl.agg.call_stats(mt.GT, mt.alleles).AF[1])

        row_filter = mt[self._row_filter].filters if self._row_filter else mt.exclude_row
        col_filter = mt[self._col_filter].filters if self._col_filter else mt.exclude_col

        pre_filter = row_filter | col_filter

        # sex warnings are for ambiguous genotypes (F_male < 0.8, F_female > 0.2) and undefined phenotypes
        mt = mt.annotate_cols(**{
            'sex_ambiguous': hl.struct(
                filters=hl.agg.filter(pre_filter == False,
                                      ((hl.agg.filter(mt.is_female == True,
                                                      impute_sex_aggregator(mt.GT, mt.aaf).f_stat)) > self._fstat_y) |
                                      ((hl.agg.filter(mt.is_female == False,
                                                      impute_sex_aggregator(mt.GT, mt.aaf).f_stat)) < self._fstat_x))
            )})

        if 'is_case' in mt.col:
            mt = mt.annotate_cols(**{
                'sex_warnings': hl.struct(
                    filters=((hl.agg.any(mt['sex_ambiguous'].filters) == True) |
                             (hl.agg.any(hl.is_missing(mt.is_case)))
                             ))})
        else:
            mt = mt.annotate_cols(**{
                'sex_warnings': hl.struct(
                    filters=((hl.agg.any(mt['sex_ambiguous'].filters) == True)
                             ))})

        return mt

    def plot(self, mt):
        pass


class geno(BaseFilter):
    def __init__(self, geno_thresh: float = 0.98, pre_row_filter: str = None,
                 pre_col_filter: str = None, data_type: str = None):
        super().__init__()
        self._geno = geno_thresh
        self._row_filter = pre_row_filter
        self._col_filter = pre_col_filter
        self._data_type = data_type

    def filter(self, mt):
        row_filter = mt[self._row_filter].filters if self._row_filter else mt.exclude_row
        col_filter = mt[self._col_filter].filters if self._col_filter else mt.exclude_col

        pre_filter = row_filter | col_filter

        mt = mt.annotate_rows(**{
            'geno_noy': hl.struct(
                filters=hl.agg.filter(((pre_filter == False) & (mt.geno_y_excluded == True)),
                                      variant_qc_aggregator(mt).call_rate) < self._geno),
            'geno_y': hl.struct(
                filters=hl.agg.filter(((pre_filter == False) & (mt.geno_y_only == True)),
                                      variant_qc_aggregator(mt).call_rate) < self._geno)})

        mt = mt.annotate_rows(**{
            'geno': hl.struct(
                filters=((hl.agg.any(mt['geno_noy'].filters) == True) |
                         (hl.agg.any(mt['geno_y'].filters) == True))
            )})

        return mt

    def plot(self, mt):
        global var_call_rate_plts
        mt = mt.annotate_rows(
            var_cr_pre=hl.agg.filter(mt.pre_geno.filters == False, agg_call_rate(mt)))

        if 'is_case' in mt.col:
            if self._data_type == "Case-only":
                mt_cases: hl.MatrixTable = mt.filter_cols(mt.is_case == True)
                cas_pre = plt_hist(mt_cases.var_cr_pre, title="Cases", threshold=self._geno, x_label='Call Rate')
                var_call_rate_plts = [cas_pre]

            if self._data_type == "Control-only":
                mt_controls: hl.MatrixTable = mt.filter_cols(mt.is_case == False)
                con_pre = plt_hist(mt_controls.var_cr_pre, title="Controls", threshold=self._geno, x_label='Call Rate')
                var_call_rate_plts = [con_pre]

            if self._data_type == "Case-Control":
                mt_controls: hl.MatrixTable = mt.filter_cols(mt.is_case == False)
                mt_cases: hl.MatrixTable = mt.filter_cols(mt.is_case == True)

                con_pre = plt_hist(mt_controls.var_cr_pre, title="Controls", threshold=self._geno, x_label='Call Rate')
                cas_pre = plt_hist(mt_cases.var_cr_pre, title="Cases", threshold=self._geno, x_label='Call Rate')

                var_call_rate_plts = [con_pre, cas_pre]

        else:
            all_cas_con = plt_hist(mt.var_cr_pre, title="Cases+Controls", threshold=self._geno, x_label='Call Rate')
            var_call_rate_plts = [all_cas_con]

        return var_call_rate_plts


class call_rate_diff(BaseFilter):
    def __init__(self, pre_col_filter: str = None, pre_row_filter: str = None,
                 initial_row_filter: str = None, cr_thresh: float = 0.02):
        super().__init__()
        self._cr_thresh = cr_thresh
        self._row_filter = pre_row_filter
        self._col_filter = pre_col_filter
        self._initial_row_filter = initial_row_filter

    def filter(self, mt):
        row_filter = mt[self._row_filter].filters if self._row_filter else mt.exclude_row
        col_filter = mt[self._col_filter].filters if self._col_filter else mt.exclude_col

        mt = mt.annotate_rows(cr=hl.or_missing(row_filter == False,
                                               hl.agg.group_by(mt.is_case,
                                                               hl.agg.filter(col_filter == False,
                                                                             variant_qc_aggregator(mt).call_rate))))

        mt = mt.annotate_rows(diff=hl.abs(mt.cr[False] - mt.cr[True]))

        mt = mt.annotate_rows(**{
            'cr_diff': hl.struct(
                filters=hl.agg.any((mt.diff > self._cr_thresh) & (mt[self._initial_row_filter].filters == False)))})

        return mt

    def plot(self, mt):
        pass


class invariant(BaseFilter):
    def __init__(self, pre_col_filter: str = None):
        super().__init__()
        self._col_filter = pre_col_filter

    def filter(self, mt):
        col_filter = mt[self._col_filter].filters if self._col_filter else mt.exclude_col

        pre_filter = col_filter

        mt = mt.annotate_rows(**{
            'monomorphic_var': hl.struct(
                filters=hl.agg.filter(pre_filter == False, hl.min(variant_qc_aggregator(mt).AC)) == 0)})

        return mt

    def plot(self, mt):
        pass


class maf(BaseFilter):
    def __init__(self, maf_thresh: float = 0.01, pre_col_filter: str = None):
        super().__init__()
        self._maf_thresh = maf_thresh
        self._col_filter = pre_col_filter

    def filter(self, mt):
        col_filter = mt[self._col_filter].filters if self._col_filter else mt.exclude_col

        pre_filter = col_filter

        mt = mt.annotate_rows(**{
            'maf': hl.struct(
                filters=hl.agg.filter(pre_filter == False, hl.min(variant_qc_aggregator(mt).AF)) < self._maf_thresh)})

        return mt

    def plot(self, mt):
        pass


class hwe_con(BaseFilter):
    def __init__(self, hwe_th_co: float = 1e-06, pre_col_filter: str = None,
                 pre_row_filter: str = None):
        super().__init__()
        self._hwe_th_co = hwe_th_co
        self._col_filter = pre_col_filter
        self._row_filter = pre_row_filter

    def filter(self, mt):
        row_filter = mt[self._row_filter].filters if self._row_filter else mt.exclude_row
        col_filter = mt[self._col_filter].filters if self._col_filter else mt.exclude_col

        pre_filter = row_filter | col_filter

        mt = mt.annotate_rows(**{
            'hwe_con_aut': hl.struct(
                filters=((row_filter == False) &
                         (hl.agg.filter(((pre_filter == False) &
                                         (mt.is_case == False) &
                                         (mt.hwe_aut == True)),
                                        variant_qc_aggregator(mt).p_value_hwe) < self._hwe_th_co))),
            'hwe_con_sex': hl.struct(
                filters=((row_filter == False) &
                         (hl.agg.filter(((pre_filter == False) &
                                         (mt.is_case == False) &
                                         (mt.hwe_sex == True)),
                                        variant_qc_aggregator(mt).p_value_hwe) < self._hwe_th_co)))
        })

        mt = mt.annotate_rows(**{
            'hwe_con': hl.struct(
                filters=((hl.agg.any(mt['hwe_con_aut'].filters) == True) |
                         (hl.agg.any(mt['hwe_con_sex'].filters) == True))
            )})

        return mt

    def plot(self, mt):
        pass


class hwe_cas(BaseFilter):
    def __init__(self, hwe_th_ca: float = 1e-10, pre_col_filter: str = None,
                 pre_row_filter: str = None):
        super().__init__()
        self._hwe_th_ca = hwe_th_ca
        self._col_filter = pre_col_filter
        self._row_filter = pre_row_filter

    def filter(self, mt):
        row_filter = mt[self._row_filter].filters if self._row_filter else mt.exclude_row
        col_filter = mt[self._col_filter].filters if self._col_filter else mt.exclude_col

        pre_filter = row_filter | col_filter

        mt = mt.annotate_rows(**{
            'hwe_cas_aut': hl.struct(
                filters=((row_filter == False) &
                         (hl.agg.filter(((pre_filter == False) &
                                         (mt.is_case == True) &
                                         (mt.hwe_aut == True)),
                                        variant_qc_aggregator(mt).p_value_hwe) < self._hwe_th_ca))),
            'hwe_cas_sex': hl.struct(
                filters=((row_filter == False) &
                         (hl.agg.filter(((pre_filter == False) &
                                         (mt.is_case == True) &
                                         (mt.hwe_sex == True)),
                                        variant_qc_aggregator(mt).p_value_hwe) < self._hwe_th_ca)))
        })

        mt = mt.annotate_rows(**{
            'hwe_cas': hl.struct(
                filters=((hl.agg.any(mt['hwe_cas_aut'].filters) == True) |
                         (hl.agg.any(mt['hwe_cas_sex'].filters) == True))
            )})

        return mt

    def plot(self, mt):
        pass


class hwe_all(BaseFilter):
    def __init__(self, hwe_th_all: float = 1e-06, pre_col_filter: str = None,
                 pre_row_filter: str = None):
        super().__init__()
        self._hwe_th_all = hwe_th_all
        self._col_filter = pre_col_filter
        self._row_filter = pre_row_filter

    def filter(self, mt):
        row_filter = mt[self._row_filter].filters if self._row_filter else mt.exclude_row
        col_filter = mt[self._col_filter].filters if self._col_filter else mt.exclude_col

        pre_filter = row_filter | col_filter

        mt = mt.annotate_rows(**{
            'hwe_all_aut': hl.struct(
                filters=((row_filter == False) &
                         (hl.agg.filter(((pre_filter == False) &
                                         (mt.hwe_aut == True)),
                                        variant_qc_aggregator(mt).p_value_hwe) < self._hwe_th_all))),
            'hwe_all_sex': hl.struct(
                filters=((row_filter == False) &
                         (hl.agg.filter(((pre_filter == False) &
                                         (mt.hwe_sex == True)),
                                        variant_qc_aggregator(mt).p_value_hwe) < self._hwe_th_all)))
        })

        mt = mt.annotate_rows(**{
            'hwe_all': hl.struct(
                filters=((hl.agg.any(mt['hwe_all_aut'].filters) == True) |
                         (hl.agg.any(mt['hwe_all_sex'].filters) == True))
            )})

        return mt

    def plot(self, mt):
        pass


class manhattan(BaseFilter):
    def __init__(self, qqtitle, mantitle):
        super().__init__()
        self._qqtitle = qqtitle
        self._mantitle = mantitle

    def filter(self, mt):
        gwas = hl.linear_regression_rows(y=mt.is_case, x=mt.GT.n_alt_alleles(), covariates=[1.0])
        n_sig_variants = gwas.filter(gwas.p_value < 5E-8).count()

        return gwas, n_sig_variants

    def plot(self, ht):
        qq_plot, lambda_gc = qqplot(ht.p_value, title=self._qqtitle)
        man_plot = manhattan_plot(ht.p_value, title=self._mantitle)

        return qq_plot, lambda_gc, man_plot
