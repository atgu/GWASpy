__author__ = 'Lindo Nkambule'

from annotations import *
from typing import Tuple, Any, Dict
from in_out import read_mt, read_vcf, read_plink
import argparse
from reports import MyDocument
import shutil
import warnings

warnings.simplefilter(action="ignore", category=RuntimeWarning)


def summary_stats(mt: hl.MatrixTable) -> Tuple[hl.MatrixTable, Dict[str, Any]]:
    results = {}

    counts = mt.aggregate_cols(hl.struct(is_case=hl.agg.counter(mt.is_case),
                                         is_female=hl.agg.counter(mt.is_female)))

    mapping = {
        True: 'case',
        False: 'control',
        None: 'unknown'
    }
    is_case_counts = {mapping[pheno]: count for pheno, count in counts['is_case'].items()}
    results['is_case_counts'] = is_case_counts

    mapping = {
        True: 'female',
        False: 'male',
        None: 'unknown'
    }
    is_female_counts = {mapping[pheno]: count for pheno, count in counts['is_female'].items()}
    results['is_female_counts'] = is_female_counts

    n_variants = mt.count_rows()
    n_samples = mt.count_cols()

    results['n_variants'] = n_variants
    results['n_samples'] = n_samples

    pheno_status = ['case', 'control', 'unknown']
    sex_status = ['female', 'male', 'unknown']

    for i in pheno_status:
        if i not in results['is_case_counts']:
            results['is_case_counts'][i] = 0

    for i in sex_status:
        if i not in results['is_female_counts']:
            results['is_female_counts'][i] = 0

    return mt, results


def preimp_qc(mt, dirname, basename, pre_geno_thresh, mind_thresh, fhet_aut, fstat_x, fstat_y, geno_thresh,
              cr_diff_thresh, maf_thresh, hwe_th_con_thresh, hwe_th_cas_thresh, report: bool = True):
    mt = mt.annotate_rows(exclude_row=False)
    mt = mt.annotate_cols(exclude_col=False)

    mt, pre_qc_counts = summary_stats(mt)

    mt = pre_geno(pre_geno_cr=pre_geno_thresh).filter(mt)
    mt = id_call_rate(mind=mind_thresh, pre_row_filter='pre_geno').filter(mt)
    mt = fhet_autosomes(pre_row_filter='pre_geno', fhet_thresh=fhet_aut).filter(mt)
    mt = fhet_sex(pre_row_filter='pre_geno', fstat_x=fstat_x, fstat_y=fstat_y).filter(mt)
    mt = fhet_sex_warnings(pre_row_filter='pre_geno', pre_col_filter='sex_violations').filter(mt)

    mt = mt.annotate_cols(**{
        'id_pass': hl.struct(
            filters=((hl.agg.any(mt['mind'].filters) == True) | (hl.agg.any(mt['fstat'].filters) == True) |
                     (hl.agg.any(mt['sex_violations'].filters) == True))
        )})

    mt = geno(pre_row_filter='pre_geno', pre_col_filter='id_pass', geno_thresh=geno_thresh).filter(mt)
    mt = call_rate_diff(pre_row_filter='geno', pre_col_filter='id_pass', initial_row_filter='pre_geno',
                        cr_thresh=cr_diff_thresh).filter(mt)
    mt = invariant(pre_col_filter='id_pass').filter(mt)
    mt = hwe_con(pre_col_filter='id_pass', pre_row_filter='geno', hwe_th_co=hwe_th_con_thresh).filter(mt)
    mt = hwe_cas(pre_col_filter='id_pass', pre_row_filter='geno', hwe_th_ca=hwe_th_cas_thresh).filter(mt)

    results = {}
    column_filters = ['mind', 'fstat', 'sex_violations', 'sex_warnings']
    row_filters = ['pre_geno', 'geno', 'cr_diff', 'monomorphic_var', 'hwe_con', 'hwe_cas']

    mt.select_entries().select_rows(*row_filters).select_cols(*column_filters).write('{}temp.mt'.format(dirname),
                                                                                     overwrite=True)
    mt_temp = hl.read_matrix_table('{}temp.mt'.format(dirname))
    column_aggregations = mt_temp.aggregate_cols(
        [hl.agg.counter(mt_temp[filter].filters) for filter in column_filters])

    row_aggregations = mt_temp.aggregate_rows(
        [hl.agg.counter(mt_temp[filter].filters) for filter in row_filters])

    shutil.rmtree('{}temp.mt'.format(dirname))

    for filt, cont in zip(column_filters, column_aggregations):
        results[filt] = cont

    for filt, cont in zip(row_filters, row_aggregations):
        results[filt] = cont

    filters = ['pre_geno', 'mind', 'fstat', 'sex_violations', 'sex_warnings', 'geno', 'cr_diff',
               'monomorphic_var', 'maf', 'hwe_con', 'hwe_cas']

    for i in filters:
        for key, value in results.items():
            if i == key:
                print(key, ': ', value)

    # report
    if report:
        print("Writing report")
        doc = MyDocument(basename=basename, pre_qc_conts=pre_qc_counts, post_qc_conts=pre_qc_counts, count_results=results)
        doc.general_info()
        fstat_fig = fhet_sex(pre_row_filter='pre_geno', fstat_x=fstat_x, fstat_y=fstat_y, figsize=(15, 20)).plot(mt)
        fstat_fig.savefig('/tmp/fstat_fig.png', dpi=300)
        doc.individual_char(mt, mind_threshold=mind_thresh, fstat_fig_path='/tmp/fstat_fig.png')
        doc.generate_pdf('report', clean_tex=False)


def main():
    parser = argparse.ArgumentParser(description='preimp_qc V1.0')
    parser.add_argument('--dirname', type=str, required=True)
    parser.add_argument('--basename', type=str, required=True)
    parser.add_argument('--input-type', type=str, required=True, choices=['vcf', 'plink', 'hail'])
    parser.add_argument('--annotations', type=str)
    parser.add_argument('--reference', type=str, default='GRCh38')
    parser.add_argument('--report', action='store_false')
    # parser.add_argument('--qc_round', type=str, required=True)

    # required for QC
    parser.add_argument('--pre-geno', type=float, default=0.95,
                        help="include only SNPs with missing-rate < NUM (before "
                             "ID filter), important for post merge of multiple "
                             "platforms")
    parser.add_argument('--mind', type=float, default=0.98, help="include only IDs with missing-rate < NUM")
    parser.add_argument('--fhet-aut', type=float, default=0.2, help="include only IDs within NUM < FHET < NUM")
    parser.add_argument('--fstat-y', type=float, default=0.5, help="include only female IDs with fhet < NUM")
    parser.add_argument('--fstat-x', type=float, default=0.5, help="include only male IDs with fhet > NUM")
    parser.add_argument('--geno', type=float, default=0.98, help="include only SNPs with missing-rate < NUM")
    parser.add_argument('--midi', type=float, default=0.02, help="include only SNPs with missing-rate-difference ("
                                                                 "case/control) < NUM")
    parser.add_argument('--withpna', type=int, default=0, choices=[0, 1], help="include monomorphic (invariant) SNPs")
    parser.add_argument('--maf', type=float, default=0.01, help="include only SNPs with MAF >= NUM")
    parser.add_argument('--hwe-th-con', type=float, default=1e-6, help="HWE_controls < NUM")
    parser.add_argument('--hwe-th-cas', type=float, default=1e-10, help="HWE_cases < NUM")

    arg = parser.parse_args()

    # read input
    if arg.input_type == 'plink':
        input_mt = read_plink(arg.dirname, arg.basename)

    elif arg.input_type == 'vcf':
        input_mt = read_vcf(arg.dirname, arg.vcf, arg.annotations)

    else:
        input_mt = read_mt(arg.dirname, arg.basename)

    print("Running QC")
    preimp_qc(input_mt, arg.dirname, arg.basename, arg.pre_geno, arg.mind, arg.fhet_aut, arg.fstat_x,
              arg.fstat_y, arg.geno, arg.midi, arg.maf, arg.hwe_th_con, arg.hwe_th_cas, report=arg.report)

    print("\nDone running QC!")


if __name__ == '__main__':
    main()


# python3 preimp_qc.py --dirname ../qcData/ --basename sim_sim2a_eur_sa_merge.miss --input-type plink

# exclude report
# python3 preimp_qc.py --dirname ../qcData/ --basename sim_sim2a_eur_sa_merge.miss --input-type plink --report