__author__ = 'Lindo Nkambule'

from gwaspy.preimp_qc.annotations import *
from typing import Tuple, Any, Dict, Union
from gwaspy.utils.read_file import read_infile
import argparse
from gwaspy.preimp_qc.report import MyDocument
import shutil
import warnings
import os
import hail as hl

warnings.simplefilter(action='ignore', category=RuntimeWarning)


def summary_stats(mt: hl.MatrixTable) -> Tuple[hl.MatrixTable, Dict[str, Any]]:
    results = {}

    if 'is_case' in mt.col:
        counts = mt.aggregate_cols(hl.struct(is_case=hl.agg.counter(mt.is_case),
                                             is_female=hl.agg.counter(mt.is_female)))

        mapping = {
            False: 'control',
            True: 'case',
            None: 'unknown'
        }
        is_case_counts = {mapping[pheno]: count for pheno, count in counts['is_case'].items()}
        results['is_case_counts'] = is_case_counts

    else:
        counts = mt.aggregate_cols(hl.struct(is_female=hl.agg.counter(mt.is_female)))

    mapping = {
        True: 'female',
        False: 'male',
        None: 'unknown'
    }
    is_female_counts = {mapping[pheno]: count for pheno, count in counts['is_female'].items()}
    results['is_female_counts'] = is_female_counts

    if 'is_case' in mt.col:
        pheno_status = ['case', 'control', 'unknown']

        for i in pheno_status:
            if i not in results['is_case_counts']:
                results['is_case_counts'][i] = 0

    sex_status = ['female', 'male', 'unknown']

    for i in sex_status:
        if i not in results['is_female_counts']:
            results['is_female_counts'][i] = 0

    n_variants = mt.count_rows()
    n_samples = mt.count_cols()

    results['n_variants'] = n_variants
    results['n_samples'] = n_samples

    return mt, results


def preimp_qc(input_type: str = None, dirname: str = None, basename: str = None, pre_geno_thresh: Union[int, float] = 0.95,
              mind_thresh: Union[int, float] = 0.98, fhet_aut: Union[int, float] = 0.2, fstat_x: Union[int, float] = 0.5,
              fstat_y: Union[int, float] = 0.5, geno_thresh: Union[int, float] = 0.98,
              cr_diff_thresh: Union[int, float] = 0.02, maf_thresh: Union[int, float] = 0.01,
              hwe_th_con_thresh: Union[int, float] = 1e-6, hwe_th_cas_thresh: Union[int, float] = 1e-10,
              hwe_th_all_thresh: Union[int, float] = 1e-06, annotations_file: str = None, report: bool = True,
              export_type: str = 'hail', out_dir: str = None, reference: str = 'GRCh38'):
    print('\nRunning QC')

    global mt, row_filters, filters, data_type, lambda_gc_pos, lambda_gc_pre, n_sig_var_pre, n_sig_var_pos, man_table_results, remove_fields

    # create temp directory for storing temp files
    if not os.path.exists('gwaspy_tmp'):
        os.makedirs('gwaspy_tmp')

    # LaTex needs full path to files
    gwaspy_dir = os.getcwd() + '/gwaspy_tmp'

    output_directory = out_dir if out_dir else dirname

    hl.init(default_reference=reference)

    # read input
    mt = read_infile(input_type=input_type, dirname=dirname, basename=basename, annotations=annotations_file)

    if 'is_case' in mt.col:
        gwas_pre, n_sig_var_pre = manhattan(qqtitle='Pre-QC QQ Plot', mantitle='Pre-QC Manhattan Plot').filter(mt)
        qqplt_pre, lambda_gc_pre, manplt_pre = manhattan(qqtitle='Pre-QC QQ Plot',
                                                         mantitle='Pre-QC Manhattan Plot').plot(gwas_pre)
        qqplt_pre.savefig('gwaspy_tmp/gwaspy_qq_pre.png', dpi=300)
        manplt_pre.savefig('gwaspy_tmp/gwaspy_man_pre.png', dpi=300)

    mt = mt.annotate_rows(exclude_row=False)
    mt = mt.annotate_cols(exclude_col=False)

    mt, pre_qc_counts = summary_stats(mt)

    if 'is_case' in mt.col:
        if (pre_qc_counts['is_case_counts']['case'] > 0) & (pre_qc_counts['is_case_counts']['control'] == 0):
            data_type = 'Case-only'
        elif (pre_qc_counts['is_case_counts']['control'] > 0) & (pre_qc_counts['is_case_counts']['case'] == 0):
            data_type = 'Control-only'
        elif (pre_qc_counts['is_case_counts']['case'] > 0) & (pre_qc_counts['is_case_counts']['control'] > 0):
            data_type = 'Case-Control'
        else:
            data_type = 'Trio'
    else:
        data_type = 'no-pheno'

    chroms = mt.aggregate_rows(hl.agg.collect_as_set(mt.locus.contig))
    if ('chrX' or 'chrY' or 'chrMT') in chroms:
        chromx, chromy, chrommt = 'chrX', 'chrY', 'chrMT'
    else:
        chromx, chromy, chrommt = 'X', 'Y', 'MT'

    # we need to compute call rate for chr1-23 and chrY separately since females have no chrY
    mt = mt.annotate_entries(
        geno_y_excluded=(hl.case()
                         .when(mt.locus.contig == chromy, False)
                         .default(True)
                         ),
        geno_y_only=(hl.case()
                     .when(mt.locus.contig == chromy, mt.is_female == False)
                     .default(False)
                     )
    )

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

    mt = geno(pre_row_filter='pre_geno', pre_col_filter='id_pass', geno_thresh=geno_thresh, data_type=data_type).filter(mt)
    mt = invariant(pre_col_filter='id_pass').filter(mt)

    # for HWE, markers in: (1) autosomes - include males+females; (2) chrX - include ONLY females; (3) exclude chrY
    mt = mt.annotate_entries(
        hwe_aut=(hl.case()
                 .when(mt.locus.contig == chromx, False)
                 .when(mt.locus.contig == chromy, False)
                 .when(mt.locus.contig == chrommt, False)
                 .default(True)
                 ),
        hwe_sex=(hl.case()
                 .when(mt.locus.contig == chromx, mt.is_female)
                 .default(False)
                 )
    )

    if 'is_case' in mt.col:
        mt = call_rate_diff(pre_row_filter='geno', pre_col_filter='id_pass', initial_row_filter='pre_geno',
                            cr_thresh=cr_diff_thresh).filter(mt)

        # check if data is case-/control-only, case-control, or trio
        # (a) Case-Only
        if data_type == 'Case-only':
            print("\n" + data_type)
            mt = hwe_cas(pre_col_filter='id_pass', pre_row_filter='geno', hwe_th_ca=1e-6).filter(mt)
            row_filters = ['pre_geno', 'geno', 'cr_diff', 'monomorphic_var', 'hwe_cas']
            filters = ['pre_geno', 'mind', 'fstat', 'sex_violations', 'sex_warnings', 'geno', 'cr_diff',
                       'monomorphic_var', 'hwe_cas']
            remove_fields = ['cr', 'diff', 'hwe_cas_aut', 'hwe_cas_sex']
        # (b) Control-Only
        elif data_type == 'Control-only':
            print("\n" + data_type)
            mt = hwe_con(pre_col_filter='id_pass', pre_row_filter='geno', hwe_th_co=1e-6).filter(mt)
            row_filters = ['pre_geno', 'geno', 'cr_diff', 'monomorphic_var', 'hwe_con']
            filters = ['pre_geno', 'mind', 'fstat', 'sex_violations', 'sex_warnings', 'geno', 'cr_diff',
                       'monomorphic_var', 'hwe_con']
            remove_fields = ['cr', 'diff', 'hwe_con_aut', 'hwe_con_sex']
        elif data_type == 'Case-Control':
            print("\n" + data_type)
            mt = hwe_cas(pre_col_filter='id_pass', pre_row_filter='geno', hwe_th_ca=hwe_th_cas_thresh).filter(mt)
            mt = hwe_con(pre_col_filter='id_pass', pre_row_filter='geno', hwe_th_co=hwe_th_con_thresh).filter(mt)
            row_filters = ['pre_geno', 'geno', 'cr_diff', 'monomorphic_var', 'hwe_con', 'hwe_cas']
            filters = ['pre_geno', 'mind', 'fstat', 'sex_violations', 'sex_warnings', 'geno', 'cr_diff',
                       'monomorphic_var', 'hwe_con', 'hwe_cas']
            remove_fields = ['cr', 'diff', 'hwe_cas_aut', 'hwe_cas_sex', 'hwe_con_aut', 'hwe_con_sex']
        else:
            # trio data
            print(data_type)
    else:
        print('Running HWE filters on whole dataset without spliting by phenotype status')
        mt = hwe_all(pre_col_filter='id_pass', pre_row_filter='geno', hwe_th_all=1e-08).filter(mt)
        row_filters = ['pre_geno', 'geno', 'monomorphic_var', 'hwe_all']
        filters = ['pre_geno', 'mind', 'fstat', 'sex_violations', 'sex_warnings', 'geno',
                   'monomorphic_var', 'hwe_all']
        remove_fields = ['hwe_all_aut', 'hwe_all_sex']

    results = {}
    column_filters = ['mind', 'fstat', 'sex_violations', 'sex_warnings']

    mt.select_entries().select_rows(*row_filters).select_cols(*column_filters).write(f'{output_directory}/temp.mt',
                                                                                     overwrite=True)
    mt_temp = hl.read_matrix_table(f'{output_directory}/temp.mt')
    column_aggregations = mt_temp.aggregate_cols(
        [hl.agg.counter(mt_temp[filter].filters) for filter in column_filters])

    row_aggregations = mt_temp.aggregate_rows(
        [hl.agg.counter(mt_temp[filter].filters) for filter in row_filters])

    for filt, cont in zip(column_filters, column_aggregations):
        results[filt] = dict(cont) # aggregate returns a frozendict, convert that back to a dict

    for filt, cont in zip(row_filters, row_aggregations):
        results[filt] = dict(cont) # # aggregate returns a frozendict, convert that back to a dict

    for i in filters:
        # some filters will have zero snps/id filtered, and there won't be a True, so add it
        if True not in results[i]:
            results[i][True] = 0
        for key, value in results.items():
            if i == key:
                print(key, ': ', value)

    if report:
        fstat_fig = fhet_sex(pre_row_filter='pre_geno', fstat_x=fstat_x, fstat_y=fstat_y, figsize=(15, 20)).plot(mt)
        fstat_fig.savefig('gwaspy_tmp/gwaspy_fstat_fig.png', dpi=300)

        id_cr_plot = id_call_rate(mind=mind_thresh, pre_row_filter='pre_geno', data_type=data_type).plot(mt)
        var_cr_plot = geno(pre_row_filter='pre_geno', pre_col_filter='id_pass', geno_thresh=geno_thresh,
                           data_type=data_type).plot(mt)

        if 'is_case' in mt.col:
            if data_type == 'Case-only':
                id_cr_plot[0].savefig('gwaspy_tmp/gwaspy_id_cas_pre.png', dpi=300)
                var_cr_plot[0].savefig('gwaspy_tmp/gwaspy_var_cas_pre.png', dpi=300)

            if data_type == 'Control-only':
                id_cr_plot[0].savefig('gwaspy_tmp/gwaspy_id_con_pre.png', dpi=300)
                var_cr_plot[0].savefig('gwaspy_tmp/gwaspy_var_con_pre.png', dpi=300)

            if data_type == 'Case-Control':
                id_cr_plot[0].savefig('gwaspy_tmp/gwaspy_id_con_pre.png', dpi=300)
                id_cr_plot[1].savefig('gwaspy_tmp/gwaspy_id_cas_pre.png', dpi=300)
                var_cr_plot[0].savefig('gwaspy_tmp/gwaspy_var_con_pre.png', dpi=300)
                var_cr_plot[1].savefig('gwaspy_tmp/gwaspy_var_cas_pre.png', dpi=300)

        else:
            id_cr_plot[0].savefig('gwaspy_tmp/gwaspy_id_cas_con_pre.png', dpi=300)
            var_cr_plot[0].savefig('gwaspy_tmp/gwaspy_var_cas_con_pre.png', dpi=300)

    # FILTER OUT ALL SNPs and IDs THAT FAIL QC
    column_filters = ['mind', 'fstat', 'sex_violations']
    for row in row_filters:
        mt = mt.filter_rows(mt[row].filters == True, keep=False)
    for col in column_filters:
        mt = mt.filter_cols(mt[col].filters == True, keep=False)

    mt_filtered, pos_qc_counts = summary_stats(mt)

    # drop entry fields we added as they will cause errors when exporting to VCF
    # e.g. Error summary: HailException: Invalid type for format field 'geno_y_excluded'. Found 'bool'.
    drop_fields = filters + ['geno_y_excluded', 'geno_y_only', 'pre_geno_noy', 'pre_geno_y', 'hwe_aut', 'hwe_sex',
                             'exclude_col', 'exclude_row', 'variant_qc', 'aaf', 'geno_noy', 'geno_y', 'sex_ambiguous',
                             'id_pass'] + remove_fields
    mt_filtered = mt_filtered.drop(*drop_fields)

    if 'is_case' in mt.col:
        gwas_pos, n_sig_var_pos = manhattan(qqtitle='Post-QC QQ Plot', mantitle='Post-QC Manhattan Plot').filter(mt)
        qqplt_pos, lambda_gc_pos, manplt_pos = manhattan(qqtitle='Post-QC QQ Plot',
                                                         mantitle='Post-QC Manhattan Plot').plot(gwas_pos)

        ncas_pre = pre_qc_counts['is_case_counts']['case']
        ncas_pos = pos_qc_counts['is_case_counts']['case']
        ncon_pre = pre_qc_counts['is_case_counts']['control']
        ncon_pos = pos_qc_counts['is_case_counts']['control']
        lambda_thous_pre = 1 + (lambda_gc_pre-1)*(1/ncas_pre+1/ncon_pre)/(1/1000+1/1000)
        lambda_thous_pos = 1 + (lambda_gc_pos-1)*(1/ncas_pos+1/ncon_pos)/(1/1000+1/1000)

        qqplt_pos.savefig('gwaspy_tmp/gwaspy_qq_pos.png', dpi=300)
        manplt_pos.savefig('gwaspy_tmp/gwaspy_man_pos.png', dpi=300)

        man_table_results = [n_sig_var_pre, n_sig_var_pos, lambda_gc_pre, lambda_gc_pos, round(lambda_thous_pre, 3),
                             round(lambda_thous_pos, 3)]

    # report
    if report:
        print('\nWriting report')
        doc = MyDocument(basename=basename)
        if 'is_case' in mt.col:
            doc.flags_table(pre_qc_counts=pre_qc_counts, pos_qc_counts=pos_qc_counts, results=results,
                            lambda_gc=lambda_gc_pos, sig_vars=n_sig_var_pos)
        else:
            doc.flags_table(pre_qc_counts=pre_qc_counts, pos_qc_counts=pos_qc_counts, results=results)
        doc.general_info(pre_qc_conts=pre_qc_counts, post_qc_conts=pos_qc_counts,
                         count_results=results, pre_filter=pre_geno_thresh, id_cr=mind_thresh, fhet_thresh=fhet_aut,
                         var_cr=geno_thresh, miss_diff=cr_diff_thresh, hwe_con=hwe_th_con_thresh,
                         hwe_cas=hwe_th_cas_thresh, hwe_all=hwe_th_all_thresh, data_type=data_type)
        if 'is_case' in mt.col:
            doc.manhattan_sec(qq_pre_path=f'{gwaspy_dir}/gwaspy_qq_pre.png', qq_pos_path=f'{gwaspy_dir}/gwaspy_qq_pos.png',
                              man_pre_path=f'{gwaspy_dir}/gwaspy_man_pre.png',
                              man_pos_path=f'{gwaspy_dir}/gwaspy_man_pos.png',
                              table_results=man_table_results)

        if data_type == 'Case-only':
            doc.individual_char(id_con_pre_path='nothing here', id_cas_pre_path=f'{gwaspy_dir}/gwaspy_id_cas_pre.png',
                                fstat_fig_path=f'{gwaspy_dir}/gwaspy_fstat_fig.png', id_all_path='nothing here',
                                data_type=data_type)
            doc.snp_char(var_con_pre_path='nothing here', var_cas_pre_path=f'{gwaspy_dir}/gwaspy_var_cas_pre.png',
                         var_all_path='nothing here', data_type=data_type)
        if data_type == 'Control-only':
            doc.individual_char(id_con_pre_path=f'{gwaspy_dir}/gwaspy_id_con_pre.png', id_cas_pre_path='nothing here',
                                fstat_fig_path=f'{gwaspy_dir}/gwaspy_fstat_fig.png', id_all_path='nothing here',
                                data_type=data_type)
            doc.snp_char(var_con_pre_path=f'{gwaspy_dir}/gwaspy_var_cas_pre.png', var_cas_pre_path='nothing here',
                         var_all_path='nothing here', data_type=data_type)
        if data_type == 'Case-Control':
            doc.individual_char(id_con_pre_path=f'{gwaspy_dir}/gwaspy_id_con_pre.png',
                                id_cas_pre_path=f'{gwaspy_dir}/gwaspy_id_cas_pre.png',
                                fstat_fig_path=f'{gwaspy_dir}/gwaspy_fstat_fig.png', id_all_path='nothing here',
                                data_type=data_type)
            doc.snp_char(var_con_pre_path=f'{gwaspy_dir}/gwaspy_var_cas_pre.png',
                         var_cas_pre_path=f'{gwaspy_dir}/gwaspy_var_cas_pre.png',
                         var_all_path='nothing here', data_type=data_type)
        if data_type == 'no-pheno':
            doc.individual_char(id_con_pre_path='nothing here', id_cas_pre_path='nothing here',
                                fstat_fig_path=f'{gwaspy_dir}/gwaspy_fstat_fig.png',
                                id_all_path=f'{gwaspy_dir}/gwaspy_id_cas_con_pre.png', data_type=data_type)
            doc.snp_char(var_con_pre_path='nothing here', var_cas_pre_path=f'nothing here',
                         var_all_path=f'{gwaspy_dir}/gwaspy_var_cas_con_pre.png', data_type=data_type)
        doc.generate_pdf(f'{gwaspy_dir}/{basename}.preimp_qc.report', clean=True, clean_tex=True)

    print('\nExporting qced file')
    if export_type:
        from gwaspy.utils.export_file import export_qced_file
        export_qced_file(mt=mt_filtered, out_dir=output_directory, basename=basename, export_type=export_type)

    if out_dir.startswith('gs://'):
        hl.hadoop_copy(f'file://{gwaspy_dir}/{basename}.preimp_qc.report.pdf',
                       f'{output_directory}GWASpy/Preimp_QC/{basename}.preimp_qc.report.pdf')
    else:
        shutil.copyfile(f'{gwaspy_dir}/{basename}.preimp_qc.report.pdf',
                        f'{output_directory}{basename}.preimp_qc.report.pdf')

    # clean-up
    print('\nCleaning up')
    shutil.rmtree('gwaspy_tmp')

    print("\nDone running QC!")


def main():
    parser = argparse.ArgumentParser(description='preimp_qc')
    parser.add_argument('--dirname', type=str, required=True)
    parser.add_argument('--basename', type=str, required=True)
    parser.add_argument('--input-type', type=str, required=True, choices=['vcf', 'plink', 'hail'])
    parser.add_argument('--export-type', type=str, default='hail', choices=['vcf', 'plink', 'hail'])
    parser.add_argument('--out-dir', type=str, default=None)
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
    parser.add_argument('--geno', type=float, default=0.98, help="include only SNPs with missing-rate > NUM")
    parser.add_argument('--midi', type=float, default=0.02, help="include only SNPs with missing-rate-difference ("
                                                                 "case/control) < NUM")
    parser.add_argument('--withpna', type=int, default=0, choices=[0, 1], help="include monomorphic (invariant) SNPs")
    parser.add_argument('--maf', type=float, default=0.01, help="include only SNPs with MAF >= NUM")
    parser.add_argument('--hwe-th-con', type=float, default=1e-06, help="HWE controls < NUM")
    parser.add_argument('--hwe-th-cas', type=float, default=1e-10, help="HWE cases < NUM")
    parser.add_argument('--hwe-th-all', type=float, default=1e-06, help="HWE cases + controls < NUM")

    arg = parser.parse_args()

    preimp_qc(input_type=arg.input_type, dirname=arg.dirname, basename=arg.basename, pre_geno_thresh=arg.pre_geno,
              mind_thresh=arg.mind, fhet_aut=arg.fhet_aut, fstat_x=arg.fstat_x, fstat_y=arg.fstat_y,
              geno_thresh=arg.geno, cr_diff_thresh=arg.midi, maf_thresh=arg.maf, hwe_th_con_thresh=arg.hwe_th_con,
              hwe_th_cas_thresh=arg.hwe_th_cas, hwe_th_all_thresh=arg.hwe_th_all, annotations_file=arg.annotations,
              report=arg.report, export_type=arg.export_type, out_dir=arg.out_dir, reference=arg.reference)


if __name__ == '__main__':
    main()

