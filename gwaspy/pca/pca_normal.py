__author__ = 'Lindo Nkambule'

import hail as hl
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from gwaspy.pca.pca_filter_snps import pca_filter_mt, relatedness_check


def pc_project(
        mt: hl.MatrixTable = None,
        loadings_ht: hl.Table = None,
        loading_location: str = 'loadings',
        af_location: str = 'pca_af') -> hl.Table:
    """
    Projects samples in `mt` on pre-computed PCs.

    :param mt: MT containing the samples to project
    :param loadings_ht: HT containing the PCA loadings and allele frequencies used for the PCA
    :param loading_location: Location of expression for loadings in `loadings_ht`
    :param af_location: Location of expression for allele frequency in `loadings_ht`
    :return: Table with scores calculated from loadings in column `scores`
    """

    n_variants = loadings_ht.count()

    mt = mt.annotate_rows(
        pca_loadings=loadings_ht[mt.row_key][loading_location],
        pca_af=loadings_ht[mt.row_key][af_location],
    )

    mt = mt.filter_rows(
        hl.is_defined(mt.pca_loadings)
        & hl.is_defined(mt.pca_af)
        & (mt.pca_af > 0)
        & (mt.pca_af < 1)
    )

    gt_norm = (mt.GT.n_alt_alleles() - 2 * mt.pca_af) / hl.sqrt(
        n_variants * 2 * mt.pca_af * (1 - mt.pca_af)
    )

    mt = mt.annotate_cols(scores=hl.agg.array_sum(mt.pca_loadings * gt_norm))

    return mt.cols().select('scores')


def plot_pca(
        in_df: pd.DataFrame,
        x_pc: str,
        y_pc: str,
        column: str):

    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(15, 15))

    if column == 'is_female':
        cat1 = in_df[in_df['is_female'] == 'female']
        label1 = 'female'
        color1 = 'blue'
        marker1 = "^"

        cat2 = in_df[in_df['is_female'] == 'male']
        label2 = 'male'
        color2 = 'orange'
        marker2 = "s"

        cat3 = in_df[in_df['is_female'] == 'unknown']
        label3 = 'unknown'
    else:
        cat1 = in_df[in_df['is_case'] == 'case']
        label1 = 'case'
        color1 = 'green'
        marker1 = "o"

        cat2 = in_df[in_df['is_case'] == 'control']
        label2 = 'control'
        color2 = 'darkorange'
        marker2 = "o"

        cat3 = in_df[in_df['is_case'] == 'unknown']
        label3 = 'unknown'

    axs.scatter(cat1[x_pc], cat1[y_pc], c=color1, label=label1, s=15, alpha=1, marker=marker1, edgecolor='black',
                linewidths=0.4)

    axs.scatter(cat2[x_pc], cat2[y_pc], c=color2, label=label2, s=15, alpha=1, marker=marker2, edgecolor='black',
                linewidths=0.4)

    axs.scatter(cat3[x_pc], cat3[y_pc], c='red', label=label3, s=15, alpha=1, marker='x', edgecolor='black',
                linewidths=0.4)

    plt.legend(bbox_to_anchor=(1, 0.5), loc='lower left', frameon=False, prop={'size': 13})
    axs.set_xlabel(xlabel=x_pc, fontsize=15)
    axs.set_ylabel(ylabel=y_pc, fontsize=15)

    plt.close()

    return fig


def run_pca_normal(
        dirname: str = None,
        basename: str = None,
        input_type: str = None,
        reference: str = 'GRCh38',
        maf: float = 0.05,
        hwe: float = 1e-3,
        call_rate: float = 0.98,
        ld_cor: float = 0.2,
        ld_window: int = 250000,
        n_pcs: int = 20,
        run_relatedness_check: bool = True,
        include_kinself: bool = False,
        relatedness_method: str = 'pc_relate',
        relatedness_thresh: float = 0.1,
        out_dir: str = None):

    print('\nReading mt')
    if reference.lower() == 'grch37':
        lifted_over = f'{dirname}{basename}.liftover.grch38.mt'
        if not hl.hadoop_exists(lifted_over):
            from gwaspy.utils.reference_liftover import liftover_to_grch38
            mt = liftover_to_grch38(dirname=dirname, basename=basename, input_type=input_type)
        else:
            print(f'\nFound lifted-over over file: {lifted_over}')
            mt = hl.read_matrix_table(lifted_over)
    else:
        from gwaspy.utils.read_file import read_infile
        mt = read_infile(input_type=input_type, dirname=dirname, basename=basename)

    print('\nFiltering mt')
    mt = pca_filter_mt(in_mt=mt, maf=maf, hwe=hwe, call_rate=call_rate, ld_cor=ld_cor, ld_window=ld_window)

    if run_relatedness_check:
        out_dir = f'{out_dir}GWASpy/PCA/{basename}/pca_normal/'
        mt, fail_samples = relatedness_check(in_mt=mt, method=relatedness_method, outdir=out_dir,
                                             kin_estimate=relatedness_thresh, include_kinself=include_kinself)
    else:
        print('Skipping relatedness checks')
        out_dir = f'{out_dir}GWASpy/PCA/{basename}/pca_normal/'
        out_filename = f'{out_dir}samples_failing_relatedness_checks.tsv'
        if not hl.hadoop_exists(out_filename):
            print('''Did not find any previously created file with samples failing relatedness checks. All samples will
            be treated as unrelateds''')
            fail_samples = []
        else:
            print(f'Related samples from file {out_filename} will be used')
            related = pd.read_csv(out_filename)
            fail_samples = related['Sample'].to_list()
            print(f'''Found {len(fail_samples)} samples that failed relatedness checks and will be projected''')

    pca_snps = mt.count_rows()
    if pca_snps > 1000000:
        import warnings
        warnings.warn(f'Too many SNPs to be used in PCA: {pca_snps}. This will make PCA run longer')

    print('\nRunning PCA on unrelated samples')

    # filter_cols will not work if fail_samples list is empty
    if len(fail_samples) > 0:
        unrelated_mt = mt.filter_cols(hl.literal(fail_samples).contains(mt['s']), keep=False)
    else:
        unrelated_mt = mt

    # run PCA on unrelated samples
    eigenvalues, pcs, loadings = hl.hwe_normalized_pca(unrelated_mt.GT, k=n_pcs, compute_loadings=True)
    unrelated_scores = pcs.transmute(**{f'PC{i}': pcs.scores[i - 1] for i in range(1, n_pcs+1)})
    unrelated_scores = unrelated_scores.annotate(Projected='No - unrelated')
    # add AF annotation
    pca_mt = unrelated_mt.annotate_rows(pca_af=hl.agg.mean(unrelated_mt.GT.n_alt_alleles()) / 2)
    loadings = loadings.annotate(pca_af=pca_mt.rows()[loadings.key].pca_af)

    if len(fail_samples) > 0:
        print(f'\nProjecting {len(fail_samples)} related samples on PCs pre-computed using unrelated samples')
        related_mt = mt.filter_cols(hl.literal(fail_samples).contains(mt['s']), keep=True)
        related_scores = pc_project(mt=related_mt, loadings_ht=loadings)
        related_scores = related_scores.transmute(**{f'PC{i}': related_scores.scores[i - 1] for i in range(1, n_pcs+1)})
        related_scores = related_scores.annotate(Projected='Yes - related')

        # merge the related scores with unrelateds
        pcs_ht = unrelated_scores.union(related_scores)
    else:
        pcs_ht = unrelated_scores

    # add phenotype and sex to the output, using information from the mt
    # first check if is_case and os_female fields exist in the mt
    all_column_field_names = list(mt.col)
    # sex status is a MUST but not phenotype status
    if 'is_case' in all_column_field_names:
        ann_cols = ['is_case', 'is_female']
    else:
        ann_cols = ['is_female']

    annotations_ht = mt.cols().select(*ann_cols)

    if 'is_case' in all_column_field_names:
        pcs_ht = pcs_ht.annotate(is_case=annotations_ht[pcs_ht.s].is_case)
    pcs_ht = pcs_ht.annotate(is_female=annotations_ht[pcs_ht.s].is_female)

    print('\nSaving PC scores file')
    out_scores_file = f'{out_dir}{basename}.pca.normal.scores.tsv'
    pcs_ht.export(out_scores_file)

    print('\nGenerating PCA plots')
    pcs_scores = pd.read_table(out_scores_file, header=0, sep='\t')

    if 'is_case' in all_column_field_names:
        pcs_scores[['is_case']] = pcs_scores[['is_case']].replace([True, False, None], ['case', 'control', 'unknown'])
    pcs_scores[['is_female']] = pcs_scores[['is_female']].replace([True, False, None], ['female', 'male', 'unknown'])

    figs_dict = {}
    for col in ann_cols:
        for i in range(1, n_pcs, 2):
            xpc = f'PC{i}'
            ypc = f'PC{i + 1}'

            figs_dict["fig{}{}".format(col, i)] = plot_pca(pcs_scores, xpc, ypc, col)

    pdf = PdfPages('/tmp/pca.normal.plots.pdf')
    for figname, figure in figs_dict.items():
        pdf.savefig(figure)
    pdf.close()
    hl.hadoop_copy('file:///tmp/pca.normal.plots.pdf',
                   f'{out_dir}{basename}.pca.normal.plots.pdf')

