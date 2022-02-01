__author__ = 'Lindo Nkambule'

import hail as hl
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from gwaspy.pca.pca_filter_snps import pca_filter_mt, relatedness_check


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
        relatedness_method: str = 'pc_relate',
        relatedness_thresh: float = 0.98,
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

    mt = relatedness_check(in_mt=mt, method=relatedness_method, outdir=out_dir, kin_estimate=relatedness_thresh)

    pca_snps = mt.count_rows()
    if pca_snps > 1000000:
        import warnings
        warnings.warn(f'Too many SNPs to be used in PCA: {pca_snps}. This will make PCA run longer')

    print('\nRunning PCA')
    eigenvalues, pcs, _ = hl.hwe_normalized_pca(mt.GT, k=n_pcs)

    pcs_ht = pcs.transmute(**{f'PC{i}': pcs.scores[i - 1] for i in range(1, n_pcs+1)})

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
    out_scores_file = f'{out_dir}GWASpy/PCA/pca_normal/{basename}.pca.normal.scores.tsv'
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

    pdf = PdfPages('/tmp/pca.no.ref.plots.pdf')
    for figname, figure in figs_dict.items():
        pdf.savefig(figure)
    pdf.close()
    hl.hadoop_copy('file:///tmp/pca.no.ref.plots.pdf',
                   f'{out_dir}GWASpy/PCA/pca_normal/{basename}.pca.no.ref.plots.pdf')

