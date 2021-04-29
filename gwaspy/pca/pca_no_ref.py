import hail as hl
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def filter_mt(
        in_mt: hl.MatrixTable,
        maf: float = 0.01):

    print("\nInitial number of SNPs before filtering: {}".format(in_mt.count_rows()))
    mt = hl.variant_qc(in_mt)
    mt_filt = mt.annotate_rows(maf=hl.min(mt.variant_qc.AF))
    mt_filt = mt_filt.filter_rows(mt_filt.maf > maf)
    # mt_filt = mt.filter_rows((mt.variant_qc.AF[0] > 0.001) & (mt.variant_qc.AF[0] < 0.999))
    print("\nNumber of SNPs after MAF filtering: {}".format(mt_filt.count_rows()))

    # MHC chr6:25-35Mb
    # chr8.inversion chr8:7-13Mb
    intervals = ['chr6:25M-35M', 'chr8:7M-13M']
    mt_filt = hl.filter_intervals(mt_filt, [hl.parse_locus_interval(x, reference_genome='GRCh38') for x in intervals],
                                  keep=False)
    print("\nNumber of SNPs after MHC and chr8 inversions filtering: {}".format(mt_filt.count_rows()))

    # This step is expensive (on local machine)
    mt_ld_prune = hl.ld_prune(mt_filt.GT, r2=0.8, bp_window_size=500000)
    mt_ld_pruned = mt.filter_rows(hl.is_defined(mt_ld_prune[mt.row_key]))
    print("\nNumber of SNPs after LD pruning: {}".format(mt_ld_pruned.count_rows()))

    return mt_ld_pruned


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
                linewidths=0.5)

    axs.scatter(cat2[x_pc], cat2[y_pc], c=color2, label=label2, s=15, alpha=1, marker=marker2, edgecolor='black',
                linewidths=0.4)

    axs.scatter(cat3[x_pc], cat3[y_pc], c='red', label=label3, s=15, alpha=1, marker='x', edgecolor='black',
                linewidths=0.4)

    plt.legend(bbox_to_anchor=(1, 0.5), loc='lower left', frameon=False, prop={'size': 13})
    axs.set_xlabel(xlabel=x_pc, fontsize=15)
    axs.set_ylabel(ylabel=y_pc, fontsize=15)

    plt.close()

    return fig


def pca_without_ref(
        dirname: str = None,
        basename: str = None,
        input_type: str = None,
        reference: str = 'GRCh38',
        n_pcs: int = 20,
        outdir: str = None):

    print("Reading mt")
    if reference.lower() == 'grch37':
        from gwaspy.utils.reference_liftover import liftover_to_grch38
        mt = liftover_to_grch38(dirname=dirname, basename=basename, input_type=input_type)
    else:
        from gwaspy.utils.read_file import read_infile
        mt = read_infile(input_type=input_type, dirname=dirname, basename=basename)

    print("Filtering mt")
    mt = filter_mt(mt)

    print("\nRunning PCA")
    eigenvalues, pcs, _ = hl.hwe_normalized_pca(mt.GT, k=n_pcs)

    pcs_ht = pcs.transmute(**{f'PC{i}': pcs.scores[i - 1] for i in range(1, n_pcs+1)})

    # add phenotype and sex to the output, using information from the mt
    ann_cols = ['is_case', 'is_female']
    annotations_ht = mt.cols().select(*ann_cols)

    pcs_ht = pcs_ht.annotate(is_case=annotations_ht[pcs_ht.s].is_case)
    pcs_ht = pcs_ht.annotate(is_female=annotations_ht[pcs_ht.s].is_female)

    print("\nSaving PC scores file")
    out_scores_file = outdir + basename + '_scores.tsv'
    pcs_ht.export(out_scores_file)

    print("\nGenerating PCA plots")
    pcs_scores = pd.read_table(out_scores_file, header=0, sep='\t')

    pcs_scores[['is_female']] = pcs_scores[['is_female']].replace([True, False, None], ['female', 'male', 'unknown'])
    pcs_scores[['is_case']] = pcs_scores[['is_case']].replace([True, False, None], ['case', 'control', 'unknown'])

    figs_dict = {}
    for col in ['is_case', 'is_female']:
        for i in range(1, n_pcs, 2):
            xpc = f'PC{i}'
            ypc = f'PC{i + 1}'

            figs_dict["fig{}{}".format(col, i)] = plot_pca(pcs_scores, xpc, ypc, col)

    pdf = PdfPages('{}{}.GWASpy.PCA.plots.pdf'.format(outdir, basename))
    for figname, figure in figs_dict.items():
        pdf.savefig(figure)
    pdf.close()
    # hl.hadoop_copy('file:///tmp/GWASpy.PCA.plots.pdf', outdir)
