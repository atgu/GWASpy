__author__ = 'Lindo Nkambule'

import hail as hl
import pandas as pd
from gwaspy.pca.pca_filter_snps import pca_filter_mt, relatedness_check
import plotly.express as px


def joint_pca(
        ref_dirname: str = 'gs://hgdp-1kg/hgdp_tgp/datasets_for_others/lindo/ds_without_outliers/',
        ref_basename: str = 'unrelated',
        in_mt: hl.MatrixTable = None,
        data_basename: str = None,
        npcs: int = 20,
        out_dir: str = None):
    """
    Merges input dataset with ref by [locus, alleles] and runs PCA on merged dataset
    :param ref_dirname: directory name where reference data is
    :param ref_basename: base filename for reference data
    :param in_mt: input data MatrixTable
    :param data_basename: base filename for input data
    :param npcs: number of principal components to be used in PCA
    :param out_dir: output directory where files are going to be saved to
    :return:
    """
    print('\nReading reference data mt')
    ref_mt = hl.read_matrix_table(f'{ref_dirname}{ref_basename}.mt')

    # We need to unkey the datasets and take only cols common between the two in order to be able to merge in Hail
    ref_mt = ref_mt.key_cols_by().key_rows_by()
    ref_downsampled = ref_mt.select_cols('s').select_rows('locus', 'alleles').select_entries('GT')
    ref_downsampled = ref_downsampled.key_cols_by('s').key_rows_by('locus', 'alleles')

    data_mt = in_mt.key_cols_by().key_rows_by()
    data_downsampled = data_mt.select_cols('s').select_rows('locus', 'alleles').select_entries('GT')
    data_downsampled = data_downsampled.key_cols_by('s').key_rows_by('locus', 'alleles')

    print('\nJoining Data with Ref by locus & alleles')
    joined = ref_downsampled.union_cols(data_downsampled)

    pca_snps = joined.count_rows()
    if pca_snps > 1000000:
        import warnings
        warnings.warn(f'Too many SNPs to be used in PCA: {pca_snps}. This will make PCA run longer')

    print(f'\nRunning PCA with {npcs} principal components')
    pca_evals, pca_scores, _ = hl.hwe_normalized_pca(joined.GT, k=npcs)

    pca_scores = pca_scores.transmute(**{f'PC{i}': pca_scores.scores[i - 1] for i in range(1, npcs+1)})
    print(f'\nExporting PCA scores to {out_dir}')
    pca_scores.export(f'{out_dir}GWASpy/PCA/{data_basename}/pca_joint/{data_basename}.1kg_hgdp.joint.pca.scores.txt.bgz')


def add_ref_superpop_labels(joint_scores: str = None, ref_info: str = None) -> pd.DataFrame:
    """
    Add SuperPop labels to reference samples so we can use them in the RF model to assign ancestry
    :param joint_scores: path to joint data+ref scores
    :param ref_info: path to information about samples in the ref scores
    :return: a pandas Dataframe of joint data+ref scores with SuperPop label for ref samples
    """
    print('\nAdding SuperPop labels to reference samples')
    joint_data = pd.read_table(joint_scores, header=0, sep='\t', compression='gzip')
    ref_info = pd.read_table(ref_info, header=0, sep='\t')
    updated_joint_data = pd.merge(left=joint_data, right=ref_info, left_on='s', right_on='Sample', how='left')

    return updated_joint_data


def plot_pca_joint(joint_scores: pd.DataFrame = None, x_pc: str = None, y_pc: str = None):
    """
    Plot PCs
    :param joint_scores: pandas dataframe of joint data+ref scores
    :param x_pc: x-axis pc
    :param y_pc: y-axis pc
    :return: pc plot
    """
    # split data and ref into separate dfs
    data = joint_scores.loc[joint_scores['SuperPop'].isnull()]
    ref = joint_scores.loc[joint_scores['SuperPop'].notnull()]

    data['Project'] = "Input Dataset"
    data = data[['s', 'pop', 'Project', x_pc, y_pc]]
    ref = ref[['s', 'SuperPop', 'Project', x_pc, y_pc]]
    ref.rename(columns={'SuperPop': 'pop'}, inplace=True)

    # concatenate the two dfs together
    concatenated = pd.concat([ref, data], axis=0)

    fig = px.scatter(concatenated, x=x_pc, y=y_pc, color='pop',
                     hover_data=['s', x_pc, y_pc, 'pop', 'Project'],
                     color_discrete_map={'AFR': "#984EA3", 'EAS': "#4DAF4A", 'EUR': "#377EB8", 'CSA': "#FF7F00",
                                         'AMR': "#E41A1C", 'MID': "#A65628", 'OCE': "#999999",
                                         'oth': "#F0E442"}).update_traces(marker_size=4)

    return fig


def run_pca_joint(
        ref_dirname: str = 'gs://hgdp-1kg/hgdp_tgp/datasets_for_others/lindo/ds_without_outliers/',
        ref_basename: str = 'unrelated',
        ref_info: str = 'gs://hgdp-1kg/hgdp_tgp/gwaspy_pca_ref/hgdp_1kg_sample_info.unrelateds.pca_outliers_removed.with_project.tsv',
        data_dirname: str = None,
        data_basename: str = None,
        out_dir: str = None,
        input_type: str = None,
        reference: str = 'GRCh38',
        npcs: int = 20,
        maf: float = 0.05,
        hwe: float = 1e-3,
        call_rate: float = 0.98,
        ld_cor: float = 0.2,
        ld_window: int = 250000,
        relatedness_method: str = 'pc_relate',
        relatedness_thresh: float = 0.98,
        prob_threshold: float = 0.8):
    """
    Project samples into predefined PCA space
    :param ref_dirname: directory name where reference data is
    :param ref_basename: base filename for reference data
    :param ref_info: reference sample information
    :param data_dirname: matrix table of data to project
    :param data_basename: matrix table of data to project
    :param out_dir: directory and filename prefix for where to put PCA projection output
    :param input_type: input file(s) type: hail, plink, or vcf
    :param reference: reference build
    :param npcs: number of principal components to be used in PCA
    :param maf: minor allele frequency threshold
    :param hwe: hardy-weinberg filter threshold
    :param call_rate: variant call rate filter threshold
    :param ld_cor: reference build
    :param ld_window: window size
    :param prob_threshold: a list of probability thresholds to use for classifying samples
    :param relatedness_method: method to use for relatedness filtering
    :param relatedness_thresh: threshold to use for filtering out related individuals
    :return: a pandas Dataframe with data PCA scores projected on the same PCA space using the Human Genome Diversity
    """
    print('\nReading data mt')
    if reference.lower() == 'grch37':
        lifted_over = f'{data_dirname}{data_basename}.liftover.grch38.mt'
        if not hl.hadoop_exists(lifted_over):
            from gwaspy.utils.reference_liftover import liftover_to_grch38
            mt = liftover_to_grch38(dirname=data_dirname, basename=data_basename, input_type=input_type)
        else:
            print(f'\nFound lifted-over over file: {lifted_over}')
            mt = hl.read_matrix_table(lifted_over)
    else:
        from gwaspy.utils.read_file import read_infile
        mt = read_infile(input_type=input_type, dirname=data_dirname, basename=data_basename)

    print("\nFiltering data mt")
    data_mt = pca_filter_mt(in_mt=mt, maf=maf, hwe=hwe, call_rate=call_rate, ld_cor=ld_cor, ld_window=ld_window)

    data_mt = relatedness_check(in_mt=data_mt, method=relatedness_method, outdir=out_dir,
                                kin_estimate=relatedness_thresh)

    joint_pca(ref_dirname=ref_dirname, ref_basename=ref_basename, in_mt=data_mt, data_basename=data_basename, npcs=npcs,
              out_dir=out_dir)

    scores_without_pop_label = f'{out_dir}GWASpy/PCA/{data_basename}/pca_joint/{data_basename}.1kg_hgdp.joint.pca.scores.txt.bgz'
    scores_with_pop_label_df = add_ref_superpop_labels(joint_scores=scores_without_pop_label, ref_info=ref_info)

    from gwaspy.pca.assign_pop_labels import assign_population_pcs
    pcs_df, clf = assign_population_pcs(pop_pc_pd=scores_with_pop_label_df, num_pcs=npcs, min_prob=prob_threshold)

    data_pops = pcs_df.loc[pcs_df['SuperPop'].isnull()]
    data_pops['pop'].value_counts()
    cols = ['s', 'pop'] + [f'prob_{i}' for i in ["AFR", "AMR", "CSA", "EAS", "EUR", "MID", "OCE"]] + [f'PC{i}' for i in
                                                                                                      range(1, npcs+1)]
    data_pops_df = data_pops[cols]

    data_pops_df.to_csv(f'{out_dir}GWASpy/PCA/{data_basename}/pca_joint/pca_sup_pops_{prob_threshold}_probs.joint.pca.txt',
                        sep='\t', index=False)

    print('\nGenerating PCA plots')

    figs_dict = {}
    # plotting more than 10 PCA plots in HTML generates wobbly, large files
    for i in range(1, 10, 2):
        xpc = f'PC{i}'
        ypc = f'PC{i + 1}'
        figs_dict["fig{}{}".format(xpc, ypc)] = plot_pca_joint(joint_scores=pcs_df,
                                                               x_pc=xpc, y_pc=ypc)

    with open('/tmp/joint.pca.plots.html', 'a') as f:
        for figname, figure in figs_dict.items():
            f.write(figure.to_html(include_plotlyjs='cdn'))

    hl.hadoop_copy('file:///tmp/joint.pca.plots.html',
                   f'{out_dir}GWASpy/PCA/{data_basename}/pca_joint/{data_basename}.joint.pca.plots.html')

