__author__ = 'Lindo Nkambule'

import hail as hl
import pandas as pd
from gwaspy.pca.pca_filter_snps import pca_filter_mt, relatedness_check
import plotly.express as px


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


def intersect_ref(
        ref_dirname: str = 'gs://hgdp-1kg/hgdp_tgp/datasets_for_others/lindo/ds_without_outliers/',
        ref_basename: str = 'unrelated',
        data_mt: hl.MatrixTable = None,
        data_basename: str = None, out_dir: str = None):
    """
    Intersects reference panel with the data and writes intersections as matrix tables
    :param ref_dirname: directory name where reference data is
    :param ref_basename: base filename for reference data
    :param data_mt: input data MatrixTable
    :param data_basename: base filename for input data
    :param out_dir: output directory where files are going to be saved to
    :return:
    """
    print('Reading reference data mt')
    ref_mt = hl.read_matrix_table(f'{ref_dirname}{ref_basename}.mt')

    # filter data to sites in ref & array data
    data_in_ref = data_mt.filter_rows(hl.is_defined(ref_mt.rows()[data_mt.row_key]))
    print('\nsites in ref and data, inds in data: {}'.format(data_in_ref.count()))
    data_in_ref.write(f'{out_dir}GWASpy/PCA/{data_basename}/pca_project/{data_basename}_intersect_{ref_basename}.mt',
                      overwrite=True)

    # filter ref to data sites
    ref_in_data = ref_mt.filter_rows(hl.is_defined(data_mt.rows()[ref_mt.row_key]))
    print('\nsites in ref and data, inds in ref: {}'.format(ref_in_data.count()))  #
    ref_in_data.write(f'{out_dir}GWASpy/PCA/{data_basename}/pca_project/{ref_basename}_intersect_{data_basename}.mt',
                      overwrite=True)


def run_ref_pca(
        mt: hl.MatrixTable = None,
        npcs: int = 20):
    """
    Run PCA on a dataset
    :param mt: dataset to run PCA on
    :param npcs: number of principal components to be used in PCA
    :return:
    """
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(mt.GT, k=npcs, compute_loadings=True)
    pca_mt = mt.annotate_rows(pca_af=hl.agg.mean(mt.GT.n_alt_alleles()) / 2)
    pca_loadings = pca_loadings.annotate(pca_af=pca_mt.rows()[pca_loadings.key].pca_af)

    pca_scores = pca_scores.transmute(**{f'PC{i}': pca_scores.scores[i - 1] for i in range(1, npcs+1)})

    return pca_scores, pca_loadings


def plot_pca_ref(data_scores, ref_scores, x_pc, y_pc):
    pcs = pd.read_table(data_scores, header=0, sep='\t')
    pcs['Project'] = "Input Dataset"
    pcs['pop'] = "input dataset"  # use this for plots instead of pop label
    pcs = pcs[['s', 'pop', 'Project', x_pc, y_pc]]

    ref_update = ref_scores.rename(columns={'SuperPop': 'pop'})
    ref_update = ref_update[['s', 'pop', 'Project', x_pc, y_pc]]

    # concatenate the two dfs together
    concatenated = pd.concat([ref_update, pcs], axis=0)

    # https://matplotlib.org/stable/tutorials/colors/colors.html
    color_map = {}
    colors = ['#8C000F', '#00FFFF', '#0343DF', '#653700', '#008000', '#ED0DD9', '#4B0082', '#008080',
              '#FF0000', '#FFD700', '#DDA0DD', '#C0C0C0']

    # get a list of unique population labels in the data
    pops = concatenated['pop'].unique().tolist()
    # update the dictionary with unique colors for each population
    for i in range(len(pops)):
        color_map[pops[i]] = colors[i]

    color_map.update({'input dataset': '#000000'})  # use black for input dataset

    fig = px.scatter(concatenated, x=x_pc, y=y_pc, color='pop',
                     hover_data=['s', x_pc, y_pc, 'pop', 'Project'],
                     color_discrete_map=color_map).update_traces(marker_size=4)

    return fig


def run_pca_project(
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
        run_relatedness_check: bool = True,
        include_kinself: bool = False,
        relatedness_method: str = 'pc_relate',
        relatedness_thresh: float = 0.1,
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
    :param hwe: hardy-weinberg fiter threshold
    :param call_rate: variant call rate filter threshold
    :param ld_cor: reference build
    :param ld_window: window size
    :param run_relatedness_check: whether or not to run relatedness checks
    :param include_kinself: whether or not to include self kinship
    :param relatedness_method: method to use for relatedness filtering
    :param relatedness_thresh: threshold to use for filtering out related individuals
    :param prob_threshold: a list of probability thresholds to use for classifying samples
    :return: a pandas Dataframe with data PCA scores projected on the same PCA space using reference data of choice
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

    print('\nFiltering data mt')
    mt = pca_filter_mt(in_mt=mt, maf=maf, hwe=hwe, call_rate=call_rate, ld_cor=ld_cor, ld_window=ld_window)

    if run_relatedness_check:
        related_out_dir = f'{out_dir}GWASpy/PCA/{data_basename}/pca_project/'
        mt, _ = relatedness_check(in_mt=mt, method=relatedness_method, outdir=related_out_dir,
                                  kin_estimate=relatedness_thresh, include_kinself=include_kinself)
    else:
        print('Skipping relatedness checks')

    # Intersect data with reference
    intersect_ref(ref_dirname=ref_dirname, ref_basename=ref_basename, data_mt=mt, data_basename=data_basename,
                  out_dir=out_dir)

    ref_in_data = hl.read_matrix_table(f'{out_dir}GWASpy/PCA/{data_basename}/pca_project/{ref_basename}_intersect_{data_basename}.mt')

    print('\nComputing reference PCs')
    ref_scores, pca_loadings = run_ref_pca(mt=ref_in_data, npcs=npcs)
    ref_scores = ref_scores.key_by('s')  # make sure we key by s so we can annotate

    # annotate ref info with SuperPop and Project information
    ref_info = hl.import_table(ref_info, key='Sample')
    ref_annotated = ref_scores.annotate(
        SuperPop=ref_info[ref_scores.s].SuperPop,
        Project=ref_info[ref_scores.s].Project)
    ref_df = ref_annotated.to_pandas()

    # project data
    project_mt = hl.read_matrix_table(f'{out_dir}GWASpy/PCA/{data_basename}/pca_project/{data_basename}_intersect_{ref_basename}.mt')

    data_scores = pc_project(mt=project_mt, loadings_ht=pca_loadings)
    data_scores = data_scores.transmute(**{f'PC{i}': data_scores.scores[i - 1] for i in range(1, npcs+1)})
    data_df = data_scores.to_pandas()

    # merge data scores with ref scores
    data_ref = pd.concat([ref_df, data_df], sort=False)

    # export merged data and ref PC scores
    df_to_export = data_ref.copy()
    df_to_export.loc[~df_to_export['SuperPop'].isnull(), 'Projected'] = 'No - reference'
    df_to_export.loc[df_to_export['SuperPop'].isnull(), 'Projected'] = 'Yes - input'
    df_to_export.to_csv(f'{out_dir}GWASpy/PCA/{data_basename}/pca_project/input_and_reference_scores.tsv',
                        sep='\t', index=False)

    from gwaspy.pca.assign_pop_labels import assign_population_pcs
    pcs_df, clf = assign_population_pcs(pop_pc_pd=data_ref, num_pcs=npcs, min_prob=prob_threshold)

    pops = pcs_df['SuperPop'].unique().tolist() # get reference population labels
    pops = [item for item in pops if  not(pd.isnull(item))] # remove NaN values (only reference samples will have the SuperPop label)
    data_pops = pcs_df.loc[pcs_df['SuperPop'].isnull()]
    cols = ['s', 'pop'] + [f'prob_{i}' for i in pops] + [f'PC{i}' for i in range(1, npcs+1)]
    data_pops_df = data_pops[cols]

    data_pops_df.to_csv(f'{out_dir}GWASpy/PCA/{data_basename}/pca_project/pca_sup_pops_{prob_threshold}_probs.project.pca.txt',
                        sep='\t', index=False)

    print("\nGenerating PCA plots")
    data_scores_prob = f'{out_dir}GWASpy/PCA/{data_basename}/pca_project/pca_sup_pops_{prob_threshold}_probs.project.pca.txt'
    ref_scores = df_to_export[df_to_export['Projected'] == 'No - reference']

    figs_dict = {}
    # plotting more than 10 PCA plots in HTML generates wobbly, large files
    for i in range(1, 10, 2):
        xpc = f'PC{i}'
        ypc = f'PC{i + 1}'
        figs_dict["fig{}{}".format(xpc, ypc)] = plot_pca_ref(data_scores=data_scores_prob,
                                                             ref_scores=ref_scores,
                                                             x_pc=xpc, y_pc=ypc)
    with open('/tmp/pca.project.plots.html', 'a') as f:
        for figname, figure in figs_dict.items():
            f.write(figure.to_html(include_plotlyjs='cdn'))

    hl.hadoop_copy('file:///tmp/pca.project.plots.html',
                   f'{out_dir}GWASpy/PCA/{data_basename}/pca_project/{data_basename}.pca.project.plots.html')



