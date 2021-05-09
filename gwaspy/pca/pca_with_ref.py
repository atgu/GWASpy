__author__ = 'Lindo Nkambule'

import hail as hl
import pandas as pd
import random
from sklearn.ensemble import RandomForestClassifier
from typing import Tuple


def pc_project(
        mt: hl.MatrixTable,
        loadings_ht: hl.Table,
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


def pca_with_ref(
        dirname: str = None,
        basename: str = None,
        pca_loadings: str = 'gs://covid19-hg-public/pca_projection/hgdp_tgp_pca_covid19hgi_snps_loadings.ht',
        outdir: str = None,
        input_type: str = None,
        reference: str = 'GRCh38') -> pd.DataFrame:
    """
    Project samples into predefined PCA space
    :param dirname: matrix table of data to project
    :param basename: matrix table of data to project
    :param pca_loadings: existing PCA space
    :param outdir: directory and filename prefix for where to put PCA projection output
    :param input_type: input file(s) type: hail, plink, or vcf
    :param reference: reference build
    :return: a pandas Dataframe with data PCA scores projected on the same PCA space using the Human Genome Diversity
    Project(HGDP) and the 1000 Genomes Project samples as reference
    """

    SAMPLE_FIELD_NAME = "s"

    print('Reading mt')
    if reference.lower() == 'grch37':
        from gwaspy.utils.reference_liftover import liftover_to_grch38
        mt = liftover_to_grch38(dirname=dirname, basename=basename, input_type=input_type)
    else:
        from gwaspy.utils.read_file import read_infile
        mt = read_infile(input_type=input_type, dirname=dirname, basename=basename)

    print('\nReading loadings')
    loadings = hl.read_table(pca_loadings)
    print('There are {} SNPs in the loading file'.format(loadings.count()))

    print('\nThe data has {} SNPs'.format(mt.count_rows()))
    mt = mt.filter_rows(hl.is_defined(loadings[mt.locus, mt.alleles]))
    print('{} SNPs will be used for projections'.format(mt.count_rows()))

    print('\nProjecting data')
    ht_projections = pc_project(mt, loadings)
    ht_projections = ht_projections.transmute(**{f"PC{i}": ht_projections.scores[i - 1] for i in range(1, 21)})

    print('\nWriting projection output')
    if outdir:
        # if specified, write out the pca score to the out directory
        ht_projections.write('{}{}_pca_project_scores.ht'.format(outdir, basename))
        ht = hl.read_table('{}{}_pca_project_scores.ht'.format(outdir, basename))
    else:
        # write out the PCA scores to the same directory where input files are
        ht_projections.write('{}{}_pca_project_scores.ht'.format(dirname, basename))
        ht = hl.read_table('{}{}_pca_project_scores.ht'.format(dirname, basename))

    df = ht.to_pandas()

    return df


def merge_data_with_ref(
        refscores: str = 'gs://covid19-hg-public/pca_projection/hgdp_tgp_pca_covid19hgi_snps_scores.txt.gz',
        ref_info: str = 'gs://covid19-hg-public/pca_projection/gnomad_meta_hgdp_tgp_v1.txt',
        data_scores: pd.DataFrame = None) -> pd.DataFrame:
    """
    Merge data with ref
    :param refscores: path to reference score
    :param ref_info: path to information about samples in the ref scores
    :param data_scores: pandas DataFrame with the data scores
    :return: a pandas Dataframe of data merged with reference
    """

    print('\nMerging data with ref')
    ref = pd.read_table(refscores, header=0, sep='\t', compression='gzip')
    ref_info = pd.read_table(ref_info, header=0, sep='\t', low_memory=False)
    ref_info = ref_info[['project_meta.sample_id', 'hgdp_tgp_meta.Population', 'hgdp_tgp_meta.Genetic.region']]

    # rename columns in ref_info
    ref_info.columns = ['s', 'Population', 'SuperPop']

    ref_merge = pd.merge(left=ref, right=ref_info, left_on='s', right_on='s', how='inner')

    data_ref = pd.concat([ref_merge, data_scores], sort=False)
    print('Done merging data with ref')

    return data_ref


def assign_population_pcs(
        pop_pc_pd: pd.DataFrame,
        num_pcs: int,
        known_col: str = 'SuperPop',
        fit: RandomForestClassifier = None,
        seed: int = 42,
        prop_train: float = 0.8,
        n_estimators: int = 100,
        min_prob: float = 0.9,
        output_col: str = 'pop',
        missing_label: str = 'oth'
) -> Tuple[pd.DataFrame, RandomForestClassifier]:
    """
    This function uses a random forest model to assign population labels based on the results of PCA.
    Default values for model and assignment parameters are those used in gnomAD.
    :param Table pop_pc_pd: Pandas dataframe containing population PCs as well as a column with population labels
    :param str known_col: Column storing the known population labels
    :param str pcs_col: Columns storing the PCs
    :param RandomForestClassifier fit: fit from a previously trained random forest model (i.e., the output from a previous RandomForestClassifier() call)
    :param int num_pcs: number of population PCs on which to train the model
    :param int seed: Random seed
    :param float prop_train: Proportion of known data used for training
    :param int n_estimators: Number of trees to use in the RF model
    :param float min_prob: Minimum probability of belonging to a given population for the population to be set (otherwise set to `None`)
    :param str output_col: Output column storing the assigned population
    :param str missing_label: Label for samples for which the assignment probability is smaller than `min_prob`
    :return: Dataframe containing sample IDs and imputed population labels, trained random forest model
    :rtype: DataFrame, RandomForestClassifier
    """

    # Expand PC column
    pc_cols = ['PC{}'.format(i + 1) for i in range(num_pcs)]
    train_data = pop_pc_pd.loc[~pop_pc_pd[known_col].isnull()]

    N = len(train_data)

    # Split training data into subsamples for fitting and evaluating
    if not fit:
        random.seed(seed)
        train_subsample_ridx = random.sample(list(range(0, N)), int(N * prop_train))
        train_fit = train_data.iloc[train_subsample_ridx]
        fit_samples = [x for x in train_fit['s']]
        evaluate_fit = train_data.loc[~train_data['s'].isin(fit_samples)]

        # Train RF
        training_set_known_labels = train_fit[known_col].values
        training_set_pcs = train_fit[pc_cols].values
        evaluation_set_pcs = evaluate_fit[pc_cols].values

        pop_clf = RandomForestClassifier(n_estimators=n_estimators, random_state=seed)
        pop_clf.fit(training_set_pcs, training_set_known_labels)
        print('Random forest feature importances are as follows: {}'.format(pop_clf.feature_importances_))

        # Evaluate RF
        predictions = pop_clf.predict(evaluation_set_pcs)
        error_rate = 1 - sum(evaluate_fit[known_col] == predictions) / float(len(predictions))
        print('Estimated error rate for RF model is {}'.format(error_rate))
    else:
        pop_clf = fit

    # Classify data
    print('Classifying data')
    pop_pc_pd[output_col] = pop_clf.predict(pop_pc_pd[pc_cols].values)
    probs = pop_clf.predict_proba(pop_pc_pd[pc_cols].values)
    probs = pd.DataFrame(probs, columns=[f'prob_{p}' for p in pop_clf.classes_])

    pop_pc_pd = pd.concat([pop_pc_pd.reset_index(drop=True), probs.reset_index(drop=True)], axis=1)

    probs['max'] = probs.max(axis=1)
    pop_pc_pd.loc[probs['max'] < min_prob, output_col] = missing_label

    return pop_pc_pd, pop_clf


# pca_with_ref(dirname="/Users/lindokuhle/Desktop/preimp_qc/data/", basename="sim_sim2a_eur_sa_merge.miss",
             # reference="GRCh37")

