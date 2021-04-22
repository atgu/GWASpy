__author__ = 'Lindo Nkambule'

import argparse
from typing import List


def pca(input_type: str = None, dirname: str = None, basename: str = None,
        ref_scores: str = 'gs://covid19-hg-public/pca_projection/hgdp_tgp_pca_covid19hgi_snps_scores.txt.gz',
        ref_info: str = 'gs://covid19-hg-public/pca_projection/gnomad_meta_hgdp_tgp_v1.txt', with_ref: bool = True,
        prob: List = None, reference: str = None, out_dir: str = None):

    if not out_dir:
        raise Exception("Output directory where files will be saved is not specified")

    if with_ref:
        print("Running PCA with a reference")

        from gwaspy.pca.pca_with_ref import pca_with_ref, merge_data_with_ref, assign_population_pcs

        data_pcs_df = pca_with_ref(dirname=dirname, basename=basename, outdir=out_dir,
                                   reference=reference, input_type=input_type)

        # data merged with ref
        data_merged_ref = merge_data_with_ref(refscores=ref_scores, ref_info=ref_info, data_scores=data_pcs_df)

        for threshold in prob:
            pcs_df, clf = assign_population_pcs(pop_pc_pd=data_merged_ref, num_pcs=20, min_prob=threshold)

            data_pops = pcs_df.loc[pcs_df['SuperPop'].isnull()]
            data_pops['pop'].value_counts()
            superpops = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID', 'OCE']
            cols = ['s', 'pop'] + [f'prob_{i}' for i in superpops] + [f'PC{i}' for i in range(1, 21)]
            data_pops_df = data_pops[cols]

            data_pops_df.to_csv('{}pca_sup_pops_{}_probs.txt'.format(out_dir, threshold),
                                sep='\t', index=False)

    else:
        print("Running PCA without a reference")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dirname', type=str, required=True)
    parser.add_argument('--basename', type=str, required=True)
    parser.add_argument('--input-type', type=str, required=True, choices=['vcf', 'plink', 'hail'])
    parser.add_argument('--ref-scores', type=str,
                        default='gs://covid19-hg-public/pca_projection/hgdp_tgp_pca_covid19hgi_snps_scores.txt.gz')
    parser.add_argument('--ref-info', type=str,
                        default='gs://covid19-hg-public/pca_projection/gnomad_meta_hgdp_tgp_v1.txt')
    parser.add_argument('--with-ref', action='store_true')
    parser.add_argument('--prob', action='append',
                        help="Minimum probability of belonging to a given population for the population to be set")
    parser.add_argument('--reference', type=str, default='grch38')
    parser.add_argument('--out-dir', type=str, required=True)

    args = parser.parse_args()

    pca(input_type=args.input_type, dirname=args.dirname, basename=args.basename, ref_scores=args.ref_scores,
        ref_info=args.ref_info, with_ref=args.with_ref, prob=args.prob, reference=args.reference, out_dir=args.out_dir)


if __name__ == '__main__':
    main()

