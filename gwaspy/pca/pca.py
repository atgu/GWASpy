__author__ = 'Lindo Nkambule'

import argparse
import hail as hl


def pca(input_type: str = None, dirname: str = None, basename: str = None,
        ref_scores: str = 'gs://covid19-hg-public/pca_projection/hgdp_tgp_pca_covid19hgi_snps_scores.txt.gz',
        ref_info: str = 'gs://covid19-hg-public/pca_projection/gnomad_meta_hgdp_tgp_v1.txt', with_ref: bool = False,
        prob: float = 0.8, reference: str = None, maf: float = 0.05, hwe: float = 1e-3, call_rate: float = 0.98,
        ld_cor: float = 0.2, ld_window: int = 250000, n_pcs: int = 20, relatedness_method: str = 'pc_relate',
        relatedness_thresh: float = 0.98, out_dir: str = None):

    if not out_dir:
        raise Exception("Output directory where files will be saved is not specified")

    if with_ref:
        print("Running PCA with a reference")

        from gwaspy.pca.pca_with_ref import pca_with_ref, merge_data_with_ref, assign_population_pcs

        data_pcs_df = pca_with_ref(dirname=dirname, basename=basename, outdir=out_dir, reference=reference,
                                   input_type=input_type, maf=maf, hwe=hwe, call_rate=call_rate, ld_cor=ld_cor,
                                   ld_window=ld_window)

        # data merged with ref
        data_merged_ref = merge_data_with_ref(refscores=ref_scores, ref_info=ref_info, data_scores=data_pcs_df)

        pcs_df, clf = assign_population_pcs(pop_pc_pd=data_merged_ref, num_pcs=20, min_prob=prob)

        data_pops = pcs_df.loc[pcs_df['SuperPop'].isnull()]
        data_pops['pop'].value_counts()
        superpops = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID', 'OCE']
        cols = ['s', 'pop'] + [f'prob_{i}' for i in superpops] + [f'PC{i}' for i in range(1, 21)]
        data_pops_df = data_pops[cols]

        data_pops_df.to_csv('{}{}_pca_sup_pops_{}_probs.txt'.format(out_dir, basename, prob),
                            sep='\t', index=False)

    else:
        print("Running PCA without a reference")
        from gwaspy.pca.pca_no_ref import pca_without_ref
        pca_without_ref(dirname=dirname, basename=basename, input_type=input_type, outdir=out_dir, reference=reference,
                        maf=maf, hwe=hwe, call_rate=call_rate, ld_cor=ld_cor, ld_window=ld_window, n_pcs=n_pcs,
                        relatedness_method=relatedness_method, relatedness_thresh=relatedness_thresh)


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
    parser.add_argument('--maf', type=float, default=0.05, help='include only SNPs with MAF >= NUM in PCA')
    parser.add_argument('--hwe', type=float, default=1e-3, help='include only SNPs with HWE >= NUM in PCA')
    parser.add_argument('--geno', type=float, default=0.98, help='include only SNPs with call-rate > NUM')
    parser.add_argument('--ld-cor', type=float, default=0.2, choices=range(0,1), metavar="[0.0-1.0]",
                        help="Squared correlation threshold (exclusive upper bound). Must be in the range [0.0, 1.0]")
    parser.add_argument('--ld-window', type=int, default=250000,
                        help="Window size in base pairs (inclusive upper bound)")
    parser.add_argument('--relatedness-method', type=str, default='pc_relate',
                        choices=['pc_relate', 'ibd', 'king'], help='Method to use for the inference of relatedness')
    parser.add_argument('--relatedness-thresh', type=float, default=0.98,
                        help='Threshold value to use in relatedness checks')
    parser.add_argument('--prob', type=float, default=0.8,
                        help="Minimum probability of belonging to a given population for the population to be set")
    parser.add_argument('--reference', type=str, default='GRCh38')
    parser.add_argument('--out-dir', type=str, required=True)

    args = parser.parse_args()

    if not args.prob:
        print("No prob value specified, {} will be used".format(args.prob))

    hl.init(default_reference=args.reference)

    pca(input_type=args.input_type, dirname=args.dirname, basename=args.basename, ref_scores=args.ref_scores,
        ref_info=args.ref_info, with_ref=args.with_ref, prob=args.prob, reference=args.reference, maf=args.maf,
        hwe=args.hwe, call_rate=args.geno, ld_cor=args.ld_cor, ld_window=args.ld_window,
        relatedness_method=args.relatedness_method, relatedness_thresh=args.relatedness_thresh, out_dir=args.out_dir)

    print("Done running PCA")


if __name__ == '__main__':
    main()

