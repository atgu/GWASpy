__author__ = 'Lindo Nkambule'

import argparse
import hail as hl


def pca(
        ref_dirname: str = 'gs://african-seq-data/hgdp_tgp/gwaspy_pca_ref/',
        ref_basename: str = 'hgdp_1kg_filtered_maf_5_GRCh38',
        ref_info: str = 'gs://african-seq-data/hgdp_tgp/gwaspy_pca_ref/hgdp_1kg_sample_info.tsv',
        reference: str = 'GRCh38', with_ref: bool = False,
        data_dirname: str = None, data_basename: str = None, input_type: str = None,
        maf: float = 0.05, hwe: float = 1e-3, call_rate: float = 0.98,
        ld_cor: float = 0.2, ld_window: int = 250000, n_pcs: int = 20, relatedness_method: str = 'pc_relate',
        relatedness_thresh: float = 0.98, prob_threshold: float = 0.8, out_dir: str = None):

    if not out_dir:
        raise Exception("Output directory where files will be saved is not specified")

    if with_ref:
        print("Running PCA with a reference")

        from gwaspy.pca.pca_with_ref import pca_with_ref
        pca_with_ref(ref_dirname=ref_dirname, ref_basename=ref_basename, ref_info=ref_info, data_dirname=data_dirname,
                     data_basename=data_basename,out_dir=out_dir, input_type=input_type, reference=reference,
                     maf=maf, hwe=hwe, call_rate=call_rate, ld_cor=ld_cor, ld_window=ld_window,
                     prob_threshold=prob_threshold)

    else:
        print("Running PCA without a reference")
        from gwaspy.pca.pca_no_ref import pca_without_ref
        pca_without_ref(dirname=data_dirname, basename=data_basename, input_type=input_type, out_dir=out_dir, reference=reference,
                        maf=maf, hwe=hwe, call_rate=call_rate, ld_cor=ld_cor, ld_window=ld_window, n_pcs=n_pcs,
                        relatedness_method=relatedness_method, relatedness_thresh=relatedness_thresh)


def main():
    parser = argparse.ArgumentParser()
    # reference args
    parser.add_argument('--ref-dirname', default='gs://african-seq-data/hgdp_tgp/gwaspy_pca_ref/')
    parser.add_argument('--ref-basename', default='hgdp_1kg_filtered_maf_5_GRCh38')
    parser.add_argument('--ref-info', default='gs://african-seq-data/hgdp_tgp/gwaspy_pca_ref/hgdp_1kg_sample_info.tsv')
    parser.add_argument('--reference', type=str, default='GRCh38')
    parser.add_argument('--with-ref', action='store_true')

    # data args
    parser.add_argument('--data-dirname', type=str, required=True)
    parser.add_argument('--data-basename', type=str, required=True)
    parser.add_argument('--input-type', type=str, required=True, choices=['vcf', 'plink', 'hail'])

    # filter args
    parser.add_argument('--maf', type=float, default=0.05, help='include only SNPs with MAF >= NUM in PCA')
    parser.add_argument('--hwe', type=float, default=1e-3, help='include only SNPs with HWE >= NUM in PCA')
    parser.add_argument('--geno', type=float, default=0.98, help='include only SNPs with call-rate > NUM')
    parser.add_argument('--ld-cor', type=float, default=0.2, choices=range(0,1), metavar="[0.0-1.0]",
                        help='Squared correlation threshold (exclusive upper bound). Must be in the range [0.0, 1.0]')
    parser.add_argument('--ld-window', type=int, default=250000,
                        help='Window size in base pairs (inclusive upper bound)')
    parser.add_argument('--relatedness-method', type=str, default='pc_relate',
                        choices=['pc_relate', 'ibd', 'king'], help='Method to use for the inference of relatedness')
    parser.add_argument('--relatedness-thresh', type=float, default=0.98,
                        help='Threshold value to use in relatedness checks')
    parser.add_argument('--prob', type=float, default=0.8,
                        help='Minimum probability of belonging to a given population for the population to be set')
    parser.add_argument('--out-dir', type=str, required=True)

    args = parser.parse_args()

    if not args.prob:
        print("No prob value specified, {} will be used".format(args.prob))

    hl.init(default_reference=args.reference)

    pca(ref_dirname=args.ref_dirname, ref_basename=args.ref_basename, ref_info=args.ref_info, reference=args.reference,
        with_ref=args.with_ref,
        input_type=args.input_type, data_dirname=args.data_dirname, data_basename=args.data_basename,
        maf=args.maf, hwe=args.hwe, call_rate=args.geno, ld_cor=args.ld_cor, ld_window=args.ld_window,
        relatedness_method=args.relatedness_method, relatedness_thresh=args.relatedness_thresh,
        prob_threshold=args.prob, out_dir=args.out_dir)

    print("Done running PCA")


if __name__ == '__main__':
    main()

