__author__ = 'Lindo Nkambule'

import argparse
import hail as hl


def pca(
        ref_dirname: str = 'gs://hgdp-1kg/hgdp_tgp/datasets_for_others/lindo/ds_without_outliers/',
        ref_basename: str = 'unrelated',
        ref_info: str = 'gs://hgdp-1kg/hgdp_tgp/gwaspy_pca_ref/hgdp_1kg_sample_info.unrelateds.pca_outliers_removed.with_project.tsv',
        reference: str = 'GRCh38', pca_type: str = None,
        data_dirname: str = None, data_basename: str = None, input_type: str = None,
        maf: float = 0.05, hwe: float = 1e-3, call_rate: float = 0.98,
        ld_cor: float = 0.2, ld_window: int = 250000, n_pcs: int = 20, run_relatedness_check: bool = True,
        include_kinself: bool = False, relatedness_method: str = 'pc_relate',
        relatedness_thresh: float = 0.1, prob_threshold: float = 0.8, out_dir: str = None):

    if not out_dir:
        raise Exception('\nOutput directory where files will be saved is not specified')

    if pca_type == 'project':
        print('\nRunning PCA using projection method')

        from gwaspy.pca.pca_project import run_pca_project
        run_pca_project(ref_dirname=ref_dirname, ref_basename=ref_basename, ref_info=ref_info,
                        data_dirname=data_dirname, data_basename=data_basename, out_dir=out_dir, input_type=input_type,
                        reference=reference, npcs=n_pcs, maf=maf, hwe=hwe, call_rate=call_rate,
                        relatedness_method=relatedness_method, run_relatedness_check=run_relatedness_check,
                        ld_cor=ld_cor, ld_window=ld_window, include_kinself=include_kinself,
                        prob_threshold=prob_threshold)

    elif pca_type == 'joint':
        print('\nRunning PCA using joint method')
        from gwaspy.pca.pca_joint import run_pca_joint
        run_pca_joint(ref_dirname=ref_dirname, ref_basename=ref_basename, ref_info=ref_info, data_dirname=data_dirname,
                      data_basename=data_basename, out_dir=out_dir, input_type=input_type, reference=reference,
                      npcs=n_pcs, maf=maf, hwe=hwe, call_rate=call_rate, ld_cor=ld_cor, ld_window=ld_window,
                      relatedness_method=relatedness_method, relatedness_thresh=relatedness_thresh,
                      prob_threshold=prob_threshold)

    else:
        print('\nRunning PCA without a reference')
        from gwaspy.pca.pca_normal import run_pca_normal
        run_pca_normal(dirname=data_dirname, basename=data_basename, input_type=input_type, out_dir=out_dir,
                       reference=reference, maf=maf, hwe=hwe, call_rate=call_rate, ld_cor=ld_cor, ld_window=ld_window,
                       n_pcs=n_pcs, run_relatedness_check=run_relatedness_check, relatedness_method=relatedness_method,
                       relatedness_thresh=relatedness_thresh, include_kinself=include_kinself)


def main():
    parser = argparse.ArgumentParser()
    # reference args
    parser.add_argument('--ref-dirname', default='gs://hgdp-1kg/hgdp_tgp/datasets_for_others/lindo/ds_without_outliers/')
    parser.add_argument('--ref-basename', default='unrelated')
    parser.add_argument('--ref-info', default='gs://hgdp-1kg/hgdp_tgp/gwaspy_pca_ref/hgdp_1kg_sample_info.unrelateds.pca_outliers_removed.with_project.tsv')
    parser.add_argument('--reference', type=str, default='GRCh38')
    parser.add_argument('--pca-type', type=str, default='normal', choices=['normal', 'project', 'joint'])

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
    parser.add_argument('--npcs', type=int, default=20, help='Number of PCs to use')
    parser.add_argument('--no-relatedness', action='store_false')
    parser.add_argument('--include-kinself', action='store_true')
    parser.add_argument('--relatedness-method', type=str, default='pc_relate',
                        choices=['pc_relate', 'ibd', 'king'], help='Method to use for the inference of relatedness')
    parser.add_argument('--relatedness-thresh', type=float, default=0.1,
                        help='Threshold value to use in relatedness checks')
    parser.add_argument('--prob', type=float, default=0.8,
                        help='Minimum probability of belonging to a given population for the population to be set')
    parser.add_argument('--out-dir', type=str, required=True)

    args = parser.parse_args()

    if not args.prob:
        print(f'No prob value specified, {args.prob} will be used')

    hl.init(default_reference=args.reference)

    pca(ref_dirname=args.ref_dirname, ref_basename=args.ref_basename, ref_info=args.ref_info, reference=args.reference,
        pca_type=args.pca_type, input_type=args.input_type, data_dirname=args.data_dirname,
        data_basename=args.data_basename, maf=args.maf, hwe=args.hwe, call_rate=args.geno, ld_cor=args.ld_cor,
        ld_window=args.ld_window, n_pcs=args.npcs, run_relatedness_check=args.no_relatedness,
        include_kinself=args.include_kinself, relatedness_method=args.relatedness_method,
        relatedness_thresh=args.relatedness_thresh, prob_threshold=args.prob, out_dir=args.out_dir)

    print('\nDone running PCA')


if __name__ == '__main__':
    main()

