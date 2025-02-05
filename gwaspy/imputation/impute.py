__author__ = 'Lindo Nkambule'

import argparse
import hailtop.batch as hb
from gwaspy.imputation.impute5_impute import impute5_imputation
from gwaspy.imputation.glimpse2_impute import glimpse_phase_impute
from typing import Union


def run_impute(backend: Union[hb.ServiceBackend, hb.LocalBackend] = None,
               input_file: str = None,
               vcf_ref: str = None,
               chromosomes: str = "all",
               software: str = 'impute5',
               output_filename: str = None,
               n_samples: int = None,
               n_panel_samples: int = 4091,
               out_dir: str = None
               ):

    if software.lower() not in ['beagle5', 'glimpse2', 'impute5']:
        raise SystemExit(f'Incorrect software {software} selected. Options are [beagle5, glimpse2, impute5]')

    b = hb.Batch(backend=backend,
                 name=f'GWASpy-Imputation-{software.upper()}')

    if vcf_ref == 'hgdp1kgp':
        print(f'\nIMPUTING GENOTYPES WITH HGDP+1KGP PANEL\n')
        ref_path = 'gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/hgdp1kgp_chrCNUMBER.filtered.SNV_INDEL.phased.shapeit5.bcf'
    else:
        print(f'\nIMPUTING GENOTYPES WITH USER-DEFINED REFERENCE PANEL\n')
        ref_path = vcf_ref

    if software == 'impute5':
        print(f'\nIMPUTING GENOTYPES USING IMPUTE5\n')
        impute5_imputation(
            batch=b,
            input_path=input_file,
            reference_path=ref_path,
            chromosomes=chromosomes,
            output_filename=output_filename,
            n_samples=n_samples,
            n_panel_samples=n_panel_samples,
            output_path=out_dir
        )
    elif software == 'glimpse2':
        glimpse_phase_impute(
            batch=b,
            bam_files=input_file,
            reference_path=ref_path,
            chromosomes = chromosomes,
            output_filename=output_filename,
            output_path=out_dir
        )
    # else: TO add BEAGLE


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-file', type=str, required=True)
    parser.add_argument('--vcf-ref', type=str, default='hgdp1kgp')
    parser.add_argument('--chromosomes', type=str, default='all')
    parser.add_argument('--local', action='store_true')
    parser.add_argument('--billing-project', required=True)
    parser.add_argument('--n-samples', type=int, required=True)
    parser.add_argument('--n-ref-samples', type=int, default=4091)
    parser.add_argument('--software', type=str, default='impute5', choices=['beagle5', 'glimpse2', 'impute5'])
    parser.add_argument('--output-filename', type=str, required=True)
    parser.add_argument('--out-dir', type=str, required=True)

    args = parser.parse_args()

    if args.local:
        backend = hb.LocalBackend()
    else:
        backend = hb.ServiceBackend(billing_project=args.billing_project,
                                    remote_tmpdir=f'{args.out_dir}/tmp/')

    run_impute(backend=backend,
               input_file=args.input_file,
               vcf_ref=args.vcf_ref,
               chromosomes=args.chromosomes,
               software=args.software,
               output_filename=args.output_filename,
               n_samples=args.n_samples,
               n_panel_samples=args.n_ref_samples,
               out_dir=args.out_dir)

    backend.close()
