__author__ = 'Lindo Nkambule'

import argparse
import hailtop.batch as hb
from gwaspy.imputation.impute5_impute import impute5_imputation
from typing import Union


def run_impute(backend: Union[hb.ServiceBackend, hb.LocalBackend] = None,
               input_vcf: str = None,
               vcf_ref: str = None,
               software: str = 'impute5',
               output_filename: str = None,
               out_dir: str = None
               ):

    if software.lower() not in ['beagle5', 'glimpse', 'impute5']:
        raise SystemExit(f'Incorrect software {software} selected. Options are [beagle5, glimpse, impute5]')

    b = hb.Batch(backend=backend,
                 name=f'GWASpy-Imputation-{software.upper()}')

    if vcf_ref == 'hgdp1kgp':
        print(f'\nIMPUTING {input_vcf} WITH HGDP+1KGP PANEL\n')
        ref_path = 'gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/hgdp1kgp_chrCNUMBER.filtered.SNV_INDEL.phased.shapeit5.bcf'
    else:
        print(f'\nIMPUTING {input_vcf} WITH USER-DEFINED REFERENCE PANEL\n')
        ref_path = vcf_ref

    if software == 'impute5':
        impute5_imputation(
            batch=b,
            input_path=input_vcf,
            reference_path=ref_path,
            output_filename=output_filename,
            output_path=out_dir
        )
    # elif: To add BEAGLE
    # else: TO add GLIMPSE (for low coverage data)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-vcf', type=str, required=True)
    parser.add_argument('--vcf-ref', type=str, default=None)
    parser.add_argument('--local', action='store_true')
    parser.add_argument('--billing-project', required=True)
    parser.add_argument('--software', type=str, default='impute5', choices=['beagle5', 'glimpse', 'impute5'])
    parser.add_argument('--output-filename', type=str, required=True)
    parser.add_argument('--out-dir', type=str, required=True)

    args = parser.parse_args()

    if args.local:
        backend = hb.LocalBackend()
    else:
        backend = hb.ServiceBackend(billing_project=args.billing_project,
                                    remote_tmpdir=f'{args.out_dir}/tmp/')

    run_impute(backend=backend,
               input_vcf=args.input_vcf,
               vcf_ref=args.vcf_ref,
               software=args.software,
               output_filename=args.output_filename,
               out_dir=args.out_dir)
