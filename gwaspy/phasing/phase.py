__author__ = 'Lindo Nkambule'

import argparse
import hailtop.batch as hb
from gwaspy.phasing.shapeit5_phase import shapeit_phasing
from typing import Union


def run_phase(backend: Union[hb.ServiceBackend, hb.LocalBackend] = None,
              input_vcf: str = None,
              vcf_ref: str = None,
              fam_file: str = None,
              data_type: str = 'array',
              software: str = 'shapeit',
              output_filename: str = None,
              out_dir: str = None):

    if data_type.lower() not in ['array', 'wgs']:
        raise SystemExit(f'Incorrect data type {data_type} selected. Options are [array, wgs]')

    if software.lower() not in ['beagle', 'shapeit']:
        raise SystemExit(f'Incorrect software {software} selected. Options are [beagle, shapeit]')

    b = hb.Batch(backend=backend,
                 name=f'GWASpy-Phasing-{software.upper()}')

    if vcf_ref:
        if vcf_ref == 'hgdp1kgp':
            print(f'\nPHASING {input_vcf} WITH HGDP+1KGP PANEL\n')
            ref_path = 'gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/hgdp1kgp_chrCNUMBER.filtered.SNV_INDEL.phased.shapeit5.bcf'
        else:
            print(f'\nPHASING {input_vcf} WITH USER-DEFINED REFERENCE PANEL\n')
            ref_path = vcf_ref
    else:
        ref_path = None
        print(f'\nPHASING {input_vcf} WITHOUT A REFERENCE PANEL\n')

    pedigree = b.read_input(fam_file) if fam_file else None

    if software == 'shapeit':
        shapeit_phasing(
            batch=b,
            input_path=input_vcf,
            reference_path=ref_path,
            fam_file=pedigree,
            data_type=data_type,
            output_filename=output_filename,
            output_path=out_dir)
    # else: To add BEAGLE


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-vcf', type=str, required=True)
    parser.add_argument('--vcf-ref', type=str, default=None)
    parser.add_argument('--pedigree', type=str, default=None)
    parser.add_argument('--local', action='store_true')
    parser.add_argument('--billing-project', required=True)
    parser.add_argument('--software', type=str, default='shapeit', choices=['beagle', 'shapeit'])
    parser.add_argument('--output-filename', type=str, required=True)
    parser.add_argument('--out-dir', type=str, required=True)

    args = parser.parse_args()

    if args.local:
        backend = hb.LocalBackend()
    else:
        backend = hb.ServiceBackend(billing_project=args.billing_project,
                                    remote_tmpdir=f'{args.out_dir}/tmp/')

    run_phase(backend=backend,
              input_vcf=args.input_vcf,
              vcf_ref=args.vcf_ref,
              fam_file=args.pedigree,
              software=args.software,
              output_filename=args.output_filename,
              out_dir=args.out_dir)
