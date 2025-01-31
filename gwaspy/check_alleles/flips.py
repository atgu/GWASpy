__author__ = 'Lindo Nkambule'

import argparse
import hailtop.batch as hb
from gwaspy.check_alleles.check_alleles import check_alleles_workflow
from typing import Union


def run_checks_fix(
        backend: Union[hb.ServiceBackend, hb.LocalBackend] = None,
        input_vcf: str = None,
        ref_path: str = None,
        step: str = "check",
        fix_mode: str = "top",
        output_filename: str = None,
        out_dir: str = None
):
    b = hb.Batch(backend=backend,
                 name=f'GWASpy-{step.upper()}-Alleles')

    check_alleles_workflow(
        batch=b,
        input_path=input_vcf,
        reference_path=ref_path,
        output_filename=output_filename,
        step=step,
        fix_mode=fix_mode,
        output_path=out_dir
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-vcf', type=str, required=True)
    parser.add_argument('--ref-fasta', type=str, default='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta')
    parser.add_argument('--local', action='store_true')
    parser.add_argument('--billing-project', required=True)
    parser.add_argument('--step', type=str, default='check', choices=['check', 'fix'])
    parser.add_argument('--mode', type=str, default='top', choices=['flip', 'flip-all', 'id', 'ref-alt', 'stats', 'swap', 'top'])
    parser.add_argument('--output-filename', type=str, required=True)
    parser.add_argument('--out-dir', type=str, required=True)

    args = parser.parse_args()

    if args.local:
        backend = hb.LocalBackend()
    else:
        backend = hb.ServiceBackend(billing_project=args.billing_project,
                                    remote_tmpdir=f'{args.out_dir}/tmp/')

    run_checks_fix(
        backend=backend,
        input_vcf=args.input_vcf,
        ref_path=args.ref_fasta,
        step=args.step,
        fix_mode=args.mode,
        output_filename=args.output_filename,
        out_dir=args.out_dir)

    backend.close()
