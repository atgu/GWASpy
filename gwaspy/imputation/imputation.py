__author__ = 'Lindo Nkambule'

import hailtop.batch as hb
import argparse


def genotype_imputation(input_vcfs: str = None,
                        females_file: str = None,
                        local: bool = False,
                        billing_project: str = None,
                        bucket: str = None,
                        memory: str = 'highmem',
                        cpu: int = 8,
                        threads: int = 7,
                        run: str = 'impute',
                        output_type: str = 'bcf',
                        out_dir: str = None):
    # Error handling
    if not out_dir:
        raise SystemExit('Output directory not specified. Specify using --out_dir if running from the command line or'
                         'out_dir argument if running inside a Python script')

    if run.lower() not in ['impute', 'concat']:
        raise SystemExit(f'Incorrect process {run} selected. Options are [impute, concat]')

    if output_type.lower() not in ['bcf', 'vcf']:
        raise SystemExit(f'Incorrect output type {run} selected. Options are [bcf, vcf]')

    if memory.lower() not in ['lowmem', 'standard', 'highmem']:
        raise SystemExit(f'Incorrect memory type {run} selected. Options are [lowmem, standard, highmem]')

    if local:
        backend = hb.LocalBackend()
    else:
        backend = hb.ServiceBackend(billing_project=billing_project,
                                    bucket=bucket)

    # impute genotypes
    if run.lower() == 'impute':
        from gwaspy.imputation.sex_aut_imp import run_impute
        run_impute(backend=backend, input_vcfs=input_vcfs, females_file=females_file, memory=memory, cpu=cpu,
                   threads=threads, out_dir=out_dir)

    # Concatenate imputed chunks
    if run.lower() == 'concat':
        from gwaspy.imputation.concat_vcfs import run_concat
        run_concat(backend=backend, input_vcfs=input_vcfs, output_type=output_type, out_dir=out_dir)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-vcfs', type=str, required=True)
    parser.add_argument('--samples-file', type=str, required=True)
    parser.add_argument('--local', action='store_true')
    parser.add_argument('--billing-project', required=True)
    parser.add_argument('--bucket', required=True)
    parser.add_argument('--memory', type=str, default='highmem', choices=['lowmem', 'standard', 'highmem'])
    parser.add_argument('--cpu', type=int, default=8)
    parser.add_argument('--threads', type=int, default=7)
    parser.add_argument('--run', type=str, default='impute', choices=['impute', 'concat'])
    parser.add_argument('--out-type', type=str, default='bcf', choices=['bcf', 'vcf'])
    parser.add_argument('--out-dir', required=True)

    args = parser.parse_args()

    genotype_imputation(input_vcfs=args.input_vcfs, females_file=args.samples_file, local=args.local,
                        billing_project=args.billing_project, bucket=args.bucket, memory=args.memory, cpu=args.cpu,
                        threads=args.threads, run=args.run, output_type=args.out_type, out_dir=args.out_dir)


if __name__ == '__main__':
    main()
