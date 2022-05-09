__author__ = 'Lindo Nkambule'

import hailtop.batch as hb
import argparse


# ONCE YOU ADD FUNCTIONALITY FOR USING DIFFERENT REF PANELS, CHANGE n_panel* parameters and include cmd args
def genotype_imputation(input_vcf: str = None,
                        females_file: str = None,
                        n_samples: int = None,
                        n_panel_samples: int = 4099,
                        buffer_region: int = 250,
                        phasing_software: str = None,
                        local: bool = False,
                        billing_project: str = None,
                        memory: str = 'highmem',
                        cpu: int = 16,
                        stages: str = 'impute,concat',
                        output_type: str = 'bcf',
                        out_dir: str = None):
    # Error handling
    if not out_dir:
        raise SystemExit('Output directory not specified. Specify using --out_dir if running from the command line or'
                         'out_dir argument if running inside a Python script')

    steps_list = stages.split(',')
    steps = [x.lower() for x in steps_list]
    unknown_steps = [i for i in steps if i not in ['impute', 'concat']]

    if len(unknown_steps) > 0:
        raise SystemExit(f'Incorrect process(es) {unknown_steps} selected. Options are [impute, concat]')

    if output_type.lower() not in ['bcf', 'vcf']:
        raise SystemExit(f'Incorrect output type {output_type} selected. Options are [bcf, vcf]')

    if memory.lower() not in ['lowmem', 'standard', 'highmem']:
        raise SystemExit(f'Incorrect memory type {memory} selected. Options are [lowmem, standard, highmem]')

    if not n_samples:
        raise SystemExit('Number of samples in input data not detected. Specify how many samples (integer), using'
                         '--n-samples if running from the command line or'
                         'n_samples argument if running inside a Python script, are in the input data')

    if local:
        backend = hb.LocalBackend()
    else:
        backend = hb.ServiceBackend(billing_project=billing_project,
                                    remote_tmpdir=f'{out_dir}/tmp/')

    # impute genotypes
    if 'impute' in steps:
        from gwaspy.imputation.sex_aut_imp import run_impute
        run_impute(backend=backend, input_vcf=input_vcf, females_file=females_file, n_samples=n_samples,
                   n_panel_samples=n_panel_samples, phasing_software=phasing_software, memory=memory,
                   buffer_region=buffer_region, out_dir=out_dir)

    # Concatenate imputed chunks
    if 'concat' in steps:
        from gwaspy.imputation.concat_vcfs import run_concat
        run_concat(backend=backend, input_vcf=input_vcf, output_type=output_type, cpu=cpu, memory=memory,
                   out_dir=out_dir)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-vcf', type=str, required=True)
    parser.add_argument('--samples-file', type=str, required=True)
    parser.add_argument('--local', action='store_true')
    parser.add_argument('--billing-project', required=True)
    parser.add_argument('--phasing-software', type=str, default='shapeit', choices=['eagle', 'shapeit'])
    parser.add_argument('--memory', type=str, default='highmem', choices=['lowmem', 'standard', 'highmem'])
    parser.add_argument('--cpu-concat', type=int, default=16)
    parser.add_argument('--n-samples', type=int, required=True)
    parser.add_argument('--buffer-region', type=int, default=250)
    parser.add_argument('--stages', type=str, default='impute,concat')
    parser.add_argument('--out-type', type=str, default='bcf', choices=['bcf', 'vcf'])
    parser.add_argument('--out-dir', required=True)

    args = parser.parse_args()

    genotype_imputation(input_vcf=args.input_vcf, females_file=args.samples_file, n_samples=args.n_samples,
                        buffer_region=args.buffer_region, phasing_software=args.phasing_software, local=args.local,
                        billing_project=args.billing_project, memory=args.memory, cpu=args.cpu_concat,
                        stages=args.stages, output_type=args.out_type, out_dir=args.out_dir)


if __name__ == '__main__':
    main()
