__author__ = 'Michael Wilson & Lindo Nkambule'

import hailtop.batch as hb
import argparse


def haplotype_phasing(input_vcf: str = None,
                      vcf_ref: str = None,
                      local: bool = False,
                      billing_project: str = None,
                      software: str = 'shapeit',
                      reference: str = 'GRCh38',
                      max_win_size_cm: float = 10.0,
                      overlap_size_cm: float = 2.0,
                      scatter_memory: int = 26,
                      cpu: int = 8,
                      threads: int = 7,
                      stages: str = 'scatter,phase,concat',
                      output_type: str = 'bcf',
                      out_dir: str = None):
    # Error handling
    if not out_dir:
        raise SystemExit('Output directory not specified. Specify using --out_dir if running from the command line or'
                         'out_dir argument if running inside a Python script')

    steps_list = stages.split(',')
    steps = [x.lower() for x in steps_list]
    unknown_steps = [i for i in steps if i not in ['scatter', 'phase', 'concat']]

    if len(unknown_steps) > 0:
        raise SystemExit(f'Incorrect process(es) {unknown_steps} selected. Options are [scatter, phase, concat]')

    if output_type.lower() not in ['bcf', 'vcf']:
        raise SystemExit(f'Incorrect output type {output_type} selected. Options are [bcf, vcf]')

    if local:
        backend = hb.LocalBackend()
    else:
        backend = hb.ServiceBackend(billing_project=billing_project,
                                    remote_tmpdir=f'{out_dir}/tmp/')

    # Scatter VCF/BCF file(s)
    if 'scatter' in steps:
        from gwaspy.phasing.scatter_vcf import run_scatter
        run_scatter(backend=backend, input_vcf=input_vcf, reference=reference, max_win_size_cm=max_win_size_cm,
                    overlap_size_cm=overlap_size_cm, scatter_memory=scatter_memory, out_dir=out_dir)

    # Phase scatterd chunks
    if 'phase' in steps:
        from gwaspy.phasing.phase_vcf import run_phase
        run_phase(backend=backend, input_vcf=input_vcf, vcf_ref_path=vcf_ref, software=software, reference=reference,
                  cpu=cpu, threads=threads, out_dir=out_dir)

    # Concatenate phased chunks
    if 'concat' in steps:
        from gwaspy.phasing.concat_vcfs import run_concat
        run_concat(backend=backend, input_vcf=input_vcf, output_type=output_type, reference=reference,
                   software=software, out_dir=out_dir)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-vcf', type=str, required=True)
    parser.add_argument('--vcf-ref', type=str, default=None)
    parser.add_argument('--local', action='store_true')
    parser.add_argument('--billing-project', required=True)
    parser.add_argument('--software', type=str, default='shapeit', choices=['eagle', 'shapeit'])
    parser.add_argument('--reference', type=str, default='GRCh38', choices=['GRCh37', 'GRCh38'])
    parser.add_argument('--max-win-size-cm', type=float, default=10.0)
    parser.add_argument('--overlap-size-cm', type=float, default=2.0)
    parser.add_argument('--cpu', type=int, default=8)
    parser.add_argument('--scatter-mem', type=int, default=26)
    parser.add_argument('--threads', type=int, default=7)
    parser.add_argument('--stages', type=str, default='scatter,phase,concat')
    parser.add_argument('--out-type', type=str, default='bcf', choices=['bcf', 'vcf'])
    parser.add_argument('--out-dir', required=True)

    args = parser.parse_args()

    haplotype_phasing(input_vcf=args.input_vcf, vcf_ref=args.vcf_ref, local=args.local,
                      billing_project=args.billing_project, software=args.software, reference=args.reference,
                      max_win_size_cm=args.max_win_size_cm, overlap_size_cm=args.overlap_size_cm,
                      scatter_memory=args.scatter_mem, cpu=args.cpu, threads=args.threads, stages=args.stages,
                      output_type=args.out_type, out_dir=args.out_dir)


if __name__ == '__main__':
    main()
