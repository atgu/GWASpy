__author__ = 'Michael Wilson & Lindo Nkambule'

import hailtop.batch as hb
import ntpath
import argparse
import pandas as pd
from typing import Union


def eagle_phasing(b: hb.batch.Batch,
                  vcf: hb.resource.ResourceFile,
                  vcf_filename_no_ext: str = None,
                  reference: str = 'GRCh38',
                  contig: Union[str, int] = None,
                  cpu: int = 8,
                  memory: str = 'standard',
                  storage: int = 50,
                  img: str = 'gcr.io/broad-mpg-gnomad/lai_phasing:latest',
                  threads: int = 16,
                  out_dir: str = None):

    output_file_name = vcf_filename_no_ext + '_' + str(contig) + '.vcf.gz'

    map_file = 'genetic_map_hg38_withX.txt.gz' if reference == 'GRCh38' else 'genetic_map_hg19_withX.txt.gz'

    phase = b.new_job(name=output_file_name)
    phase.cpu(cpu)
    phase.memory(memory)
    phase.storage(f'{storage}Gi')
    phase.image(img)
    phase.declare_resource_group(ofile={'vcf': '{root}.vcf.gz'})
    cmd = f'''
    Eagle_v2.4.1/eagle \
        --geneticMapFile Eagle_v2.4.1/tables/{map_file} \
        --numThreads {threads} \
        --chrom {contig} \
        --outPrefix {phase.ofile} \
        --vcfOutFormat z \
        --vcf {vcf}
    '''

    phase.command(cmd)

    b.write_output(phase.ofile, f'{out_dir}/GWASpy/Phasing/{vcf_filename_no_ext}/{output_file_name}')

    return phase


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-vcfs', required=True)
    parser.add_argument('--local', action='store_true')
    parser.add_argument('--software', type=str, default='eagle', choices=['eagle', 'shapeit'])
    parser.add_argument('--reference', type=str, default='GRCh38', choices=['GRCh37', 'GRCh38'])
    parser.add_argument('--cpu', type=int, default=8)
    parser.add_argument('--memory', type=str, default='standard', choices=['lowmem', 'standard', 'highmem'])
    parser.add_argument('--storage', type=int, default=50)
    parser.add_argument('--threads', type=int, default=16)
    parser.add_argument('--out-dir', required=True)

    args = parser.parse_args()

    if args.local:
        backend = hb.LocalBackend()
    else:
        backend = hb.ServiceBackend()

    phasing = hb.Batch(backend=backend,
                       name='haplotype-phasing')

    vcf_paths = pd.read_csv(args.input_vcfs, sep='\t', header=None)

    for index, row in vcf_paths.iterrows():
        vcf = row[0]
        in_vcf = phasing.read_input(vcf)
        vcf_name = ntpath.basename(vcf)
        if vcf_name.endswith('.gz'):
            file_no_ext = vcf_name[:-7]
        elif vcf_name.endswith('.bgz'):
            file_no_ext = vcf_name[:-8]
        else:
            file_no_ext = vcf_name[:-4]

        for i in range(1, 24):
            chrom = f'chr{i}' if args.reference == 'GRCh38' else i

            if args.software == 'eagle':
                phasing_job = eagle_phasing(b=phasing, vcf=in_vcf, vcf_filename_no_ext=file_no_ext,
                                            reference=args.reference, contig=chrom, cpu=args.cpu, memory=args.memory,
                                            storage=args.storage, threads=args.threads, out_dir=args.out_dir)

            else:
                print("Support for SHAPEIT coming")

    phasing.run()


if __name__ == '__main__':
    main()
