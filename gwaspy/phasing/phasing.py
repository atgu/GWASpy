__author__ = 'Michael Wilson & Lindo Nkambule'

import hailtop.batch as hb
import ntpath
import argparse
import pandas as pd
from typing import Union


# phasing without a reference panel (this is fine for big dataset)
# we want to use HGDP+1KG as reference panel, but it hasn't been phased yet.
def eagle_phasing(b: hb.batch.Batch,
                  vcf: hb.resource.ResourceFile,
                  vcf_filename_no_ext: str = None,
                  ref_vcf: hb.resource.ResourceFile = None,
                  reference: str = 'GRCh38',
                  contig: Union[str, int] = None,
                  cpu: int = 8,
                  memory: str = 'standard',
                  storage: int = 50,
                  img: str = 'docker.io/lindonkambule/gwaspy:v1',
                  threads: int = 16,
                  out_dir: str = None):

    output_file_name = vcf_filename_no_ext + '_' + str(contig) + '.phased.eagle'

    map_file = f'/opt/genetic_maps_eagle/hg38/genetic_map_hg38_chr{contig}_withX.txt.gz' if reference == 'GRCh38' else f'/opt/genetic_maps_eagle/hg19/genetic_map_hg19_chr{contig}_withX.txt.gz'

    phase = b.new_job(name=output_file_name)
    phase.cpu(cpu)
    phase.memory(memory)
    phase.storage(f'{storage}Gi')
    phase.image(img)
    # phase.declare_resource_group(ofile={'vcf': '{root}.vcf.gz'})

    if ref_vcf:
        cmd = f'''
        eagle \
            --geneticMapFile {map_file} \
            --numThreads {threads} \
            --chrom {contig} \
            --outPrefix {output_file_name} \
            --vcfOutFormat z \
            --vcfRef {ref_vcf} \
            --vcfTarget {vcf}
        '''

    else:
        cmd = f'''
        eagle \
            --geneticMapFile {map_file} \
            --numThreads {threads} \
            --chrom {contig} \
            --outPrefix {output_file_name} \
            --vcfOutFormat z \
            --vcf {vcf}
        '''

    phase.command(cmd)

    phase.command(f'mv {output_file_name}.vcf.gz {phase.ofile}')
    b.write_output(phase.ofile, f'{out_dir}/GWASpy/Phasing/{vcf_filename_no_ext}/{output_file_name}.vcf.gz')

    return phase


def shapeit_phasing(b: hb.batch.Batch,
                    vcf: hb.resource.ResourceFile,
                    vcf_filename_no_ext: str = None,
                    ref_vcf: hb.resource.ResourceFile = None,
                    reference: str = 'GRCh38',
                    contig: Union[str, int] = None,
                    cpu: int = 8,
                    memory: str = 'standard',
                    storage: int = 50,
                    img: str = 'docker.io/lindonkambule/gwaspy:v1',
                    threads: int = 16,
                    out_dir: str = None):

    output_file_name = vcf_filename_no_ext + '_' + str(contig) + '.phased.shapeit.vcf.gz'

    map_file = f'/opt/genetic_maps_shapeit/hg38/genetic_map_hg38_chr{contig}_withX.txt.gz' if reference == 'GRCh38' else f'/opt/genetic_maps_shapeit/hg19/genetic_map_hg19_chr{contig}_withX.txt.gz'

    phase = b.new_job(name=output_file_name)
    phase.cpu(cpu)
    phase.memory(memory)
    phase.storage(f'{storage}Gi')
    phase.image(img)

    if ref_vcf:
        # shapeit requires that the VCF be indexed
        phase.command(f'bcftools index {ref_vcf}')
        cmd = f'''
        shapeit4.2 \
            --input {vcf} \
            --map {map_file} \
            --region {contig} \
            --reference {ref_vcf} \
            --output {output_file_name} \
            --thread {threads}
        '''

    else:
        cmd = f'''
        shapeit4.2 \
            --input {vcf} \
            --map {map_file} \
            --region {contig} \
            --output {output_file_name} \
            --thread {threads}
        '''

    # shapeit requires that the VCF be indexed
    phase.command(f'bcftools index {vcf}')
    phase.command(cmd)

    phase.command(f'mv {output_file_name} {phase.ofile}')
    b.write_output(phase.ofile, f'{out_dir}/GWASpy/Phasing/{vcf_filename_no_ext}/{output_file_name}')

    return phase


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-vcfs', type=str, required=True)
    parser.add_argument('--vcf-ref', type=str, default=None)
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

    if args.vcf_ref:
        vcf_ref = phasing.read_input(args.vcf_ref)
        print('RUNNING PHASING WITH A REFERENCE PANEL\n')
    else:
        vcf_ref = None
        print('RUNNING PHASING WITHOUT A REFERENCE PANEL\n')

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
                phasing_job = eagle_phasing(b=phasing, vcf=in_vcf, vcf_filename_no_ext=file_no_ext, ref_vcf=vcf_ref,
                                            reference=args.reference, contig=chrom, cpu=args.cpu, memory=args.memory,
                                            storage=args.storage, threads=args.threads, out_dir=args.out_dir)

            else:
                phasing_job = shapeit_phasing(b=phasing, vcf=in_vcf, vcf_filename_no_ext=file_no_ext, ref_vcf=vcf_ref,
                                              reference=args.reference, contig=chrom, cpu=args.cpu, memory=args.memory,
                                              storage=args.storage, threads=args.threads, out_dir=args.out_dir)

    phasing.run()


if __name__ == '__main__':
    main()
