__author__ = 'Michael Wilson & Lindo Nkambule'

import hailtop.batch as hb
import hail as hl
import argparse
import pandas as pd
from gwaspy.phasing.get_filebase import get_vcf_filebase
from gwaspy.phasing.scatter_vcf import create_windows_bed, vcf_scatter


def eagle_phasing(b: hb.batch.Batch,
                  vcf_file: str = None,
                  ref_vcf: hb.resource.ResourceFile = None,
                  reference: str = 'GRCh38',
                  cpu: int = 8,
                  memory: str = 'standard',
                  storage: int = 50,
                  img: str = 'docker.io/lindonkambule/gwaspy:v1',
                  threads: int = 16,
                  out_dir: str = None):

    vcf_filename_no_ext = get_vcf_filebase(vcf_file)
    output_file_name = f'{vcf_filename_no_ext}.phased.eagle'
    vcf = b.read_input(vcf_file)

    map_file = '/opt/genetic_map_hg38_withX.txt.gz' if reference == 'GRCh38' else '/opt/genetic_map_hg19_withX.txt.gz'

    phase = b.new_job(name=output_file_name)
    phase.cpu(cpu)
    phase.memory(memory)
    phase.storage(f'{storage}Gi')
    phase.image(img)

    if ref_vcf:
        cmd = f'''
        eagle \
            --geneticMapFile {map_file} \
            --numThreads {threads} \
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
            --outPrefix {output_file_name} \
            --vcfOutFormat z \
            --vcf {vcf}
        '''

    phase.command(cmd)

    phase.command(f'mv {output_file_name}.vcf.gz {phase.ofile}')
    b.write_output(phase.ofile, f'{out_dir}/{output_file_name}.vcf.gz')

    return phase


def shapeit_phasing(b: hb.batch.Batch,
                    vcf_file: str = None,
                    ref_vcf: hb.resource.ResourceFile = None,
                    reference: str = 'GRCh38',
                    cpu: int = 8,
                    memory: str = 'standard',
                    storage: int = 50,
                    img: str = 'docker.io/lindonkambule/gwaspy:v1',
                    threads: int = 16,
                    out_dir: str = None):

    vcf_filename_no_ext = get_vcf_filebase(vcf_file)
    output_file_name = f'{vcf_filename_no_ext}.phased.shapeit.vcf.gz'
    vcf = b.read_input(vcf_file)

    map_file = '/opt/genetic_map_hg38_withX.txt.gz' if reference == 'GRCh38' else '/opt/genetic_map_hg19_withX.txt.gz'

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
            --reference {ref_vcf} \
            --output {output_file_name} \
            --thread {threads}
        '''

    else:
        cmd = f'''
        shapeit4.2 \
            --input {vcf} \
            --map {map_file} \
            --output {output_file_name} \
            --thread {threads}
        '''

    # shapeit requires that the VCF be indexed
    phase.command(f'bcftools index {vcf}')
    phase.command(cmd)

    phase.command(f'mv {output_file_name} {phase.ofile}')
    b.write_output(phase.ofile, f'{out_dir}/{output_file_name}')

    return phase


def haplotype_phasing(input_vcfs: str = None,
                      vcf_ref: str = None,
                      local: bool = False,
                      billing_project: str = None,
                      bucket: str = None,
                      software: str = 'shapeit',
                      reference: str = 'GRCh38',
                      max_win_size_cm: float = 10.0,
                      overlap_size_cm: float = 2.0,
                      scatter_memory: int = 26,
                      cpu: int = 8,
                      memory: str = 'standard',
                      storage: int = 50,
                      threads: int = 16,
                      out_dir: str = None):
    # Error handling
    if software.lower() not in ['eagle', 'shapeit']:
        raise SystemExit(f'Incorrect software {software} selected. Options are [eagle, shapeit]')

    if reference not in ['GRCh37', 'GRCh38']:
        raise SystemExit(f'Incorrect reference genome build {reference} selected. Options are [GRCh37, GRCh38]')

    if memory not in ['lowmem', 'standard', 'highmem']:
        raise SystemExit(f'Incorrect memory type {memory} selected. Options are [lowmem, standard, highmem]')

    if not out_dir:
        raise SystemExit('Output directory not specified. Specify using --out_dir if running from the command line or'
                         'out_dir argument if running inside a Python script')

    if local:
        backend = hb.LocalBackend()
    else:
        backend = hb.ServiceBackend(billing_project=billing_project,
                                    bucket=bucket)

    ###################################### SCATTER VCF FILE ######################################
    create_windows_bed(reference=reference, max_win_size_cm=max_win_size_cm, out_dir=out_dir,
                       overlap_size_cm=overlap_size_cm)

    scatter_job = hb.Batch(backend=backend, name='scatter-vcf')

    vcf_paths = pd.read_csv(input_vcfs, sep='\t', header=None)

    for index, row in vcf_paths.iterrows():
        vcf = row[0]

        vcf_scatter(b=scatter_job, vcf_file=vcf, intervals_bed=f'{out_dir}/GWASpy/Phasing/gwaspy.refscatter.bed',
                    memory=scatter_memory, out_dir=out_dir)

    scatter_job.run()

    ######################################  PHASING ######################################

    phasing = hb.Batch(backend=backend,
                       name=f'haplotype-phasing-{software}')

    if vcf_ref:
        vcf_ref = phasing.read_input(vcf_ref)
        print('RUNNING PHASING WITH A REFERENCE PANEL\n')
    else:
        vcf_ref = None
        print('RUNNING PHASING WITHOUT A REFERENCE PANEL\n')

    vcf_paths = pd.read_csv(input_vcfs, sep='\t', header=None)

    for index, row in vcf_paths.iterrows():
        vcf = row[0]
        vcf_filebase = get_vcf_filebase(vcf)
        scatter_vcfs_paths = hl.utils.hadoop_ls(f'{out_dir}/GWASpy/Phasing/{vcf_filebase}/scatter_vcfs')

        vcfs = []
        for i in scatter_vcfs_paths:
            vcfs.append(i['path'])

        phased_vcf_out_dir = f'{out_dir}/GWASpy/Phasing/{vcf_filebase}/phased_scatter'

        for file in vcfs:
            if software == 'eagle':
                eagle_phasing(b=phasing, vcf_file=file, ref_vcf=vcf_ref, reference=reference, cpu=cpu, memory=memory,
                              storage=storage, threads=threads, out_dir=phased_vcf_out_dir)

            else:
                shapeit_phasing(b=phasing, vcf_file=file, ref_vcf=vcf_ref, reference=reference, cpu=cpu,
                                memory=memory, storage=storage, threads=threads, out_dir=phased_vcf_out_dir)

    phasing.run()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-vcfs', type=str, required=True)
    parser.add_argument('--vcf-ref', type=str, default=None)
    parser.add_argument('--local', action='store_true')
    parser.add_argument('--billing-project', required=True)
    parser.add_argument('--bucket', required=True)
    parser.add_argument('--software', type=str, default='shapeit', choices=['eagle', 'shapeit'])
    parser.add_argument('--reference', type=str, default='GRCh38', choices=['GRCh37', 'GRCh38'])
    parser.add_argument('--max-win-size-cm', type=float, default=10.0)
    parser.add_argument('--overlap-size-cm', type=float, default=2.0)
    parser.add_argument('--cpu', type=int, default=8)
    parser.add_argument('--scatter-mem', type=int, default=26)
    parser.add_argument('--memory', type=str, default='standard', choices=['lowmem', 'standard', 'highmem'])
    parser.add_argument('--storage', type=int, default=50)
    parser.add_argument('--threads', type=int, default=16)
    parser.add_argument('--out-dir', required=True)

    args = parser.parse_args()

    haplotype_phasing(input_vcfs=args.input_vcfs, vcf_ref=args.vcf_ref, local=args.local, bucket=args.bucket,
                      billing_project=args.billing_project, software=args.software, reference=args.reference,
                      max_win_size_cm=args.max_win_size_cm, overlap_size_cm=args.overlap_size_cm,
                      scatter_memory=args.scatter_mem, cpu=args.cpu, memory=args.memory, storage=args.storage,
                      threads=args.threads, out_dir=args.out_dir)


if __name__ == '__main__':
    main()
