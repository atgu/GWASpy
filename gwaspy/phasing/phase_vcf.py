__author__ = 'Lindo Nkambule'

import hailtop.batch as hb
import hail as hl
import pandas as pd
from gwaspy.phasing.get_filebase import get_vcf_filebase
from gwaspy.utils.get_file_size import bytes_to_gb
from typing import Union


def eagle_phasing(b: hb.batch.Batch,
                  vcf_file: str = None,
                  ref_vcf_file: hb.ResourceFile = None,
                  ref_size: Union[float, int] = None,
                  reference: str = 'GRCh38',
                  cpu: int = 8,
                  img: str = 'docker.io/lindonkambule/gwaspy:v1',
                  threads: int = 7,
                  out_dir: str = None):

    vcf_size = bytes_to_gb(vcf_file)
    mem = 'highmem' if vcf_size > 1 else 'standard'
    disk_size = round(5.0 + 3.0 * vcf_size) + round(5.0 + 3.0 * ref_size) if ref_vcf_file else round(5.0 + 3.0 * vcf_size)

    vcf_filename_no_ext = get_vcf_filebase(vcf_file)
    output_file_name = f'{vcf_filename_no_ext}.phased.withref.eagle' if ref_vcf_file else \
        f'{vcf_filename_no_ext}.phased.eagle'
    vcf = b.read_input(vcf_file)

    map_file = '/opt/genetic_map_hg38_withX.txt.gz' if reference == 'GRCh38' else '/opt/genetic_map_hg19_withX.txt.gz'

    phase = b.new_job(name=f'phase-{vcf_file}')
    phase.cpu(cpu)
    phase.memory(mem)
    phase.storage(f'{disk_size}Gi')
    phase.image(img)

    if ref_vcf_file:
        phase.command('echo INDEXING REFERENCE FILE')
        phase.command(f'bcftools index {ref_vcf_file}')
        cmd = f'''
        eagle \
            --geneticMapFile {map_file} \
            --numThreads {threads} \
            --outPrefix {output_file_name} \
            --vcfOutFormat b \
            --vcfRef {ref_vcf_file} \
            --vcfTarget {vcf}
        '''

    else:
        cmd = f'''
        eagle \
            --geneticMapFile {map_file} \
            --numThreads {threads} \
            --outPrefix {output_file_name} \
            --vcfOutFormat b \
            --vcf {vcf}
        '''

    phase.command('echo INDEXING INPUT FILE')
    phase.command(f'bcftools index {vcf}')
    phase.command('echo RUNNING PHASING')
    phase.command(cmd)
    # index the output
    phase.command(f'bcftools index {output_file_name}.bcf')

    phase.command(f'mv {output_file_name}.bcf {phase.ofile}')
    phase.command(f'mv {output_file_name}.bcf.csi {phase.csi}')
    b.write_output(phase.ofile, f'{out_dir}/{output_file_name}.bcf')
    b.write_output(phase.csi, f'{out_dir}/{output_file_name}.bcf.csi')


def shapeit_phasing(b: hb.batch.Batch,
                    vcf_file: str = None,
                    ref_vcf_file: hb.ResourceFile = None,
                    ref_size: Union[float, int] = None,
                    reference: str = 'GRCh38',
                    region: Union[str, int] = None,
                    map_chromosome: str = None,
                    cpu: int = 8,
                    img: str = 'docker.io/lindonkambule/gwaspy:v1',
                    threads: int = 7,
                    out_dir: str = None):

    vcf_size = bytes_to_gb(vcf_file)
    mem = 'highmem' if vcf_size > 1 else 'standard'
    disk_size = round(5.0 + 3.0 * vcf_size) + round(5.0 + 3.0 * ref_size) if ref_vcf_file else round(5.0 + 3.0 * vcf_size)

    vcf_filename_no_ext = get_vcf_filebase(vcf_file)
    output_file_name = f'{vcf_filename_no_ext}.phased.withref.shapeit.bcf'if ref_vcf_file else \
        f'{vcf_filename_no_ext}.phased.shapeit.bcf'
    vcf = b.read_input(vcf_file)

    chrom = map_chromosome if reference == 'GRCh38' else f'chr{map_chromosome}'

    if reference == 'GRCh37':
        if chrom == 'chr23':
            chrom = 'chrX'

    map_file = f'/shapeit4/maps/b38/{chrom}.b38.gmap.gz' if reference == 'GRCh38'\
        else f'/shapeit4/maps/b37/{chrom}.b37.gmap.gz'

    phase = b.new_job(name=f'phase-{vcf_file}')
    phase.cpu(cpu)
    phase.memory(mem)
    phase.storage(f'{disk_size}Gi')
    phase.image(img)

    if ref_vcf_file:
        phase.command('echo INDEXING REFERENCE FILE')
        phase.command(f'bcftools index {ref_vcf_file}')
        cmd = f'''
        shapeit4.2 \
            --input {vcf} \
            --map {map_file} \
            --region {region} \
            --reference {ref_vcf_file} \
            --output {output_file_name} \
            --thread {threads}
        '''

    else:
        cmd = f'''
        shapeit4.2 \
            --input {vcf} \
            --map {map_file} \
            --region {region} \
            --output {output_file_name} \
            --thread {threads}
        '''

    phase.command('echo INDEXING INPUT FILE')
    phase.command(f'bcftools index {vcf}')
    phase.command('echo RUNNING PHASING')
    phase.command(cmd)
    # index the output
    phase.command(f'bcftools index {output_file_name}')

    phase.command(f'mv {output_file_name} {phase.ofile}')
    phase.command(f'mv {output_file_name}.csi {phase.csi}')
    b.write_output(phase.ofile, f'{out_dir}/{output_file_name}')
    b.write_output(phase.csi, f'{out_dir}/{output_file_name}.csi')


def run_phase(backend: Union[hb.ServiceBackend, hb.LocalBackend] = None,
              input_vcf: str = None,
              vcf_ref_path: str = None,
              software: str = 'shapeit',
              reference: str = 'GRCh38',
              cpu: int = 8,
              threads: int = 7,
              out_dir: str = None):

    # error handling
    global ref_path, ref_chr, chrom
    if software.lower() not in ['eagle', 'shapeit']:
        raise SystemExit(f'Incorrect software {software} selected. Options are [eagle, shapeit]')

    if reference not in ['GRCh37', 'GRCh38']:
        raise SystemExit(f'Incorrect reference genome build {reference} selected. Options are [GRCh37, GRCh38]')

    vcf_filebase = get_vcf_filebase(input_vcf)
    phasing_b = hb.Batch(backend=backend,
                         name=f'phasing-{vcf_filebase}-{software}')

    if vcf_ref_path:
        if vcf_ref_path == 'hgdp_1kg':
            print(f'\n2. PHASING {input_vcf} WITH HGDP + 1000 GENOMES REFERENCE PANEL\n')
            ref_path = 'gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes/hgdp.tgp.gwaspy.merged.chrCNUMBER.merged.bcf'
        else:
            print(f'\n2. PHASING {input_vcf} WITH USER-DEFINED REFERENCE PANEL\n')
            ref_path = vcf_ref_path
    else:
        print(f'\n2. PHASING {input_vcf} WITHOUT A REFERENCE PANEL\n')

    # get the regions so we can map each file to its specific region
    regions = pd.read_csv(f'{out_dir}/GWASpy/{vcf_filebase}/Phasing/regions.lines', sep='\t', names=['reg', 'ind'])
    regions_dict = pd.Series(regions.reg.values, index=regions.ind).to_dict()

    scatter_vcfs_paths = hl.utils.hadoop_ls(f'{out_dir}/GWASpy/{vcf_filebase}/Phasing/scatter_vcfs/*.bcf')

    vcfs = []
    for i in scatter_vcfs_paths:
        vcfs.append(i['path'])

    phased_vcf_out_dir = f'{out_dir}/GWASpy/{vcf_filebase}/Phasing/phased_scatter'

    for i in range(1, 24):
        if reference == 'GRCh38':
            if i == 23:
                chrom = 'chrX'
                ref_chr = 'X'
            else:
                chrom = f'chr{i}'
                ref_chr = i

        else:
            chrom = str(i)
            ref_chr = i

        ref_chrom_path = ref_path.replace('CNUMBER', str(ref_chr)) if vcf_ref_path else None
        ref_vcf_chrom_file = phasing_b.read_input(ref_chrom_path) if vcf_ref_path else None
        ref_vcf_chrom_file_size = bytes_to_gb(ref_chrom_path) if vcf_ref_path else None

        for file in vcfs:
            # get specific region for file using regions.line file
            vcf_basename = get_vcf_filebase(file)
            file_index = int(vcf_basename.split('.')[-1])
            file_region = regions_dict[file_index]
            map_chrom = file_region.split(':')[0]

            if map_chrom == chrom:
                if software == 'eagle':
                    # check if file exists to avoid re-doing things
                    out_phased_filename = f'{vcf_basename}.phased.withref.eagle.bcf' if ref_vcf_chrom_file else \
                        f'{vcf_basename}.phased.eagle.bcf'

                    if hl.hadoop_exists(f'{phased_vcf_out_dir}/{out_phased_filename}'):
                        continue
                    else:
                        eagle_phasing(b=phasing_b, vcf_file=file, ref_vcf_file=ref_vcf_chrom_file,
                                      reference=reference, ref_size=ref_vcf_chrom_file_size, cpu=cpu, threads=threads,
                                      out_dir=phased_vcf_out_dir)

                else:
                    out_phased_filename = f'{vcf_basename}.phased.withref.shapeit.bcf' if ref_vcf_chrom_file else \
                        f'{vcf_basename}.phased.shapeit.bcf'

                    # check if file exists to avoid re-doing things
                    if hl.hadoop_exists(f'{phased_vcf_out_dir}/{out_phased_filename}'):
                        continue
                    else:
                        shapeit_phasing(b=phasing_b, vcf_file=file, ref_vcf_file=ref_vcf_chrom_file,
                                        ref_size=ref_vcf_chrom_file_size, reference=reference, region=file_region,
                                        map_chromosome=map_chrom, cpu=cpu, threads=threads, out_dir=phased_vcf_out_dir)

    phasing_b.run()


