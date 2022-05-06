__author__ = 'Michael Wilson & Lindo Nkambule'

import hailtop.batch as hb
import hail as hl
import pandas as pd
from gwaspy.utils.get_file_size import bytes_to_gb
from gwaspy.phasing.get_filebase import get_vcf_filebase
from typing import List
from typing import Union


def concat_vcfs(b: hb.batch.Batch,
                vcf_basename: str = None,
                vcfs_to_merge: List = None,
                output_type: str = 'bcf',
                software: str = None,
                chrom: str = None,
                docker_img: str = 'docker.io/lindonkambule/gwaspy:v1',
                cpu: int = 8,
                out_dir: str = None):

    global index_cmd

    out_type = 'b' if output_type == 'bcf' else 'z'
    threads = cpu - 1
    vcfs_sizes_sum = 0
    merge_vcf_i = ''

    out_filename = f'{vcf_basename}.{chrom}.phased.{software}.bcf' if output_type == 'bcf' else \
        f'{vcf_basename}.{chrom}.phased.{software}.vcf.gz'
    out_index_name = f'{vcf_basename}.{chrom}.phased.{software}.bcf.csi' if output_type == 'bcf' else \
        f'{vcf_basename}.{chrom}.phased.{software}.vcf.gz.csi'

    for line in vcfs_to_merge:
        vcfs_sizes_sum += 1 + bytes_to_gb(line)

    mem = 'highmem' if vcfs_sizes_sum > 2 else 'standard'
    disk_size = 10 + vcfs_sizes_sum

    concat = b.new_job(name=f'concat-{vcf_basename}')
    concat.memory(mem)
    concat.storage(f'{disk_size}Gi')
    concat.image(docker_img)
    concat.cpu(cpu)

    for line in vcfs_to_merge:
        input_vcf = b.read_input_group(vcf=line,
                                       ind=f'{line}.csi')
        merge_vcf_i += f'{input_vcf.vcf} '

    cmd = f'''
        bcftools concat \
            --no-version \
            --output-type {out_type} \
            --output {out_filename} \
            --threads {threads} \
            --ligate \
            {merge_vcf_i}
    '''

    concat.command(cmd)
    # index the merged output
    concat.command(f'bcftools index {out_filename}')

    concat.command(f'mv {out_filename} {concat.ofile}')
    concat.command(f'mv {out_index_name} {concat.idx}')
    b.write_output(concat.ofile, f'{out_dir}/GWASpy/{vcf_basename}/Phasing/phased_merged/{out_filename}')
    b.write_output(concat.idx, f'{out_dir}/GWASpy/{vcf_basename}/Phasing/phased_merged/{out_index_name}')


def run_concat(backend: Union[hb.ServiceBackend, hb.LocalBackend] = None,
               input_vcf: str = None,
               output_type: str = 'bcf',
               reference: str = 'GRCh38',
               software: str = None,
               out_dir: str = None):

    print(f'\n3. CONCAT {input_vcf}\n')
    vcf_filebase = get_vcf_filebase(input_vcf)

    # get the regions so we can map each file to its specific region
    regions = pd.read_csv(f'{out_dir}/GWASpy/{vcf_filebase}/Phasing/regions.lines', sep='\t', names=['reg', 'ind'])
    regions_dict = pd.Series(regions.reg.values, index=regions.ind).to_dict()

    concat_b = hb.Batch(backend=backend, name=f'concat-phased-chunks-{vcf_filebase}')

    if software == 'shapeit':
        phased_vcf_chunks = hl.utils.hadoop_ls(f'{out_dir}/GWASpy/{vcf_filebase}/Phasing/phased_scatter/*.shapeit.bcf')
    else:
        phased_vcf_chunks = hl.utils.hadoop_ls(f'{out_dir}/GWASpy/{vcf_filebase}/Phasing/phased_scatter/*.eagle.bcf')

    for i in range(1, 24):
        if reference == 'GRCh38':
            if i == 23:
                chrom = 'chrX'
            else:
                chrom = f'chr{i}'

            out_chrom_name = chrom
        else:
            chrom = str(i)
            out_chrom_name = f'chr{chrom}'

        chrom_phased_files_to_concat = []

        for file in phased_vcf_chunks:
            f = file['path']
            vcf_basename = get_vcf_filebase(f)
            file_index = int(vcf_basename.split('.')[-3])
            file_region = regions_dict[file_index]
            map_chrom = file_region.split(':')[0]
            if map_chrom == chrom:
                chrom_phased_files_to_concat.append(f)

        # naturally sort the list of files to merge
        from gwaspy.utils.natural_sort import natural_keys
        chrom_phased_files_to_concat.sort(key=natural_keys)

        # checkpoint to see if file already exists to avoid redoing things
        chrom_out_filename = f'{vcf_filebase}.{out_chrom_name}.phased.{software}.bcf' if output_type == 'bcf' else \
            f'{vcf_filebase}.{out_chrom_name}.phased.{software}.vcf.gz'
        chrom_out_filname_path = f'{out_dir}/GWASpy/{vcf_filebase}/Phasing/phased_merged/{chrom_out_filename}'

        if hl.hadoop_exists(chrom_out_filname_path):
            continue
        else:
            concat_vcfs(b=concat_b, vcfs_to_merge=chrom_phased_files_to_concat, vcf_basename=vcf_filebase,
                        output_type=output_type, software=software, chrom=out_chrom_name, out_dir=out_dir)

    concat_b.run()
