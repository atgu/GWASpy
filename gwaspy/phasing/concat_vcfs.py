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
                docker_img: str = 'docker.io/lindonkambule/gwaspy:v1',
                memory: str = 'standard',
                cpu: int = 4,
                out_dir: str = None):

    out_type = 'u' if output_type == 'bcf' else 'z'
    threads = cpu - 1
    vcfs_sizes_sum = 0
    merge_vcf_i = ''

    for line in vcfs_to_merge:
        vcfs_sizes_sum += bytes_to_gb(line)
        input_vcf = b.read_input(line)
        merge_vcf_i += f'{input_vcf} \t'

    out_filename = f'{vcf_basename}.merged.bcf' if output_type == 'bcf' else f'{vcf_basename}.merged.vcf.gz'

    concat = b.new_job(name=f'concat-{vcf_basename}')
    concat.memory(f'{memory}Gi')
    concat.storage(f'{vcfs_sizes_sum}Gi')
    concat.image(docker_img)
    concat.cpu(cpu)

    cmd = f'''
        bcftools concat \
            --no-version \
            --output-type {out_type} \
            --output {out_filename} \
            --threads {threads} \
            "{merge_vcf_i}"
    '''

    concat.command(cmd)

    concat.command(f'mv {out_filename} {concat.ofile}')
    b.write_output(concat.ofile, f'{out_dir}/GWASpy/Phasing/{vcf_basename}/phased_merged')


def run_concat(backend: Union[hb.ServiceBackend, hb.LocalBackend] = None,
               input_vcfs: str = None,
               output_type: str = 'bcf',
               out_dir: str = None):

    concat_b = hb.Batch(backend=backend, name=f'concat-phased-chunks')

    vcf_paths = pd.read_csv(input_vcfs, sep='\t', header=None)

    for index, row in vcf_paths.iterrows():
        vcf = row[0]
        vcf_filebase = get_vcf_filebase(vcf)

        phased_vcfs_chunks = hl.utils.hadoop_ls(f'{out_dir}/GWASpy/Phasing/{vcf_filebase}/phased_scatter')

        vcfs = []
        for i in phased_vcfs_chunks:
            vcfs.append(i['path'])

        concat_vcfs(b=concat_b, vcfs_to_merge=vcfs, vcf_basename=vcf_filebase, output_type=output_type, out_dir=out_dir)

    concat_b.run()
