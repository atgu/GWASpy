__author__ = 'Michael Wilson & Lindo Nkambule'

import hailtop.batch as hb
from typing import List
from gwaspy.utils.get_file_size import bytes_to_gb


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
