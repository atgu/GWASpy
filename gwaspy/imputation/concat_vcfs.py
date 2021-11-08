__author__ = 'Lindo Nkambule'

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
                output_type: str = 'vcf',
                chrom: str = None,
                docker_img: str = 'docker.io/lindonkambule/gwaspy:v1',
                cpu: int = 8,
                out_dir: str = None):

    global index_cmd

    out_type = 'b' if output_type == 'bcf' else 'z'
    threads = cpu - 1
    vcfs_sizes_sum = 0
    merge_vcf_i = ''

    out_filename = f'{vcf_basename}.{chrom}.merged.bcf' if output_type == 'bcf' else \
        f'{vcf_basename}.{chrom}.merged.vcf.gz'
    out_index_name = f'{vcf_basename}.{chrom}.merged.bcf.csi' if output_type == 'bcf' else \
        f'{vcf_basename}.{chrom}.merged.vcf.gz.csi'

    for line in vcfs_to_merge:
        vcfs_sizes_sum += 2 + bytes_to_gb(line)

    mem = 'highmem' if vcfs_sizes_sum > 2 else 'standard'

    concat = b.new_job(name=f'concat-{vcf_basename}')
    concat.memory(mem)
    concat.storage(f'{vcfs_sizes_sum}Gi')
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
    b.write_output(concat.ofile, f'{out_dir}/GWASpy/Imputation/{vcf_basename}/imputed_merged/{out_filename}')
    b.write_output(concat.idx, f'{out_dir}/GWASpy/Imputation/{vcf_basename}/imputed_merged/{out_index_name}')


def run_concat(backend: Union[hb.ServiceBackend, hb.LocalBackend] = None,
               input_vcfs: str = None,
               output_type: str = 'vcf',
               out_dir: str = None):

    concat_b = hb.Batch(backend=backend, name=f'concat-phased-chunks')

    vcf_paths = pd.read_csv(input_vcfs, sep='\t', header=None)

    # get the regions so we can map each file to its specific region
    regions = pd.read_csv(f'{out_dir}/GWASpy/Phasing/regions.lines', sep='\t', names=['reg', 'ind'])
    regions_dict = pd.Series(regions.reg.values, index=regions.ind).to_dict()

    for index, row in vcf_paths.iterrows():
        vcf = row[0]
        vcf_filebase = get_vcf_filebase(vcf)

        imputed_vcfs_chunks = hl.utils.hadoop_ls(f'{out_dir}/GWASpy/Imputation/{vcf_filebase}/imputed_chunks/*.vcf.gz')

        for i in range(1, 24):
            if i == 23:
                chrom = 'chrX'
            else:
                chrom = f'chr{i}'

            chrom_phased_files_to_concat = []

            for file in imputed_vcfs_chunks:
                f = file['path']
                vcf_basename = get_vcf_filebase(f)
                file_index = int(vcf_basename.split('.')[-4])
                file_region = regions_dict[file_index]
                map_chrom = file_region.split(':')[0]
                if map_chrom == chrom:
                    chrom_phased_files_to_concat.append(f)

            concat_vcfs(b=concat_b, vcfs_to_merge=chrom_phased_files_to_concat, vcf_basename=vcf_filebase,
                        output_type=output_type, chrom=chrom, out_dir=out_dir)

    concat_b.run()
