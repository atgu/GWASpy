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
                cpu: int = 16,
                memory: str = 'standard',
                docker_img: str = 'docker.io/lindonkambule/gwaspy:v1',
                out_dir: str = None):

    global index_cmd

    out_type = 'b' if output_type == 'bcf' else 'z'
    vcfs_sizes_sum = 0
    merge_vcf_i = ''

    out_filename = f'{vcf_basename}.{chrom}.merged.bcf' if output_type == 'bcf' else \
        f'{vcf_basename}.{chrom}.merged.vcf.gz'
    out_index_name = f'{vcf_basename}.{chrom}.merged.bcf.csi' if output_type == 'bcf' else \
        f'{vcf_basename}.{chrom}.merged.vcf.gz.csi'

    for line in vcfs_to_merge:
        vcfs_sizes_sum += 2 + bytes_to_gb(line)

    disk_size = int(round(10 + (2 * vcfs_sizes_sum)))
    threads = cpu - 1

    concat = b.new_job(name=f'concat-{vcf_basename}')
    concat.memory(memory)
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
            {merge_vcf_i}
    '''

    concat.command(cmd)
    # index the merged output
    concat.command(f'bcftools index --force {out_filename}')

    concat.command(f'mv {out_filename} {concat.ofile}')
    concat.command(f'mv {out_index_name} {concat.idx}')
    b.write_output(concat.ofile, f'{out_dir}/GWASpy/{vcf_basename}/Imputation/imputed_merged/{out_filename}')
    b.write_output(concat.idx, f'{out_dir}/GWASpy/{vcf_basename}/Imputation/imputed_merged/{out_index_name}')


def run_concat(backend: Union[hb.ServiceBackend, hb.LocalBackend] = None,
               input_vcf: str = None,
               output_type: str = 'vcf',
               cpu: int = 16,
               memory: str = 'standard',
               out_dir: str = None):

    print(f'\n2. CONCAT {input_vcf}\n')
    vcf_filebase = get_vcf_filebase(input_vcf)
    concat_b = hb.Batch(backend=backend, name=f'concat-imputed-chunks-{vcf_filebase}')

    # get the regions so we can map each file to its specific region
    regions = pd.read_csv(f'{out_dir}/GWASpy/{vcf_filebase}/Imputation/imputation.regions', sep='\t', names=['reg', 'ind'])
    regions_dict = pd.Series(regions.reg.values, index=regions.ind).to_dict()

    imputed_vcfs_chunks = hl.utils.hadoop_ls(f'{out_dir}/GWASpy/{vcf_filebase}/Imputation/imputed_chunks/*.bcf')

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

        # naturally sort the list of files to merge
        from gwaspy.utils.natural_sort import natural_keys
        chrom_phased_files_to_concat.sort(key=natural_keys)

        concat_vcfs(b=concat_b, vcfs_to_merge=chrom_phased_files_to_concat, vcf_basename=vcf_filebase,
                    output_type=output_type, chrom=chrom, cpu=cpu, memory=memory, out_dir=out_dir)

    concat_b.run()
