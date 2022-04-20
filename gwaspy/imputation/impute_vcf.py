__author__ = 'Lindo Nkambule'

import hailtop.batch as hb
import hail as hl
import pandas as pd
from typing import Union
from gwaspy.phasing.get_filebase import get_vcf_filebase
from gwaspy.utils.get_file_size import bytes_to_gb


def imputation(b: hb.batch.Batch,
               vcf: str = None,
               vcf_filename_no_ext: str = None,
               ref: hb.ResourceGroup = None,
               ref_size: Union[int, float] = None,
               region: str = None,
               chromosome: str = None,
               cpu: int = 8,
               memory: str = 'highmem',
               img: str = 'docker.io/lindonkambule/gwaspy:v1',
               threads: int = 7,
               out_dir: str = None):

    # in_vcf = b.read_input(vcf)
    in_vcf = b.read_input_group(**{'bcf': vcf,
                                'bcf.csi': f'{vcf}.csi'})
    vcf_size = bytes_to_gb(vcf)

    output_file_name = vcf_filename_no_ext + '.imputed.bcf'
    file_dir = vcf_filename_no_ext.split('.')[0]

    disk_size = ref_size + (vcf_size * 4)

    map_file = f'/shapeit4/maps/b38/{chromosome}.b38.gmap.gz'

    impute = b.new_job(name=output_file_name)
    impute.cpu(cpu)
    impute.memory(memory)
    impute.storage(f'{disk_size}Gi')
    impute.image(img)

    cmd = f'''
        impute5_1.1.5_static \
            --h {ref.bcf} \
            --m {map_file} \
            --g {in_vcf.bcf} \
            --r {region} \
            --out-gp-field \
            --o {output_file_name} \
            --threads {threads}
    '''

    impute.command(cmd)
    # index file to use when merging
    impute.command(f'bcftools index {output_file_name}')

    impute.command(f'mv {output_file_name} {impute.ofile}')
    impute.command(f'mv {output_file_name}.csi {impute.ind}')
    b.write_output(impute.ofile, f'{out_dir}/GWASpy/Imputation/{file_dir}/imputed_chunks/{output_file_name}')
    b.write_output(impute.ind, f'{out_dir}/GWASpy/Imputation/{file_dir}/imputed_chunks/{output_file_name}.csi')


def run_impute(backend: Union[hb.ServiceBackend, hb.LocalBackend] = None,
               input_vcfs: str = None,
               phasing_software: str = None,
               memory: str = 'highmem',
               cpu: int = 8,
               threads: int = 7,
               out_dir: str = None):

    print(f'RUNNING IMPUTATION ON FILES PHASED WITH {phasing_software.upper()}')
    impute_b = hb.Batch(backend=backend, name=f'impute-phased-chunks')

    vcf_paths = pd.read_csv(input_vcfs, sep='\t', header=None)

    # get the regions so we can map each file to its specific region
    regions = pd.read_csv(f'{out_dir}/GWASpy/Phasing/regions.lines', sep='\t', names=['reg', 'ind'])
    regions_dict = pd.Series(regions.reg.values, index=regions.ind).to_dict()

    for index, row in vcf_paths.iterrows():
        vcf = row[0]
        vcf_filebase = get_vcf_filebase(vcf)

        if phasing_software == 'shapeit':
            phased_vcfs_chunks = hl.utils.hadoop_ls(f'{out_dir}/GWASpy/Phasing/{vcf_filebase}/phased_scatter/*.shapeit.bcf')
        else:
            phased_vcfs_chunks = hl.utils.hadoop_ls(f'{out_dir}/GWASpy/Phasing/{vcf_filebase}/phased_scatter/*.eagle.bcf')

        for i in range(1, 24):
            if i == 23:
                chrom = 'chrX'
            else:
                chrom = f'chr{i}'

            ref_bcf = f'gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes/hgdp.tgp.gwaspy.merged.{chrom}.merged.bcf'
            ref_ind = f'gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes/hgdp.tgp.gwaspy.merged.{chrom}.merged.bcf.csi'
            ref_size = bytes_to_gb(ref_bcf)
            ref = impute_b.read_input_group(**{'bcf': ref_bcf,
                                               'bcf.csi': ref_ind})

            for file in phased_vcfs_chunks:
                f = file['path']
                vcf_basename = get_vcf_filebase(f)
                file_index = int(vcf_basename.split('.')[-3])
                file_region = regions_dict[file_index]
                map_chrom = file_region.split(':')[0]

                imp_out_filename = f'{vcf_basename}.imputed.bcf'
                file_dir = vcf_basename.split('.')[0]
                output_filepath_name = f'{out_dir}/GWASpy/Imputation/{file_dir}/imputed_chunks/{imp_out_filename}'

                if map_chrom == chrom:
                    # check if imputed file already exists
                    if hl.hadoop_exists(output_filepath_name):
                        continue
                    else:
                        imputation(b=impute_b, vcf=f, vcf_filename_no_ext=vcf_basename, ref=ref, ref_size=ref_size,
                                   region=file_region, chromosome=chrom, cpu=cpu, memory=memory,
                                   threads=threads, out_dir=out_dir)

    impute_b.run()
