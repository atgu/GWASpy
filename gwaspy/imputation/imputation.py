__author__ = 'Lindo Nkambule'

import hailtop.batch as hb
import hail as hl
import argparse
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

    in_vcf = b.read_input(vcf)
    vcf_size = bytes_to_gb(vcf)

    output_file_name = vcf_filename_no_ext + '.imputed.vcf.gz'
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
            --g {in_vcf} \
            --r {region} \
            --out-gp-field \
            --o {output_file_name} \
            --threads {threads}
    '''

    impute.command(cmd)

    impute.command(f'mv {output_file_name} {impute.ofile}')
    b.write_output(impute.ofile, f'{out_dir}/GWASpy/Imputation/{vcf_filename_no_ext}/{output_file_name}')

    return impute


def run_impute(backend: Union[hb.ServiceBackend, hb.LocalBackend] = None,
               input_vcfs: str = None,
               memory: str = 'highmem',
               cpu: int = 8,
               threads: int = 7,
               out_dir: str = None):

    print('RUNNING IMPUTATION')
    impute_b = hb.Batch(backend=backend, name=f'impute-phased-chunks')

    vcf_paths = pd.read_csv(input_vcfs, sep='\t', header=None)

    # get the regions so we can map each file to its specific region
    regions = pd.read_csv(f'{out_dir}/GWASpy/Phasing/regions.lines', sep='\t', names=['reg', 'ind'])
    regions_dict = pd.Series(regions.reg.values, index=regions.ind).to_dict()

    for index, row in vcf_paths.iterrows():
        vcf = row[0]
        vcf_filebase = get_vcf_filebase(vcf)

        phased_vcfs_chunks = hl.utils.hadoop_ls(f'{out_dir}/GWASpy/Phasing/{vcf_filebase}/phased_scatter')

        for i in range(1, 24):
            if i == 23:
                chrom = 'chrX'
            else:
                chrom = f'chr{i}'

            ref_bcf = f'gs://hgdp-1kg/hgdp_tgp_phasing/vcf/hgdp.tgp.gwaspy.merged.{chrom}.merged.bcf'
            ref_ind = f'gs://hgdp-1kg/hgdp_tgp_phasing/vcf/hgdp.tgp.gwaspy.merged.{chrom}.merged.bcf.csi'
            ref_size = bytes_to_gb(ref_bcf)
            ref = impute_b.read_input_group(**{'bcf': ref_bcf,
                                               'bcf.csi': ref_ind})

            for file in phased_vcfs_chunks:
                f = file['path']
                vcf_basename = get_vcf_filebase(f)
                file_index = int(vcf_basename.split('.')[-3])
                file_region = regions_dict[file_index]
                map_chrom = file_region.split(':')[0]

                if map_chrom == chrom:
                    imputation(b=impute_b, vcf=f, vcf_filename_no_ext=vcf_filebase, ref=ref, ref_size=ref_size,
                               region=file_region[3:], chromosome=chrom, cpu=cpu, memory=memory,
                               threads=threads, out_dir=out_dir)

    impute_b.run()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-vcfs', type=str, required=True)
    parser.add_argument('--local', action='store_true')
    parser.add_argument('--memory', type=str, default='highmem', choices=['lowmem', 'standard', 'highmem'])
    parser.add_argument('--cpu', type=int, default=8)
    parser.add_argument('--threads', type=int, default=7)
    parser.add_argument('--out-dir', required=True)

    args = parser.parse_args()

    if args.local:
        backend = hb.LocalBackend()
    else:
        backend = hb.ServiceBackend()

    run_impute(backend=backend, input_vcfs=args.input_vcfs, memory=args.memory, cpu=args.cpu, threads=args.threads,
               out_dir=args.out_dir)


if __name__ == '__main__':
    main()
