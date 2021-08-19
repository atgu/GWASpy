__author__ = 'Lindo Nkambule'

import hailtop.batch as hb
import argparse
from typing import Union


def imputation(b: hb.batch.Batch,
               vcf: hb.resource.ResourceFile,
               vcf_filename_no_ext: str = None,
               ref_imp5: hb.resource.ResourceFile = None,
               reference: str = 'GRCh38',
               contig: Union[str, int] = None,
               cpu: int = 8,
               memory: str = 'standard',
               storage: int = 50,
               img: str = 'docker.io/lindonkambule/gwaspy:v1',
               threads: int = 16,
               out_dir: str = None):

    output_file_name = vcf_filename_no_ext + '_' + str(contig) + '.imputed.vcf.gz'

    if reference == 'GRCh38':
        if contig == 'chr23':
            in_contig = 'chrX'
        else:
            in_contig = contig

    else:
        if contig == 23:
            in_contig = 'X'
        else:
            in_contig = contig

    map_file = f'/shapeit4/maps/b38/{in_contig}.b38.gmap.gz' if reference == 'GRCh38' else f'/shapeit4/maps/b37/chr{in_contig}.b37.gmap.gz'

    impute = b.new_job(name=output_file_name)
    impute.cpu(cpu)
    impute.memory(memory)
    impute.storage(f'{storage}Gi')
    impute.image(img)

    cmd = f'''
        impute5_1.1.5_static \
            --h {ref_imp5} \
            --m {map_file} \
            --g {vcf} \
            --r {in_contig} \
            --o {output_file_name} \
            --threads {threads}
    '''

    impute.command(cmd)

    impute.command(f'mv {output_file_name} {impute.ofile}')
    b.write_output(impute.ofile, f'{out_dir}/GWASpy/Imputation/{vcf_filename_no_ext}/{output_file_name}')

    return impute


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

    imputation()


if __name__ == '__main__':
    main()
