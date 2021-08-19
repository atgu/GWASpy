__author__ = 'Lindo Nkambule'

import hailtop.batch as hb
import ntpath
import argparse
import pandas as pd
from typing import Union


def imputation(b: hb.batch.Batch,
               vcf: hb.resource.ResourceFile,
               vcf_filename_no_ext: str = None,
               ref_imp5: hb.resource.ResourceGroup = None,
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
            --h {ref_imp5.imp5} \
            --m {map_file} \
            --g {vcf} \
            --r {in_contig} \
            --out-gp-field \
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
    parser.add_argument('--local', action='store_true')
    parser.add_argument('--reference', type=str, default='GRCh38', choices=['GRCh37', 'GRCh38'])
    parser.add_argument('--cpu', type=int, default=8)
    parser.add_argument('--memory', type=str, default='standard', choices=['lowmem', 'standard', 'highmem'])
    parser.add_argument('--storage', type=int, default=50)
    parser.add_argument('--threads', type=int, default=16)
    parser.add_argument('--out-dir', required=True)

    args = parser.parse_args()

    if args.local:
        backend = hb.LocalBackend()
    else:
        backend = hb.ServiceBackend()

    phasing = hb.Batch(backend=backend,
                       name='genotype-imputation')

    vcf_paths = pd.read_csv(args.input_vcfs, sep='\t', header=None)

    for index, row in vcf_paths.iterrows():
        vcf = row[0]
        in_vcf = phasing.read_input(vcf)
        vcf_name = ntpath.basename(vcf)
        if vcf_name.endswith('.gz'):
            file_no_ext = vcf_name[:-7]
        elif vcf_name.endswith('.bgz'):
            file_no_ext = vcf_name[:-8]
        else:
            file_no_ext = vcf_name[:-4]

        for i in range(1, 24):
            chrom = f'chr{i}' if args.reference == 'GRCh38' else i

            # AFTER PHASING THE 1KG+HGDP DATA, CONVERT THE FILES TO IMP5 THEN CHANGE THE PATHS BELOW !!!
            ref_imp_chrom_file = f'gs://path/to/reference_chr{i}.imp5'
            ref_imp_chrom_file_idx = f'gs://path/to/reference_chr{i}.imp5.idx'
            ref_chrom_files = phasing.read_input_group(**{'imp5': ref_imp_chrom_file,
                                                          'imp5.idx': ref_imp_chrom_file_idx})

            imputation(b=phasing, vcf=in_vcf, vcf_filename_no_ext=file_no_ext, ref_imp5=ref_chrom_files,
                       reference=args.reference, contig=chrom, cpu=args.cpu, memory=args.memory,
                       storage=args.storage, threads=args.threads, out_dir=args.out_dir)


if __name__ == '__main__':
    main()
