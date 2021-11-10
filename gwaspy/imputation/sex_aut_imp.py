__author__ = 'Lindo Nkambule'

import hailtop.batch as hb
from typing import Union
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

    global cmd
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

    if int(region.split(':')[1].split('-')[1]) < 2781479:
        # run diploid imputation
        print('This chunk is in the PAR1 region, we will run diploid imputation')
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

    elif int(region.split(':')[1].split('-')[0]) >= 155701383:
        # run diploid imputation
        print('This chunk is in the PAR2 region, we will run diploid imputation')
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

    elif (int(region.split(':')[1].split('-')[0]) >= 2781479) & (
            int(region.split(':')[1].split('-')[1]) <= 155701382):
        print('This chunk is only in the NON-PAR region, so we will split it by sex before running imputation')
        print(region)
        # (1) split by sex NB: THE USER SHOULD SUPPLY A FILE CONTAINING EITHER ONLY FEMALES OR MALE SAMPLE IDS
        # WITH ONE SAMPLE ID PER LINE
        # e.g. subset the VCF to only INCLUDE samples in females.txt
        # bcftools view -S females.txt hgdp.tgp.gwaspy.merged.19.vcf.gz > females.vcf.gz
        # e.g. subset the VCF to EXCLUDE samples in females.txt, just put a "^" sign before the file to exclude
        # bcftools view -S ^females.txt hgdp.tgp.gwaspy.merged.19.vcf.gz > males.vcf.gz

        # (2) run imputation separately for females and males
        cmd = f'''
                impute5_1.1.5_static --h {ref.bcf} --m {map_file} --g females.vcf.gz --r {region} --out-gp-field \
                --o {output_file_name} --threads {threads}
                impute5_1.1.5_static --h {ref.bcf} --m {map_file} --g males.vcf.gz --r {region} --out-gp-field \
                --o {output_file_name} --threads {threads} -- haploid
            '''
        impute.command(cmd)

        # (3) merge the imputed files back together into one chunk
        # bcftools index females.imputed.vcf.gz
        # bcftools index males.imputed.vcf.gz
        # echo -e "females.imputed.vcf.gz\nmales.imputed.vcf.gz" > merge.txt
        # bcftools merge -l merge.txt -Oz -o merge.vcf.gz
    else:
        # split the chunk into PAR AND non-PAR, then split non-PAR data by sex
        print('This chunk is in both PAR and NON-PAR regions, so we will split it by region, then the NON-PAR by sex')
        print(region)

    impute.command(cmd)

    impute.command(f'mv {output_file_name} {impute.ofile}')
    b.write_output(impute.ofile, f'{out_dir}/GWASpy/Imputation/{vcf_filename_no_ext}/{output_file_name}')