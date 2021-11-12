__author__ = 'Lindo Nkambule'

import hailtop.batch as hb
import hail as hl
import pandas as pd
from typing import Union
from gwaspy.utils.get_file_size import bytes_to_gb
from gwaspy.phasing.get_filebase import get_vcf_filebase


def aut_impute(b: hb.batch.Batch,
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

    in_vcf = b.read_input_group(**{'bcf': vcf,
                                'bcf.csi': f'{vcf}.csi'})
    vcf_size = bytes_to_gb(vcf)

    out_file_name = vcf_filename_no_ext + '.imputed.vcf.gz'
    file_dir = vcf_filename_no_ext.split('.')[0]

    disk_size = ref_size + (vcf_size * 4)

    map_file = f'/shapeit4/maps/b38/{chromosome}.b38.gmap.gz'

    impute = b.new_job(name=out_file_name)
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
            --o {out_file_name} \
            --threads {threads}
    '''

    impute.command(cmd)
    # index file to use when merging
    impute.command(f'bcftools index {out_file_name}')

    impute.command(f'mv {out_file_name} {impute.ofile}')
    impute.command(f'mv {out_file_name}.csi {impute.idx}')
    b.write_output(impute.ofile, f'{out_dir}/GWASpy/Imputation/{file_dir}/imputed_chunks/{out_file_name}')
    b.write_output(impute.idx, f'{out_dir}/GWASpy/Imputation/{file_dir}/imputed_chunks/{out_file_name}.csi')


def sex_impute(b: hb.batch.Batch,
               vcf: str = None,
               females_list: str = None,
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
    in_females = b.read_input(females_list)
    vcf_size = bytes_to_gb(vcf)

    out_file_name = vcf_filename_no_ext + '.imputed.bcf'
    file_dir = vcf_filename_no_ext.split('.')[0]

    disk_size = ref_size + (vcf_size * 4)

    map_file = f'/shapeit4/maps/b38/{chromosome}.b38.gmap.gz'

    impute = b.new_job(name=out_file_name)
    impute.cpu(cpu)
    impute.memory(memory)
    impute.storage(f'{disk_size}Gi')
    impute.image(img)

    if (int(region.split(':')[1].split('-')[1]) <= 2781479) | (int(region.split(':')[1].split('-')[0]) >= 155701383):
        # run diploid imputation
        cmd = f'''
            impute5_1.1.5_static \
                --h {ref.bcf} \
                --m {map_file} \
                --g {in_vcf} \
                --r {region} \
                --out-gp-field \
                --o {out_file_name} \
                --threads {threads}
        '''

        impute.command(cmd)
        # index file to use when merging
        impute.command(f'bcftools index {out_file_name}')

        impute.command(f'mv {out_file_name} {impute.ofile}')
        impute.command(f'mv {out_file_name}.csi {impute.idx}')
        b.write_output(impute.ofile, f'{out_dir}/GWASpy/Imputation/{file_dir}/imputed_chunks/{out_file_name}')
        b.write_output(impute.idx, f'{out_dir}/GWASpy/Imputation/{file_dir}/imputed_chunks/{out_file_name}.csi')

    elif (int(region.split(':')[1].split('-')[0]) >= 2781479) & (
            int(region.split(':')[1].split('-')[1]) <= 155701382):
        # (1) split by sex NB: THE USER SHOULD SUPPLY A FILE CONTAINING EITHER ONLY FEMALES OR MALE SAMPLE IDS
        # WITH ONE SAMPLE ID PER LINE
        cmd_split = f'''
                echo This chunk is only in the NON-PAR region, so we will split it by sex before running imputation
                bcftools view -S {in_females} {in_vcf} --output-type b --output females.bcf
                bcftools view -S ^{in_females} {in_vcf} --output-type b --output males.bcf
        '''
        impute.command(cmd_split)

        # (2) run imputation separately for females and males
        cmd_impute = f'''
                impute5_1.1.5_static --h {ref.bcf} --m {map_file} --g females.bcf --r {region} --out-gp-field \
                    --o females.imputed.bcf --threads {threads}
                impute5_1.1.5_static --h {ref.bcf} --m {map_file} --g males.bcf --r {region} --out-gp-field \
                    --o males.imputed.bcf --threads {threads} -- haploid
        '''
        impute.command(cmd_impute)

        # (3) merge the imputed files back together into one chunk
        # (3a) index females and males files to be merged and create a file with these for bcftools
        cmd_index_chunks = f'''
                bcftools index females.imputed.bcf \
                bcftools index males.imputed.bcf \
                echo -e "females.imputed.bcf\nmales.imputed.bcf" > merge.txt
        '''
        impute.command(cmd_index_chunks)

        # (3b) merge the files into one
        cmd_merge = f'''
                bcftools merge --file-list merge.txt --output-type b --output {out_file_name}
        '''
        impute.command(cmd_merge)

        impute.command(f'bcftools index {out_file_name}')

        impute.command(f'mv {out_file_name} {impute.ofile}')
        impute.command(f'mv {out_file_name}.csi {impute.idx}')
        b.write_output(impute.ofile, f'{out_dir}/GWASpy/Imputation/{file_dir}/imputed_chunks/{out_file_name}')
        b.write_output(impute.idx, f'{out_dir}/GWASpy/Imputation/{file_dir}/imputed_chunks/{out_file_name}.csi')

    else:
        # (1) split the chunk into PAR AND non-PAR, then split non-PAR data by sex
        if (int(region.split(':')[1].split('-')[0]) <= 2781479) & (int(region.split(':')[1].split('-')[1]) >= 2781479):
            cmd_split = f'''
                    echo SPLITTING OUT THE PAR1 REGION
                    bcftools view {in_vcf} --regions chrX:10001-2781479 --output-type b --output par.bcf
                    echo SPLITTING OUT THE NON-PAR REGION
                    bcftools view {in_vcf} --regions chrX:2781479-155701382 --output-type b --output nonpar.bcf
                    echo SPLITTING THE NON-PAR REGION BY SEX
                    bcftools view -S {in_females} nonpar.bcf --output-type b --output females.bcf
                    bcftools view -S ^{in_females} nonpar.bcf --output-type b --output males.bcf
                    
            '''

        else:
            cmd_split = f'''
                    echo SPLITTING OUT THE PAR2 REGION
                    bcftools view {in_vcf} --regions chrX:155701383-156030895 --output-type b --output par.bcf
                    echo SPLITTING OUT THE NON-PAR REGION
                    bcftools view {in_vcf} --regions chrX:2781479-155701382 --output-type b --output nonpar.bcf
                    echo SPLITTING THE NON-PAR REGION BY SEX
                    bcftools view -S {in_females} nonpar.bcf --output-type b --output females.bcf
                    bcftools view -S ^{in_females} nonpar.bcf --output-type b --output males.bcf
            '''

        impute.command(cmd_split)

        # (2) run imputation separately for par, females, and males
        cmd_impute = f'''
                impute5_1.1.5_static --h {ref.bcf} --m {map_file} --g par.bcf --r {region} --out-gp-field \
                    --o par.imputed.bcf --threads {threads}
                impute5_1.1.5_static --h {ref.bcf} --m {map_file} --g females.bcf --r {region} --out-gp-field \
                    --o females.imputed.bcf --threads {threads}
                impute5_1.1.5_static --h {ref.bcf} --m {map_file} --g males.bcf --r {region} --out-gp-field \
                    --o males.imputed.bcf --threads {threads} -- haploid
        '''
        impute.command(cmd_impute)

        # (3) merge the imputed files back together into one chunk
        # (3a) index females and males files to be merged and create a file with these for bcftools
        cmd_index_chunks = f'''
                    bcftools index par.imputed.bcf \
                    bcftools index females.imputed.bcf \
                    bcftools index males.imputed.bcf \
                    echo -e "females.imputed.bcf\nmales.imputed.bcf" > merge_sex.txt
        '''
        impute.command(cmd_index_chunks)

        # (3b) merge the sex files into one
        cmd_merge_sex = f'''
                bcftools merge --file-list merge_sex.txt --output-type b --output nonpar.imputed.bcf
        '''
        impute.command(cmd_merge_sex)

        # (3c) merge the non-par and par regions together
        if (int(region.split(':')[1].split('-')[0]) <= 2781479) & (int(region.split(':')[1].split('-')[1]) >= 2781479):
            cmd_merge_regs = f'''
                    echo -e "par.imputed.bcf\nnonpar.imputed.bcf" > merge_regions.txt
                    bcftools concat --naive --file-list merge_regions.txt --output-type b --output {out_file_name}
            '''
        else:
            cmd_merge_regs = f'''
                    echo -e "nonpar.imputed.bcf\npar.imputed.bcf" > merge_regions.txt
                    bcftools concat --naive --file-list merge_regions.txt --output-type b --output {out_file_name}
            '''
        impute.command(cmd_merge_regs)

        impute.command(f'bcftools index {out_file_name}')

        impute.command(f'mv {out_file_name} {impute.ofile}')
        impute.command(f'mv {out_file_name}.csi {impute.idx}')
        b.write_output(impute.ofile, f'{out_dir}/GWASpy/Imputation/{file_dir}/imputed_chunks/{out_file_name}')
        b.write_output(impute.idx, f'{out_dir}/GWASpy/Imputation/{file_dir}/imputed_chunks/{out_file_name}.csi')


def run_impute(backend: Union[hb.ServiceBackend, hb.LocalBackend] = None,
               input_vcfs: str = None,
               females_file: str = None,
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

        phased_vcfs_chunks = hl.utils.hadoop_ls(f'{out_dir}/GWASpy/Phasing/{vcf_filebase}/phased_scatter/*.bcf')

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
                    if chrom == 'chrX':
                        sex_impute(b=impute_b, vcf=f, females_list=females_file, vcf_filename_no_ext=vcf_basename,
                                   ref=ref, ref_size=ref_size, region=file_region, chromosome=chrom, cpu=cpu,
                                   memory=memory, threads=threads, out_dir=out_dir)
                    else:
                        aut_impute(b=impute_b, vcf=f, vcf_filename_no_ext=vcf_basename, ref=ref, ref_size=ref_size,
                                   region=file_region, chromosome=chrom, cpu=cpu, memory=memory,
                                   threads=threads, out_dir=out_dir)

    impute_b.run()

