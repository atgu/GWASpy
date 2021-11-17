__author__ = 'Lindo Nkambule'

import hailtop.batch as hb
import hail as hl
import pandas as pd
from typing import Union
from gwaspy.utils.get_file_size import bytes_to_gb
from gwaspy.phasing.get_filebase import get_vcf_filebase


def aut_impute(b: hb.batch.Batch,
               vcf: hb.ResourceGroup = None,
               vcf_filename_no_ext: str = None,
               ref: hb.ResourceGroup = None,
               region: str = None,
               chromosome: str = None,
               buffer: int = 250,
               storage: int = None,
               memory: str = None,
               cpu: int = None,
               img: str = 'docker.io/lindonkambule/gwaspy:v1',
               out_dir: str = None):

    out_file_name = vcf_filename_no_ext + '.imputed.bcf'
    file_dir = vcf_filename_no_ext.split('.')[0]

    map_file = f'/shapeit4/maps/b38/{chromosome}.b38.gmap.gz'

    threads = cpu - 1

    impute = b.new_job(name=out_file_name)
    impute.cpu(cpu)
    impute.memory(memory)
    impute.storage(f'{storage}Gi')
    impute.image(img)

    cmd = f'''
        impute5_1.1.5_static \
            --h {ref.bcf} \
            --m {map_file} \
            --g {vcf.bcf} \
            --r {region} \
            --out-gp-field \
            --o {out_file_name} \
            --b {buffer} \
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
               vcf: hb.ResourceGroup = None,
               vcf_filename_no_ext: str = None,
               females_list: hb.ResourceFile = None,
               ref: hb.ResourceGroup = None,
               region: str = None,
               chromosome: str = None,
               buffer: int = 250,
               storage: int = None,
               memory: str = None,
               cpu: int = None,
               img: str = 'docker.io/lindonkambule/gwaspy:v1',
               out_dir: str = None):

    global cmd

    in_females = females_list

    out_file_name = vcf_filename_no_ext + '.imputed.bcf'
    file_dir = vcf_filename_no_ext.split('.')[0]

    threads = cpu - 1

    map_file = f'/shapeit4/maps/b38/{chromosome}.b38.gmap.gz'

    impute = b.new_job(name=out_file_name)
    impute.cpu(cpu)
    impute.memory(memory)
    impute.storage(f'{storage}Gi')
    impute.image(img)

    if (int(region.split(':')[1].split('-')[1]) <= 2781479) | (int(region.split(':')[1].split('-')[0]) >= 155701383):
        # run diploid imputation
        cmd = f'''
            impute5_1.1.5_static \
                --h {ref.bcf} \
                --m {map_file} \
                --g {vcf.bcf} \
                --r {region} \
                --out-gp-field \
                --o {out_file_name} \
                --b {buffer} \
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
                echo THIS CHUNK IS ONLY IN THE NON-PAR REGION, SO WE WILL SPLIT IT BY SEX BEFORE RUNNING IMPUTATION
                bcftools view -S {in_females} {vcf.bcf} --output-type b --output females.bcf
                bcftools view -S ^{in_females} {vcf.bcf} --output-type b --output males.bcf
                echo INDEXING THE FILES BEFORE RUNNING IMPUTATION
                bcftools index females.bcf
                bcftools index males.bcf    
        '''
        impute.command(cmd_split)

        # (2) run imputation separately for females and males
        cmd_impute = f'''
                impute5_1.1.5_static --h {ref.bcf} --m {map_file} --g females.bcf --r {region} --out-gp-field \
                    --o females.imputed.bcf --b {buffer} --threads {threads}
                impute5_1.1.5_static --h {ref.bcf} --m {map_file} --g males.bcf --r {region} --out-gp-field \
                    --o males.imputed.bcf --b {buffer} --threads {threads} -- haploid
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
                    bcftools view {vcf.bcf} --regions chrX:10001-2781479 --output-type b --output par.bcf
                    bcftools view par.bcf | grep -v "#" | head -n1
                    echo SPLITTING OUT THE NON-PAR REGION
                    bcftools view {vcf.bcf} --regions chrX:2781479-155701382 --output-type b --output nonpar.bcf
                    echo SPLITTING THE NON-PAR REGION BY SEX
                    bcftools view -S {in_females} nonpar.bcf --output-type b --output females.bcf
                    bcftools view females.bcf | grep -v "#" | head -n1
                    bcftools view -S ^{in_females} nonpar.bcf --output-type b --output males.bcf
                    bcftools view males.bcf | grep -v "#" | head -n1
                    echo INDEXING THE FILES BEFORE RUNNING IMPUTATION
                    bcftools index par.bcf
                    bcftools index females.bcf
                    bcftools index males.bcf
                    
            '''

        else:
            cmd_split = f'''
                    echo SPLITTING OUT THE PAR2 REGION
                    bcftools view {vcf.bcf} --regions chrX:155701383-156030895 --output-type b --output par.bcf
                    bcftools view par.bcf | grep -v "#" | head -n1
                    echo SPLITTING OUT THE NON-PAR REGION
                    bcftools view {vcf.bcf} --regions chrX:2781479-155701382 --output-type b --output nonpar.bcf
                    echo SPLITTING THE NON-PAR REGION BY SEX
                    bcftools view -S {in_females} nonpar.bcf --output-type b --output females.bcf
                    bcftools view females.bcf | grep -v "#" | head -n1
                    bcftools view -S ^{in_females} nonpar.bcf --output-type b --output males.bcf
                    bcftools view males.bcf | grep -v "#" | head -n1
                    echo INDEXING THE FILES BEFORE RUNNING IMPUTATION
                    bcftools index par.bcf
                    bcftools index females.bcf
                    bcftools index males.bcf
            '''

        impute.command(cmd_split)

        # (2) run imputation separately for par, females, and males
        cmd_impute = f'''
                impute5_1.1.5_static --h {ref.bcf} --m {map_file} --g par.bcf --r {region} --out-gp-field \
                    --o par.imputed.bcf --b {buffer} --threads {threads}
                impute5_1.1.5_static --h {ref.bcf} --m {map_file} --g females.bcf --r {region} --out-gp-field \
                    --o females.imputed.bcf --b {buffer} --threads {threads}
                impute5_1.1.5_static --h {ref.bcf} --m {map_file} --g males.bcf --r {region} --out-gp-field \
                    --o males.imputed.bcf --b {buffer} --threads {threads} -- haploid
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
               n_samples: int = None,
               n_panel_samples: int = 4099,
               memory: str = 'highmem',
               buffer_region: int = 250,
               out_dir: str = None):

    print('RUNNING IMPUTATION')
    impute_b = hb.Batch(backend=backend, name=f'impute-phased-chunks')

    vcf_paths = pd.read_csv(input_vcfs, sep='\t', header=None)

    # use regions file to update the regions for imputation so that there's no overlaps like in phasing
    regions = pd.read_csv(f'{out_dir}/GWASpy/Phasing/gwaspy.refscatter.bed', delim_whitespace=True,
                          names=['chrom', 'start', 'end'])
    chroms_dfs = []

    for chrom, df_group in regions.groupby('chrom'):
        # print(df_group.loc[df_group.index[0], 'end'])
        df_group.loc[df_group.index[0], 'stop'] = df_group.loc[df_group.index[0], 'end']

        for i in range(1, len(df_group)):
            df_group.loc[df_group.index[i], 'stop'] = df_group.loc[df_group.index[i - 1], 'end'] + 1

        df_group['stop'] = df_group['stop'].astype(int)

        # add index column
        df_group['ind'] = df_group.index

        # update the first line to start at 1
        df_group.loc[df_group.index[0], 'stop'] = 1

        # combine the chromosome, start, and end positions into one
        df_group['reg'] = df_group['chrom'].astype(str) + ":" + df_group['stop'].astype(str) + "-" + df_group[
            'end'].astype(str)

        # select only the two needed columns
        regions_to_import_group = df_group[['reg', 'ind']]

        chroms_dfs.append(regions_to_import_group)

    regions_to_import = pd.concat(chroms_dfs, axis=0)
    regions_to_import = regions_to_import.sort_values('ind')
    regions_to_import.to_csv(f'{out_dir}/GWASpy/Imputation/imputation.regions', sep='\t', header=False,
                             index=False)

    regions_dict = pd.Series(regions_to_import.reg.values, index=regions_to_import.ind).to_dict()

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
            ref_size = bytes_to_gb(ref_bcf)
            ref = impute_b.read_input_group(**{'bcf': ref_bcf,
                                               'bcf.csi': f'{ref_bcf}.csi'})

            phased_bcf = f'{out_dir}/GWASpy/Phasing/{vcf_filebase}/phased_merged/{vcf_filebase}.{chrom}.phased.bcf'
            in_vcf = impute_b.read_input_group(**{'bcf': phased_bcf,
                                                  'bcf.csi': f'{phased_bcf}.csi'})
            vcf_size = bytes_to_gb(vcf)

            disk_size = int(round(10.0 + 3.0 * vcf_size + ((1.0 + 2.0 * n_samples/n_panel_samples) * ref_size)))
            job_memory = memory
            job_cpu = 16 if job_memory == 'highmem' else 8

            for file in phased_vcfs_chunks:
                f = file['path']
                vcf_basename = get_vcf_filebase(f)
                file_index = int(vcf_basename.split('.')[-3])
                file_region = regions_dict[file_index]
                map_chrom = file_region.split(':')[0]

                if map_chrom == chrom:
                    if chrom == 'chrX':
                        females_in = impute_b.read_input(females_file)

                        sex_impute(b=impute_b, vcf=in_vcf, females_list=females_in, vcf_filename_no_ext=vcf_basename,
                                   ref=ref, region=file_region, chromosome=chrom, buffer=buffer_region,
                                   storage=disk_size, memory=job_memory, cpu=job_cpu, out_dir=out_dir)
                    else:
                        aut_impute(b=impute_b, vcf=in_vcf, vcf_filename_no_ext=vcf_basename, ref=ref,
                                   region=file_region, chromosome=chrom, buffer=buffer_region, storage=disk_size,
                                   memory=job_memory, cpu=job_cpu, out_dir=out_dir)

    impute_b.run()

