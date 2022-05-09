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
    b.write_output(impute.ofile, f'{out_dir}/GWASpy/{file_dir}/Imputation/imputed_chunks/{out_file_name}')
    b.write_output(impute.idx, f'{out_dir}/GWASpy/{file_dir}/Imputation/imputed_chunks/{out_file_name}.csi')


def sex_impute(b: hb.batch.Batch,
               vcf: hb.ResourceGroup = None,
               vcf_filename_no_ext: str = None,
               females_list: hb.ResourceFile = None,
               ref: hb.ResourceGroup = None,
               region: str = None,
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

    impute = b.new_job(name=out_file_name)
    impute.cpu(cpu)
    impute.memory(memory)
    impute.storage(f'{storage}Gi')
    impute.image(img)

    start = int(region.split(':')[1].split('-')[0])
    end = int(region.split(':')[1].split('-')[1])

    # A. PAR1 REGION ONLY
    if end <= 2781479:
        # run diploid imputation
        impute.command('echo THIS CHUNK IS IN PAR1 REGION, SO WE WILL RUN DIPLOID IMPUTATION')
        map_file = '/shapeit4/maps/b38/chrX_par1.b38.gmap.gz'

        if start < 10001:
            start_new = 10001
            end_new = end
        else:
            start_new = start
            end_new = end

        new_imp_region = f'chrX:{start_new}-{end_new}'

        cmd = f'''
            impute5_1.1.5_static \
                --h {ref.bcf} \
                --m {map_file} \
                --g {vcf.bcf} \
                --r {new_imp_region} \
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
        b.write_output(impute.ofile, f'{out_dir}/GWASpy/{file_dir}/Imputation/imputed_chunks/{out_file_name}')
        b.write_output(impute.idx, f'{out_dir}/GWASpy/{file_dir}/Imputation/imputed_chunks/{out_file_name}.csi')

    # B. PAR2 REGION ONLY
    elif start >= 155701383:
        # run diploid imputation
        impute.command('echo THIS CHUNK IS IN PAR2 REGION, SO WE WILL RUN DIPLOID IMPUTATION')
        map_file = '/shapeit4/maps/b38/chrX_par2.b38.gmap.gz'

        if end > 156030895:
            end_new = 156030895
            start_new = start
        else:
            start_new = start
            end_new = end

        new_imp_region = f'chrX:{start_new}-{end_new}'

        cmd = f'''
            impute5_1.1.5_static \
                --h {ref.bcf} \
                --m {map_file} \
                --g {vcf.bcf} \
                --r {new_imp_region} \
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
        b.write_output(impute.ofile, f'{out_dir}/GWASpy/{file_dir}/Imputation/imputed_chunks/{out_file_name}')
        b.write_output(impute.idx, f'{out_dir}/GWASpy/{file_dir}/Imputation/imputed_chunks/{out_file_name}.csi')

    # C. NON-PAR REGION ONLY
    elif (start >= 2781479) & (end <= 155701382):
        map_file = '/shapeit4/maps/b38/chrX.b38.gmap.gz'
        start_new = start
        end_new = end
        new_imp_region = f'chrX:{start_new}-{end_new}'

        # (1) split by sex NB: THE USER SHOULD SUPPLY A FILE CONTAINING EITHER ONLY FEMALES OR MALE SAMPLE IDS
        # WITH ONE SAMPLE ID PER LINE
        cmd_split = f'''
                echo THIS CHUNK IS ONLY IN THE NON-PAR REGION, SO WE WILL SPLIT IT BY SEX BEFORE RUNNING IMPUTATION
                echo GETTING THE ORDER OF SAMPLES
                bcftools query -l {vcf.bcf} > samples_order.txt
                echo SPLITTING SAMPLES BY SEX
                bcftools view -S {in_females} {vcf.bcf} --output-type b --output females.bcf
                bcftools view -S ^{in_females} {vcf.bcf} --output-type b --output males.bcf
                echo INDEXING THE FILES BEFORE RUNNING IMPUTATION
                bcftools index females.bcf
                bcftools index males.bcf
        '''
        impute.command(cmd_split)

        # (2) run imputation separately for females and males
        cmd_impute = f'''
                impute5_1.1.5_static --h {ref.bcf} --m {map_file} --g females.bcf --r {new_imp_region} --out-gp-field \
                    --o females.imputed.bcf --b {buffer} --threads {threads}
                impute5_1.1.5_static --h {ref.bcf} --m {map_file} --g males.bcf --r {new_imp_region} --out-gp-field \
                    --o males.imputed.bcf --b {buffer} --threads {threads} -- haploid
        '''
        impute.command(cmd_impute)

        # (3) merge the imputed files back together into one chunk
        # (3a) index females and males files to be merged and create a file with these for bcftools
        cmd_index_chunks = f'''
                bcftools index females.imputed.bcf
                bcftools index males.imputed.bcf
                rm females.bcf* males.bcf*
        '''
        impute.command(cmd_index_chunks)

        # (3b) sometimes there are duplicates (raw and imputed with flipped alleles) and this causes an error when
        # merging, so we remove any duplicates
        cmd_sort = f'''
                echo CHECKING FOR DUPLICATE VARIANTS AND REMOVING THEM
                bcftools norm -d any females.imputed.bcf --output-type b --output females.imputed.sorted.bcf
                rm females.imputed.bcf*
                bcftools norm -d any males.imputed.bcf --output-type b --output males.imputed.sorted.bcf
                rm males.imputed.bcf*
                bcftools index females.imputed.sorted.bcf
                bcftools index males.imputed.sorted.bcf
                echo -e "females.imputed.sorted.bcf\nmales.imputed.sorted.bcf" > merge.txt
        '''

        impute.command(cmd_sort)

        # (3c) merge the files into one
        cmd_merge = f'''
                echo MERGING THE MALES AND FEMALES FILE BACK TOGETHER
                bcftools merge --file-list merge.txt --output-type b --output merged.sex.bcf
        '''
        impute.command(cmd_merge)

        # (3d) reorder samples back to how they were initially
        # splitting and re-merging samples will result in change in order compared to the initial file
        # so we have to use the initial samples order for the concat step by chromosome of imputation
        cmd_order = f'''
                echo ORDERING SAMPLES TO HOW THEY WERE INITIALLY
                bcftools view -S samples_order.txt merged.sex.bcf --output-type b --output {out_file_name}
        '''
        impute.command(cmd_order)

        impute.command(f'bcftools index {out_file_name}')

        impute.command(f'mv {out_file_name} {impute.ofile}')
        impute.command(f'mv {out_file_name}.csi {impute.idx}')
        b.write_output(impute.ofile, f'{out_dir}/GWASpy/{file_dir}/Imputation/imputed_chunks/{out_file_name}')
        b.write_output(impute.idx, f'{out_dir}/GWASpy/{file_dir}/Imputation/imputed_chunks/{out_file_name}.csi')

    # D. MIXED REGIONS
    else:
        if (start <= 2781479) & (end >= 2781479):
            map_file_par = '/shapeit4/maps/b38/chrX_par1.b38.gmap.gz'
            par_region = f'chrX:{start}-{2781479}'
            non_par_region = f'chrX:{2781479}-{end}'

            # to be used in step 3d
            cmd_merge_regs = f'''
                    echo ORDERING PAR1 SAMPLES TO HOW THEY WERE INITIALLY
                    bcftools view -S samples_order.txt par.imputed.bcf --output-type b --output par.imputed.sort.bcf
                    bcftools index par.imputed.sort.bcf
                    rm par.imputed.bcf*
                    echo MERGING THE PAR1 AND NON-PAR FILES BACK TOGETHER
                    echo -e "par.imputed.sort.bcf\nnonpar.imputed.sort.bcf" > merge_regions.txt
                    bcftools concat --naive --file-list merge_regions.txt --output-type b --output {out_file_name}
            '''

        else:
            map_file_par = '/shapeit4/maps/b38/chrX_par2.b38.gmap.gz'
            if end > 156030895:
                par_region = f'chrX:{155701383}-{156030895}'
            else:
                par_region = f'chrX:{155701383}-{end}'

            non_par_region = f'chrX:{start}-{155701382}'

            # to be used in step 3c
            cmd_merge_regs = f'''
                    echo ORDERING PAR2 SAMPLES TO HOW THEY WERE INITIALLY
                    bcftools view -S samples_order.txt par.imputed.bcf --output-type b --output par.imputed.sort.bcf
                    bcftools index par.imputed.sort.bcf
                    rm par.imputed.bcf*
                    echo MERGING THE PAR2 AND NON-PAR FILES BACK TOGETHER
                    echo -e "nonpar.imputed.sort.bcf\npar.imputed.sort.bcf" > merge_regions.txt
                    bcftools concat --naive --file-list merge_regions.txt --output-type b --output {out_file_name}
            '''

        cmd_split = f'''
                echo THIS CHUNK IS IN PAR1/2 and NON-PAR REGION, SO WE WILL SPLIT THESE TWO REGIONS
                echo GETTING THE ORDER OF SAMPLES
                bcftools query -l {vcf.bcf} > samples_order.txt
                echo SPLITTING OUT THE PAR REGION
                bcftools view {vcf.bcf} --regions {par_region} --output-type b --output par.bcf
                echo SPLITTING OUT THE NON-PAR REGION
                bcftools view {vcf.bcf} --regions {non_par_region} --output-type b --output nonpar.bcf
                echo SPLITTING THE NON-PAR REGION BY SEX
                bcftools view -S {in_females} nonpar.bcf --output-type b --output females.bcf
                bcftools view -S ^{in_females} nonpar.bcf --output-type b --output males.bcf
                echo INDEXING THE FILES BEFORE RUNNING IMPUTATION
                bcftools index par.bcf
                bcftools index females.bcf
                bcftools index males.bcf
                rm {vcf.bcf}
                rm nonpar.bcf
        '''

        impute.command(cmd_split)

        # (2) run imputation separately for par, females, and males
        map_file = '/shapeit4/maps/b38/chrX.b38.gmap.gz'
        cmd_impute = f'''
                impute5_1.1.5_static --h {ref.bcf} --m {map_file_par} --g par.bcf --r {par_region} --out-gp-field \
                    --o par.imputed.bcf --b {buffer} --threads {threads}
                impute5_1.1.5_static --h {ref.bcf} --m {map_file} --g females.bcf --r {non_par_region} --out-gp-field \
                    --o females.imputed.bcf --b {buffer} --threads {threads}
                impute5_1.1.5_static --h {ref.bcf} --m {map_file} --g males.bcf --r {non_par_region} --out-gp-field \
                    --o males.imputed.bcf --b {buffer} --threads {threads} -- haploid
        '''
        impute.command(cmd_impute)

        # (3) merge the imputed files back together into one chunk
        # (3a) index females and males files to be merged and create a file with these for bcftools
        cmd_index_chunks = f'''
                    bcftools index par.imputed.bcf
                    bcftools index females.imputed.bcf
                    bcftools index males.imputed.bcf
        '''
        impute.command(cmd_index_chunks)

        # (3b) sometimes there are duplicates (raw and imputed with flipped alleles) and this causes an error when
        # merging, so we remove any duplicates
        # bcftools view Nepal_PTSD_GSA_Updated_May2021_qced.chrX.phased.bcf | grep "#"
        cmd_sort = f'''
                echo CHECKING FOR DUPLICATE VARIANTS AND REMOVING THEM
                bcftools norm -d any females.imputed.bcf --output-type b --output females.imputed.sorted.bcf
                rm females.imputed.bcf*
                bcftools norm -d any males.imputed.bcf --output-type b --output males.imputed.sorted.bcf
                rm males.imputed.bcf*
                bcftools index females.imputed.sorted.bcf
                bcftools index males.imputed.sorted.bcf
                echo -e "females.imputed.sorted.bcf\nmales.imputed.sorted.bcf" > merge_sex.txt
        '''
        impute.command(cmd_sort)

        # (3c) merge the sex files into one
        cmd_merge_sex = f'''
                echo MERGING THE MALES AND FEMALES FILE BACK TOGETHER
                bcftools merge --file-list merge_sex.txt --output-type b --output nonpar.imputed.bcf
                bcftools index nonpar.imputed.bcf
                echo ORDERING NON-PAR SAMPLES TO HOW THEY WERE INITIALLY
                bcftools view -S samples_order.txt nonpar.imputed.bcf --output-type b --output nonpar.imputed.sort.bcf
                bcftools index nonpar.imputed.sort.bcf
                rm nonpar.imputed.bcf* females.imputed.sorted.bcf* males.imputed.sorted.bcf
        '''
        impute.command(cmd_merge_sex)

        # (3d) merge the non-par and par regions together
        impute.command(cmd_merge_regs)

        impute.command(f'bcftools index {out_file_name}')

        impute.command(f'mv {out_file_name} {impute.ofile}')
        impute.command(f'mv {out_file_name}.csi {impute.idx}')
        b.write_output(impute.ofile, f'{out_dir}/GWASpy/{file_dir}/Imputation/imputed_chunks/{out_file_name}')
        b.write_output(impute.idx, f'{out_dir}/GWASpy/{file_dir}/Imputation/imputed_chunks/{out_file_name}.csi')


def run_impute(backend: Union[hb.ServiceBackend, hb.LocalBackend] = None,
               input_vcf: str = None,
               females_file: str = None,
               n_samples: int = None,
               n_panel_samples: int = 4099,
               phasing_software: str = None,
               memory: str = 'highmem',
               buffer_region: int = 250,
               out_dir: str = None):

    global phased_bcf
    print(f'\n1. IMPUTATION ON {input_vcf} PHASED CHUNKS\n')
    vcf_filebase = get_vcf_filebase(input_vcf)

    impute_b = hb.Batch(backend=backend, name=f'impute-phased-chunks-{vcf_filebase}')

    # use regions file to update the regions for imputation so that there's no overlaps like in phasing
    regions = pd.read_csv(f'{out_dir}/GWASpy/{vcf_filebase}/Phasing/refscatter.bed', delim_whitespace=True,
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
    regions_to_import.to_csv(f'{out_dir}/GWASpy/{vcf_filebase}/Imputation/imputation.regions', sep='\t', header=False,
                             index=False)

    regions_dict = pd.Series(regions_to_import.reg.values, index=regions_to_import.ind).to_dict()

    if phasing_software == 'shapeit':
        phased_vcfs_chunks = hl.utils.hadoop_ls(f'{out_dir}/GWASpy/{vcf_filebase}/Phasing/phased_scatter/*.shapeit.bcf')
    else:
        phased_vcfs_chunks = hl.utils.hadoop_ls(f'{out_dir}/GWASpy/{vcf_filebase}/Phasing/phased_scatter/*.eagle.bcf')

    for i in range(1, 24):
        if i == 23:
            chrom = 'chrX'
        else:
            chrom = f'chr{i}'

        ref_bcf = f'gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes/hgdp.tgp.gwaspy.merged.{chrom}.merged.bcf'
        ref_size = bytes_to_gb(ref_bcf)
        ref = impute_b.read_input_group(**{'bcf': ref_bcf,
                                           'bcf.csi': f'{ref_bcf}.csi'})

        # output is not always bcf
        phased_filename = f'{out_dir}/GWASpy/{vcf_filebase}/Phasing/phased_merged/{vcf_filebase}.{chrom}.phased.{phasing_software}'
        if hl.hadoop_exists(f'{phased_filename}.bcf'):
            phased_bcf = f'{phased_filename}.bcf'
        elif hl.hadoop_exists(f'{phased_filename}.vcf.gz'):
            phased_bcf = f'{phased_filename}.vcf.gz'

        in_vcf = impute_b.read_input_group(**{'bcf': phased_bcf,
                                              'bcf.csi': f'{phased_bcf}.csi'})
        vcf_size = bytes_to_gb(input_vcf)

        disk_size = int(round(10.0 + 3.0 * vcf_size + ((1.0 + 2.0 * n_samples/n_panel_samples) * ref_size)))
        job_memory = memory
        job_cpu = 16 if job_memory == 'highmem' else 8

        for file in phased_vcfs_chunks:
            f = file['path']
            vcf_basename = get_vcf_filebase(f)
            file_index = int(vcf_basename.split('.')[-3])
            file_region = regions_dict[file_index]
            map_chrom = file_region.split(':')[0]

            imp_out_filename = f'{vcf_basename}.imputed.bcf'
            # file_dir = vcf_basename.split('.')[0]
            output_filepath_name = f'{out_dir}/GWASpy/{vcf_filebase}/Imputation/imputed_chunks/{imp_out_filename}'

            if map_chrom == chrom:
                # check if imputed file already exists
                if hl.hadoop_exists(output_filepath_name):
                    continue

                else:
                    if chrom == 'chrX':
                        females_in = impute_b.read_input(females_file)

                        sex_impute(b=impute_b, vcf=in_vcf, females_list=females_in, vcf_filename_no_ext=vcf_basename,
                                   ref=ref, region=file_region, buffer=buffer_region,
                                   storage=disk_size, memory=job_memory, cpu=job_cpu, out_dir=out_dir)
                    else:
                        aut_impute(b=impute_b, vcf=in_vcf, vcf_filename_no_ext=vcf_basename, ref=ref,
                                   region=file_region, chromosome=chrom, buffer=buffer_region, storage=disk_size,
                                   memory=job_memory, cpu=job_cpu, out_dir=out_dir)

    impute_b.run()

