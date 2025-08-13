__author__ = 'Mary T. Yohannes, Toni Boltz, Lindo Nkambule'

import hailtop.batch as hb
import hailtop.fs as hfs
import os
import pandas as pd
from hailtop.batch.job import Job
from typing import List


def size(file: str):
    """
    Convert the size from bytes to GiB
    :param file: path to file, str
    :return: file size in GiB
    """
    file_info = hfs.stat(file)   # returns a named tuple
    size_gigs = file_info.size / (1024 * 1024 * 1024)

    return size_gigs


def glimpse_phase_impute(
        batch: hb.Batch = None,
        input_file: str = None,
        reference_path: str = None,
        chromosomes: str = "all",
        output_filename: str = None,
        output_path: str = None):

    # 1. Convert the reference panel into GLIMPSE2’s binary file format.
    def create_binary_ref(
            b: hb.batch.Batch,
            ref_bcf: hb.ResourceGroup = None,
            input_region: str = None,
            output_region: str = None,
            img: str = 'docker.io/lindonkambule/gwaspy_glimpse2:with-bcftools-and-updated-info-score',
            ncpu: int = 4,
            memory: str = 'standard',
            storage: int = None
    ) -> Job:
        """Convert the reference panel into GLIMPSE2’s binary file format"""
        j = b.new_job(name=f'create binary ref panel: {input_region}')
        chrom = input_region.split(":")[0]
        bin_reg_out = input_region.replace(":", "_").replace("-", "_")

        j.image(img)
        j.cpu(ncpu)
        j.memory(memory)
        j.regions(['us-central1'])
        j.storage(f'{storage}Gi')

        j.command(f"""
                    GLIMPSE2_split_reference --reference {ref_bcf['bcf']} \
                    --map /gwaspy/resources/maps/b38/{chrom}.b38.gmap.gz \
                    --input-region {input_region} \
                    --output-region {output_region} \
                    --threads {ncpu} \
                    --output ref_bin_out
                    mv ref_bin_out_{bin_reg_out}.bin {j.ref_bin}
                    """
                  )
        return j

    # 2. Impute and phase chunks
    def glimpse2_impute_phase(
            b: hb.batch.Batch,
            in_bams_list: List[hb.ResourceGroup] = None,
            in_gls_vcf: hb.ResourceGroup = None,
            ref_binary: hb.ResourceGroup = None,
            region: str = None,
            img: str = 'docker.io/lindonkambule/gwaspy_glimpse2:with-bcftools-and-updated-info-score',
            ncpu: int = 4,
            memory: str = 'standard',
            storage: int = None
    ) -> Job:
        """Impute and phase data"""
        j = b.new_job(name=f'GLIMPSE2_phase: {region}')

        j.image(img)
        j.cpu(ncpu)
        j.memory(memory)
        j.regions(['us-central1'])
        j.storage(f'{storage}Gi')

        j.declare_resource_group(
            imputed_phased={
                'bcf': '{root}.bcf',
                'bcf.csi': '{root}.bcf.csi'
            }
        )

        if in_bams_list is not None:
            bams_in = '\n'.join([f'{v["bam"]}' for v in in_bams_list])
            j.command(f'echo "{bams_in}" > bams.txt')
            j.command(f"""
                        GLIMPSE2_phase --bam-list bams.txt \
                        --reference {ref_binary} \
                        --ne 20000 \
                        --threads 4 \
                        --output {j.imputed_phased['bcf']}""")
        else:
            j.command(f"""
                        GLIMPSE2_phase --input-gl {in_gls_vcf['vcf']} \
                        --reference {ref_binary} \
                        --ne 20000 \
                        --threads 4 \
                        --output {j.imputed_phased['bcf']}""")

        return j

    # 3. Ligate imputed chunks
    def ligate_chunks(
            b: hb.batch.Batch = None,
            imputed_chunk_list: List[hb.ResourceGroup] = None,
            chrom: str = None,
            img: str = 'docker.io/lindonkambule/gwaspy_glimpse2:with-bcftools-and-updated-info-score',
            ncpu: int = 4,
            memory: str = 'highmem',
            storage: int = None,
            output_vcf_name: str = None,
            out_dir: str = None
    ) -> Job:
        """Ligate imputed and phased data"""
        j = b.new_job(name=f'GLIMPSE2_ligate: {chrom}')

        j.image(img)
        j.cpu(ncpu)
        j.memory(memory)
        j.regions(['us-central1'])
        j.storage(f'{storage}Gi')

        j.declare_resource_group(
            ligated={
                'bcf': '{root}.bcf',
                'bcf.csi': '{root}.bcf.csi'
            }
        )

        imp_chunk_bcf_names = '\n'.join([f'{v["bcf"]}' for v in imputed_chunk_list])
        j.command(f'echo "{imp_chunk_bcf_names}" > chunk_list.txt')

        j.command(f"""GLIMPSE2_ligate --input chunk_list.txt --output {j.ligated['bcf']} --threads {ncpu}""")

        b.write_output(j.ligated,
                       f'{out_dir}/glimpse/{output_vcf_name}_{chrom}.phased.imputed')

        return j

    chroms = chromosomes.replace(" ", "") # remove spaces if there are any
    chroms = [i for i in range(1, 23)] if chroms == "all" else chroms.split(",")

    input_ext = os.path.splitext(input_file)[1]

    for i in chroms:
        # reference
        ref_chrom_path = reference_path.replace('CNUMBER', str(i))
        ref_idx = f'{ref_chrom_path}.tbi' if hfs.exists(f'{ref_chrom_path}.tbi') else f'{ref_chrom_path}.csi'
        ref_vcf = batch.read_input_group(**{'bcf': ref_chrom_path,
                                            'index': ref_idx})
        ref_size = round(size(ref_chrom_path))

        # input is BAM files
        if input_ext in ['.txt', '.tsv', '.csv']:
            bams = pd.read_csv(input_file, sep='\t', header=None, names=['path'])
            bams = bams['path'].tolist()
            gls_vcf = None

            # read the BAMs
            bam_list = []
            input_size = 0
            for bam_file in bams:
                bam_idx = f'{bam_file}.bai' if hfs.exists(f'{bam_file}.bai') else f'{bam_file}.crai'
                sample_bam = batch.read_input_group(**{'bam': bam_file,
                                                       'idx': bam_idx})
                bam_list.append(sample_bam)
                input_size += size(bam_file)
        else: # input is a BCF/VCF file with genotype likelihoods
            gls_chrom_path = input_file.replace('CNUMBER', str(i))
            gls_idx = f'{gls_chrom_path}.tbi' if hfs.exists(f'{gls_chrom_path}.tbi') else f'{gls_chrom_path}.csi'
            gls_vcf = batch.read_input_group(**{'vcf': gls_chrom_path,
                                                'index': gls_idx})
            input_size = size(gls_chrom_path)
            bam_list = None

        # create reference panel binary format to speed things up
        chunks_file = pd.read_csv(
            f'https://raw.githubusercontent.com/odelaneau/shapeit5/main/resources/chunks/b38/4cM/chunks_chr{i}.txt',
            sep='\t', header=None,
            names=['index', 'chrom', 'irg', 'org', 'col5', 'col6', 'col7', 'col8'])
        regions = [(irg, org) for irg, org in zip(chunks_file['irg'], chunks_file['org'])]

        ref_bin_scattered = [
            create_binary_ref(
                b=batch,
                ref_bcf=ref_vcf,
                input_region=regions[reg][0],
                output_region=regions[reg][1],
                storage=round(ref_size*0.2)
            ).ref_bin
            for reg in range(len(regions))
        ]

        # impute and phase using GLIMPSE2
        imputed_phased_chunks = [
            glimpse2_impute_phase(
                b=batch,
                in_bams_list=bam_list,
                in_gls_vcf=gls_vcf,
                ref_binary=ref_bin_scattered[reg],
                region=regions[reg][1],
                storage=round(ref_size*0.2 + input_size) # binary format files are smaller than originals
            ).imputed_phased
            for reg in range(len(regions))
        ]

        ligate_chunks(
            b=batch,
            imputed_chunk_list=imputed_phased_chunks,
            chrom=f'chr{i}',
            storage=round(ref_size*2),   # imputed calls should be slightly bigger than ref panel
            output_vcf_name=output_filename,
            out_dir=output_path
        )

    batch.run()
