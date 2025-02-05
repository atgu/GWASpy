__author__ = 'Mary T. Yohannes, Toni Boltz, Lindo Nkambule'

import hailtop.batch as hb
import hailtop.fs as hfs
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
        bam_files: str = None,
        reference_path: str = None,
        chromosomes: str = "all",
        output_filename: str = None,
        output_path: str = None):

    # 1. Perform a basic QC step by keeping only SNPs and remove multi-allelic records from the reference panel
    def qc_ref(
            b: hb.batch.Batch,
            ref_bcf: hb.ResourceGroup = None,
            chrom: str = None,
            ncpu: int = 4,
            memory: str = 'standard',
            storage: int = None,
            img: str = 'docker.io/lindonkambule/gwaspy_phase_impute:latest'
    ) -> Job:
        j = b.new_job(name=f'reference panel basic qc: {chrom}')  # define job

        j.image(img)
        j.cpu(ncpu)
        j.memory(memory)
        j.regions(['us-central1'])
        j.storage(f'{storage}Gi')

        j.declare_resource_group(
            qced_ref={
                'bcf': '{root}.bcf',
                'bcf.csi': '{root}.bcf.csi'
            }
        )

        j.command(f"""bcftools norm -m -any {ref_bcf['vcf']} -Ou --threads {ncpu} | bcftools view -m 2 -M 2 -v snps --threads {ncpu} -Ou -o {j.qced_ref['bcf']}""")
        j.command(f"""bcftools index {j.qced_ref['bcf']} --output {j.qced_ref['bcf.csi']} --threads {ncpu}""")

        return j

    # 2. Extracting variable sites in the reference panel (using QC'ed file from step 1)
    def extract_ref_sites(
            b: hb.batch.Batch,
            ref_bcf: hb.ResourceGroup = None,
            chrom: str = None,
            ncpu: int = 4,
            memory: str = 'standard',
            storage: int = None,
            img: str = 'docker.io/lindonkambule/gwaspy_glimpse2:with-bcftools-and-updated-info-score'
    ) -> Job:
        j = b.new_job(name=f'extract reference panel sites: {chrom}')  # define job

        j.image(img)
        j.cpu(ncpu)
        j.memory(memory)
        j.regions(['us-central1'])
        j.storage(f'{storage}Gi')

        j.declare_resource_group(
            ref_sites={
                'bcf': '{root}.bcf',
                'bcf.csi': '{root}.vcf.gz.csi',
                'tsv.gz': '{root}.tsv.gz',
                'tsv.gz.tbi': '{root}.tsv.gz.tbi'
            }
        )

        # extract sites in the reference panel, drop the genotypes
        j.command(f"""bcftools view -G -Ou -o {j.ref_sites['bcf']} {ref_bcf['bcf']} """)
        j.command(f"""bcftools index {j.ref_sites['bcf']} --output {j.ref_sites['bcf.csi']} --threads {ncpu}""")

        j.command(f"""bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {j.ref_sites['bcf']} | bgzip -c > {j.ref_sites['tsv.gz']}""")
        j.command(f"""tabix -s1 -b2 -e2 {j.ref_sites['tsv.gz']}""")

        return j

    # Step 3) Computing genotype likelihoods (GLs) for a single sample at specific positions
    def compute_gls(
            b: hb.batch.Batch,
            sample_bam_path: str = None,
            ref_sites: hb.ResourceGroup = None,
            ref_genome_fasta: hb.ResourceGroup = None,
            sample_id: str = None,
            chrom: str = None,
            ncpu: int = 4,
            memory: str = 'standard',
            storage: int = None,
            img: str = 'docker.io/lindonkambule/gwaspy_phase_impute:latest'
    ) -> Job:
        j = b.new_job(name=f'compute GLs:{sample_id}-{chrom}')

        bam_idx = f'{sample_bam_path}.bai' if hfs.exists(f'{sample_bam_path}.bai') else f'{sample_bam_path}.crai'
        sample_bam = b.read_input_group(**{'bam': sample_bam_path,
                                           'idx': bam_idx})

        j.image(img)
        j.cpu(ncpu)
        j.memory(memory)
        j.regions(['us-central1'])
        j.storage(f'{storage}Gi')

        j.declare_resource_group(
            gl_vcf={
                'bcf': '{root}.bcf',
                'bcf.csi': '{root}.bcf.csi'
            }
        )

        # compute GLs
        j.command(f"""bcftools mpileup -f {ref_genome_fasta['fasta']} -I -E -a 'FORMAT/DP' -T {ref_sites['bcf']} -r {chrom} {sample_bam['bam']} -Ou | bcftools call -Aim -C alleles -T {ref_sites['tsv.gz']} -Ou -o {j.gl_vcf['bcf']}""")

        j.command(f"""bcftools index {j.gl_vcf['bcf']} --output {j.gl_vcf['bcf.csi']} --threads {ncpu}""")

        return j

    # 4. Merging GLs across multiple BAMS per CHROM
    def merge_gl(
            b: hb.batch.Batch,
            gl_samples_list: List[hb.ResourceGroup] = None,
            chrom: str = None,
            ncpu: int = 4,
            memory: str = 'standard',
            storage: int = None,
            img: str = 'docker.io/lindonkambule/gwaspy_phase_impute:latest'
    ) -> Job:
        """Generate a single file containing genotype likelihoods for all target samples for a particular chromosome"""
        j = b.new_job(name=f'merge GLs-{chrom}')

        j.image(img)
        j.cpu(ncpu)
        j.memory(memory)
        j.regions(['us-central1'])
        j.storage(f'{storage}Gi')

        j.declare_resource_group(
            merged={
                'bcf': '{root}.bcf',
                'bcf.csi': '{root}.bcf.csi'
            }
        )

        gl_file_names = '\n'.join([f'{v["bcf"]}' for v in gl_samples_list])
        j.command(f'echo "{gl_file_names}" > list_merge.txt')

        j.command(f"""bcftools merge -m none -r {chrom} -Ou -o {j.merged['bcf']} -l list_merge.txt""")
        j.command(f"""bcftools index {j.merged['bcf']} --output {j.merged['bcf.csi']} --threads {ncpu}""")

        return j

    # 5. Convert the reference panel into GLIMPSE2’s binary file format.
    # Uses QC'ed reference files from step 1
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
                    mv ref_bin_out.bin {j.ref_bin}
                    """
                  )
        return j

    # 6. Impute and phase the chunks that are in a binary format
    # use the genotype likelihoods using --input-gl option instead of --bam-list
    def glimpse2_impute_phase(
            b: hb.batch.Batch,
            gl_vcf: hb.ResourceGroup = None,
            ref_binary: hb.ResourceGroup = None,
            region: str = None,
            img: str = 'docker.io/lindonkambule/gwaspy_glimpse2:with-bcftools-and-updated-info-score',
            ncpu: int = 4,
            memory: str = 'highmem',
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

        j.command(f"""
                    GLIMPSE2_phase --input-gl {gl_vcf['bcf']} \
                    --reference {ref_binary} \
                    --ne 20000 \
                    --threads 4 \
                    --output {j.imputed_phased['bcf']}""")

        j.command(f"""bcftools index {j.imputed_phased['bcf']} --output {j.imputed_phased['bcf.csi']} --threads {ncpu}""")

        return j

    # 7. Ligate imputed chunks
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
        j.command(f"""bcftools index {j.ligated['bcf']} --output {j.ligated['bcf.csi']} --threads {ncpu}""")

        b.write_output(j.ligated,
                       f'{out_dir}/glimpse/{output_vcf_name}_{chrom}.phased.imputed')

        return j

    ref_genome = batch.read_input_group(
        fasta='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta',
        idx='gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai'
    )
    ref_genome_size = round(size('gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta'))

    bams = pd.read_csv(bam_files, sep='\t', header=None, names=['sample', 'path'])
    bam_list = [(sample_id, bam_path) for sample_id, bam_path in zip(bams['sample'], bams['path'])]
    bam_sizes = [size(bam_list[s][1]) for s in range(len(bam_list))]

    chroms = chromosomes.replace(" ", "") # remove spaces if there are any
    chroms = [i for i in range(1, 23)] if chroms == "all" else chroms.split(",")

    for i in chroms:
        ref_chrom_path = reference_path.replace('CNUMBER', str(i))
        ref_idx = f'{ref_chrom_path}.tbi' if hfs.exists(f'{ref_chrom_path}.tbi') else f'{ref_chrom_path}.csi'
        ref_vcf = batch.read_input_group(**{'vcf': ref_chrom_path,
                                            'index': ref_idx})
        ref_size = round(size(ref_chrom_path))

        qced_ref = qc_ref(
            b=batch,
            ref_bcf=ref_vcf,
            chrom=f'chr{i}',
            storage=ref_size+2
        ).qced_ref

        # extract reference panel sites
        ref_panel_sites = extract_ref_sites(
            b=batch,
            ref_bcf=qced_ref,
            chrom=f'chr{i}',
            storage=ref_size+2
        ).ref_sites

        # compules GLs for each BAM
        bams_gls = [
            compute_gls(
                b=batch,
                sample_bam_path=bam_list[s][1],
                sample_id=bam_list[s][0],
                ref_sites=ref_panel_sites,
                ref_genome_fasta=ref_genome,
                chrom=f'chr{i}',
                storage=round(ref_genome_size+bam_sizes[s]+ref_size*0.3+5)
            ).gl_vcf
            for s in range(len(bam_list))
        ]

        # merge GLs into one file
        bams_merged_gls = merge_gl(
            b=batch,
            gl_samples_list=bams_gls,
            chrom=f'chr{i}',
            storage=round(sum(bam_sizes)*0.3)
        ).merged

        # create reference panel binary format to speed things up
        chunks_file = pd.read_csv(
            f'https://raw.githubusercontent.com/odelaneau/shapeit5/main/resources/chunks/b38/4cM/chunks_chr{i}.txt',
            sep='\t', header=None,
            names=['index', 'chrom', 'irg', 'org', 'col5', 'col6', 'col7', 'col8'])
        regions = [(irg, org) for irg, org in zip(chunks_file['irg'], chunks_file['org'])]

        ref_bin_scattered = [
            create_binary_ref(
                b=batch,
                ref_bcf=qced_ref,
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
                gl_vcf=bams_merged_gls,
                ref_binary=ref_bin_scattered[reg],
                region=regions[reg][1],
                storage=round(ref_size*0.2 + sum(bam_sizes)*0.3) # binary format and gl files are smaller than originals
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
