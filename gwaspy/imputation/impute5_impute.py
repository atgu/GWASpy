__author__ = 'Lindo Nkambule'

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


def impute5_imputation(
        batch: hb.Batch = None,
        input_path: str = None,
        reference_path: str = None,
        output_filename: str = None,
        output_path: str = None):

    def imputation(
            b: hb.batch.Batch,
            vcf: hb.ResourceGroup = None,
            reference_vcf: hb.ResourceGroup = None,
            region: str = None,
            buffer_region: str = None,
            ncpu: int = 8,
            memory: str = 'highmem',
            storage: int = None,
            img: str = 'docker.io/lindonkambule/gwaspy_phase_impute:latest',
    ) -> Job:

        j = b.new_job(name=f'impute: {region}')
        chrom = region.split(":")[0]

        j.declare_resource_group(
            imputed_chunk={
                'chunk.bcf': '{root}.bcf',
                'chunk.bcf.csi': '{root}.bcf.csi'
            }
        )

        j.cpu(ncpu)
        j.memory(memory)
        j.regions(['us-central1'])
        j.storage(f'{storage}Gi')
        j.image(img)

        j.command(f"""
                    impute5_v1.2.0_static \
                    --h {reference_vcf['vcf']} \
                    --m /root/gwaspy/resources/maps/b38/{chrom}.b38.gmap.gz \
                    --g {vcf['vcf']} \
                    --r {region} \
                    --buffer-region {buffer_region} \
                    --o {j.imputed_chunk['chunk.bcf']} \
                    --threads {ncpu}
                    """
                  )

        return j

    def concatenate_imputed_chunks(
            b: hb.batch.Batch,
            chunks_list: List[hb.ResourceGroup] = None,
            output_vcf_name: str = None,
            chrom: str = None,
            memory: str = 'standard',
            out_dir: str = None,
            storage: int = None,
            ncpu: int = 4,
            img: str = 'docker.io/lindonkambule/gwaspy_phase_impute:latest',
    ) -> Job:
        """Concatenate imputed chunks by chromosome"""
        j = b.new_job(name=f'concatenate_imputed: {chrom}')

        imputed_chunks_files = '\n'.join([f'{v["chunk.bcf"]}' for v in chunks_list])

        j.declare_resource_group(
            concatenated_chrom={
                'bcf': '{root}.bcf',
                'bcf.csi': '{root}.bcf.csi'
            }
        )

        j.image(img)
        j.cpu(ncpu)
        j.memory(memory)
        j.regions(['us-central1'])
        j.storage(f'{storage}Gi')
        j.command(f'echo "{imputed_chunks_files}" > list_concatenate.txt')
        j.command(f"""
                    bcftools concat -n -f list_concatenate.txt -o {j.concatenated_chrom['bcf']}
                    """
                  )

        j.command(f"""
                    bcftools index {j.concatenated_chrom['bcf']} \
                    --output {j.concatenated_chrom['bcf.csi']} \
                    --threads {ncpu-2}
                    """
                  )

        b.write_output(j.concatenated_chrom,
                       f'{out_dir}/imputation/{output_vcf_name}_{chrom}.impute5.imputed')

        return j

    for i in range(1, 23):
        # read chrom input files
        if "CNUMBER" in input_path:     # input VCF is already split by chromosome
            vcf_path = input_path.replace('CNUMBER', str(i))
            input_idx = f'{vcf_path}.tbi' if hfs.exists(f'{vcf_path}.tbi') else f'{vcf_path}.csi'
            chrom_vcf = batch.read_input_group(**{'vcf': vcf_path,
                                                  'index': input_idx})
        else:
            vcf_path = input_path
            input_idx = f'{vcf_path}.tbi' if hfs.exists(f'{vcf_path}.tbi') else f'{vcf_path}.csi'
            chrom_vcf = batch.read_input_group(**{'vcf': input_path,
                                                  'index': input_idx})

        if not hfs.exists(input_idx):
            raise SystemExit('GWASpy requires the input file to be indexed (.tbi or .csi). Found none, exiting')

        vcf_size = round(size(vcf_path))

        ref_chrom_path = reference_path.replace('CNUMBER', str(i))
        ref_idx = f'{ref_chrom_path}.tbi' if hfs.exists(f'{ref_chrom_path}.tbi') else f'{ref_chrom_path}.csi'
        ref_vcf = batch.read_input_group(**{'vcf': ref_chrom_path,
                                            'index': ref_idx})
        ref_size = round(size(ref_chrom_path))

        imputation_chunks = pd.read_csv(
            f'https://raw.githubusercontent.com/odelaneau/shapeit5/main/resources/chunks/b38/4cM/chunks_chr{i}.txt',
            sep='\t', header=None,
            names=['index', 'chrom', 'irg', 'org', 'col5', 'col6', 'col7', 'col8'])
        imputation_chunks = imputation_chunks[['irg', 'org']]
        imp_chunks = list(imputation_chunks.itertuples(index=False, name=None))
        # imp_chunks_no_buffer = imputation_chunks['org'].tolist()  # 4th column (with no buffer between chunks)

        # Impute genotypes
        imputed_chunks = [
            imputation(
                b=batch,
                vcf=chrom_vcf,
                reference_vcf=ref_vcf,
                region=imp_chunks[i][1],
                buffer_region=imp_chunks[i][0],
                storage=round(vcf_size + ref_size + 5)
            ).imputed_chunk
            for i in range(len(imp_chunks))
        ]

        # Concatenate imputed chunks
        concatenate_imputed_chunks(
            b=batch,
            chunks_list=imputed_chunks,
            output_vcf_name=output_filename,
            chrom=f'chr{i}',
            out_dir=output_path,
            storage=round(vcf_size + ref_size + 10)
        )

    batch.run()
