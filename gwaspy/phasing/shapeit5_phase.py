__author__ = 'Lindo Nkambule'

import hailtop.batch as hb
import hailtop.fs as hfs
import pandas as pd
from hailtop.batch.job import Job
from typing import List, Union, Optional


def size(file: str):
    """
    Convert the size from bytes to GiB
    :param file: path to file, str
    :return: file size in GiB
    """

    file_info = hfs.stat(file)   # returns a named tuple
    size_gigs = file_info.size / (1024 * 1024 * 1024)

    return size_gigs


def shapeit_phasing(
        batch: hb.Batch = None,
        input_path: str = None,
        reference_path: Optional[str] = None,
        fam_file: Optional[hb.ResourceFile] = None,
        data_type: str = 'array',
        output_filename: str = None,
        output_path: str = None):

    def phase_common(
            b: hb.batch.Batch = None,
            vcf: hb.ResourceGroup = None,
            ref_vf: Optional[hb.ResourceGroup] = None,
            maf: float = 0.001,
            region: str = None,
            pedigree: Optional[hb.ResourceFile] = None,
            ncpu: int = 8,
            memory: str = 'standard',
            storage: int = None,
            img: str = 'docker.io/lindonkambule/gwaspy_phase_impute:latest',
            output_vcf_name: str = None,
            out_dir: str = None
    ) -> Job:
        """Phase common variants by chunks"""
        j = b.new_job(name=f'phase_common: {region}')
        chrom = region.split(":")[0]

        j.declare_resource_group(
            phased_common_chunk={
                'bcf': '{root}.bcf',
                'bcf.csi': '{root}.bcf.csi',
                'log': '{root}.log'
            }
        )

        j.image(img)
        j.cpu(ncpu)
        j.memory(memory)
        j.regions(['us-central1'])
        j.storage(f'{storage}Gi')

        # phase common variants
        j.command(f"""
                    phase_common_static --input {vcf['vcf']} \
                    {f"--reference {ref_vf['vcf']}" if ref_vf else ''} \
                    --map /root/gwaspy/resources/maps/b38/{chrom}.b38.gmap.gz \
                    --output {j.phased_common_chunk['bcf']} \
                    --thread {ncpu-1} \
                    --log {j.phased_common_chunk['log']} \
                    --filter-maf {maf} \
                    {f'--pedigree {pedigree}' if pedigree else ''} \
                    --region {region}
                    """
                  )

        if out_dir:
            b.write_output(j.phased_common_chunk,
                           f'{out_dir}/phase_common/{output_vcf_name}_{chrom}.array.shapeit5_common')

        return j

    def ligate_common_chunks(
            b: hb.batch.Batch,
            common_variants_chunks_list: List[hb.ResourceGroup] = None,
            pedigree: Optional[hb.ResourceFile] = None,
            output_vcf_name: str = None,
            chrom: str = None,
            memory: str = 'standard',
            out_dir: str = None,
            storage: int = None,
            ncpu: int = 4,
            img: str = 'docker.io/lindonkambule/gwaspy_phase_impute:latest',
    ) -> Job:
        # requires a VCF/BCF with its index
        j = b.new_job(name=f'ligate_common: {chrom}')

        phased_common_chunks = '\n'.join([f'{v["bcf"]}' for v in common_variants_chunks_list])

        j.declare_resource_group(
            ligated_chrom={
                'bcf': '{root}.bcf',
                'bcf.csi': '{root}.bcf.csi',
                'log': '{root}.log'
            }
        )

        j.image(img)
        j.cpu(ncpu)
        j.memory(memory)
        j.regions(['us-central1'])
        j.storage(f'{storage}Gi')
        j.command(f'echo "{phased_common_chunks}" > common_chunks_list_ligate.txt')
        j.command(f"""
                    ligate_static --input common_chunks_list_ligate.txt \
                    {f'--pedigree {pedigree}' if pedigree else ''} \
                    --output {j.ligated_chrom['bcf']} \
                    --thread {ncpu-1} \
                    --log {j.ligated_chrom['log']} \
                    --index
                    """
                  )

        b.write_output(j.ligated_chrom,
                       f'{out_dir}/phase_common/{output_vcf_name}_{chrom}.shapeit5_common')

        return j

    def phase_rare(
            b: hb.batch.Batch = None,
            vcf: hb.ResourceGroup = None,
            scaffold_vcf: hb.ResourceGroup = None,
            scaffold_region: str = None,
            input_region: str = None,
            pedigree: Optional[hb.ResourceFile] = None,
            ncpu: int = 4,
            memory: str = 'standard',
            storage: int = None,
            img: str = 'docker.io/lindonkambule/gwaspy_phase_impute:latest',
    ) -> Job:
        """Phase common variants by chunks"""
        j = b.new_job(name=f'phase_rare: {scaffold_region}')
        chrom = scaffold_region.split(":")[0]

        j.declare_resource_group(
            phased_rare_chunk={
                'chunk.bcf': '{root}.bcf',
                'chunk.bcf.csi': '{root}.bcf.csi',
                'chunk.log': '{root}.log'
            }
        )

        j.image(img)
        j.cpu(ncpu)
        j.memory(memory)
        j.regions(['us-central1'])
        j.storage(f'{storage}Gi')

        j.command(f"""
                    phase_rare_static \
                    --input {vcf['vcf']} --input-region {input_region} \
                    --scaffold {scaffold_vcf['bcf']} --scaffold-region {scaffold_region} \
                    --map /root/gwaspy/resources/maps/b38/{chrom}.b38.gmap.gz \
                    {f'--pedigree {pedigree}' if pedigree else ''} \
                    --output {j.phased_rare_chunk['chunk.bcf']} \
                    --thread {ncpu-1} \
                    --log {j.phased_rare_chunk['chunk.log']}
                    """
                  )

        return j

    def concatenate_rare_chunks(
            b: hb.batch.Batch,
            rare_variants_chunks_list: List[hb.ResourceGroup] = None,
            output_vcf_name: str = None,
            chrom: str = None,
            memory: str = 'standard',
            out_dir: str = None,
            storage: int = None,
            ncpu: int = 4,
            img: str = 'docker.io/lindonkambule/gwaspy_phase_impute:latest',
    ) -> Job:
        j = b.new_job(name=f'concatenate_rare: {chrom}')

        phased_rare_chunks = '\n'.join([f'{v["chunk.bcf"]}' for v in rare_variants_chunks_list])

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
        j.command(f'echo "{phased_rare_chunks}" > list_concatenate.txt')
        j.command(f"""
                    bcftools concat -n -f list_concatenate.txt -o {j.concatenated_chrom['bcf']}
                    """
                  )

        j.command(f"""
                    bcftools index {j.concatenated_chrom['bcf']} \
                    --output {j.concatenated_chrom['bcf.csi']} \
                    --threads {ncpu-1}
                    """
                  )

        b.write_output(j.concatenated_chrom,
                       f'{out_dir}/phase_rare/{output_vcf_name}_{chrom}.full.shapeit5_rare')

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

        if reference_path:
            ref_chrom_path = reference_path.replace('CNUMBER', str(i))
            ref_idx = f'{ref_chrom_path}.tbi' if hfs.exists(f'{ref_chrom_path}.tbi') else f'{ref_chrom_path}.csi'
            ref_vcf = batch.read_input_group(**{'vcf': ref_chrom_path,
                                                'index': ref_idx})
            ref_size = round(size(ref_chrom_path))
        else:
            ref_vcf = None
            ref_size = 0

        if data_type == 'array':
            phase_common(
                b=batch,
                vcf=chrom_vcf,
                ref_vf=ref_vcf,
                maf=0.0,
                pedigree=fam_file,
                region=f'chr{i}',
                storage=round(vcf_size*1.5 + ref_size + 2),
                output_vcf_name=output_filename,
                out_dir=f'{output_path}/shapeit5'
            )
        else:   # WGS
            # 1A. Phase common chunks
            common_chunks_file = pd.read_csv(
                f'https://raw.githubusercontent.com/odelaneau/shapeit5/main/resources/chunks/b38/20cM/chunks_chr{i}.txt',
                sep='\t', header=None,
                names=['index', 'chrom', 'irg', 'org'])
            common_regions = common_chunks_file['irg'].values.tolist()
            common_chunks_phased = [
                phase_common(
                    b=batch,
                    vcf=chrom_vcf,
                    ref_vf=ref_vcf,
                    pedigree=fam_file,
                    region=common_regions[i],
                    storage=round(vcf_size*1.5 + ref_size + 2)
                ).phased_common_chunk
                for i in range(len(common_regions))
            ]

            # 1B. Ligate phased common chunks into one chromosome
            common_phased_ligated_scaffold = ligate_common_chunks(
                b=batch,
                common_variants_chunks_list=common_chunks_phased,
                pedigree=fam_file,
                output_vcf_name=output_filename,
                chrom=f'chr{i}',
                out_dir=f'{output_path}/shapeit5',
                storage=round(vcf_size*0.1)
            ).ligated_chrom

            # 2A. Phase rare chunks
            rare_chunks_file = pd.read_csv(
                f'https://raw.githubusercontent.com/odelaneau/shapeit5/main/resources/chunks/b38/4cM/chunks_chr{i}.txt',
                sep='\t', header=None,
                names=['index', 'chrom', 'irg', 'org', 'col5', 'col6', 'col7', 'col8'])
            # irg (3rd col) is SCAFFOLD_REG and org (4th col) is INPUT_REG
            # https://github.com/odelaneau/shapeit5/tree/main/resources/chunks/b38
            rare_regions = [(irg, org) for irg, org in zip(rare_chunks_file['irg'], rare_chunks_file['org'])]
            rare_chunks_phased = [
                phase_rare(
                    b=batch,
                    vcf=chrom_vcf,
                    scaffold_vcf=common_phased_ligated_scaffold,
                    scaffold_region=rare_regions[i][0],  # irg (3rd col in chunks)
                    input_region=rare_regions[i][1],  # org (4th col in chunks)
                    pedigree=fam_file,
                    storage=round(vcf_size*1.5*1.5)  # we have two input files (unphased VCF+scaffold) and one output
                ).phased_rare_chunk
                for i in range(len(rare_regions))
            ]

            # 2B. Concatenate phased rare chunks
            concatenate_rare_chunks(
                b=batch,
                rare_variants_chunks_list=rare_chunks_phased,
                output_vcf_name=output_filename,
                chrom=f'chr{i}',
                out_dir=f'{output_path}/shapeit5',
                storage=round(vcf_size*0.1)
            )

    batch.run()
