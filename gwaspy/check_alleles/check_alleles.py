_author__ = 'Lindo Nkambule'

import hailtop.batch as hb
import hailtop.fs as hfs

from hailtop.batch.job import Job


def size(file: str):
    """
    Convert the size from bytes to GiB
    :param file: path to file, str
    :return: file size in GiB
    """
    file_info = hfs.stat(file)   # returns a named tuple
    size_gigs = file_info.size / (1024 * 1024 * 1024)

    return size_gigs


def check_alleles_workflow(
        batch: hb.Batch = None,
        input_path: str = None,
        reference_path: str = None,
        output_filename: str = None,
        step: str = "check",
        fix_mode: str = "top",
        output_path: str = None):

    def get_stats(
            b: hb.batch.Batch,
            job_name: str = None,
            vcf: hb.ResourceGroup = None,
            ref_fasta: hb.ResourceGroup = None,
            output_name: str = None,
            out_dir: str = None,
            ncpu: int = 8,
            memory: str = 'standard',
            storage: int = None,
            img: str = 'docker.io/lindonkambule/gwaspy_phase_impute:latest',
    ) -> Job:
        j = b.new_job(name=f'Check alleles: {job_name}')

        j.image(img)
        j.memory(memory)
        j.cpu(ncpu)
        j.storage(f'{storage}Gi')

        j.command(
            f"""
            bcftools +fixref {vcf['vcf']} -- -f {ref_fasta['ref_fasta']} > stats.txt
            mv stats.txt {j.ofile}
            """
        )

        b.write_output(j.ofile,
                       f'{out_dir}/check_alleles/{output_name}.stats.txt')

        return j

    def fix_alleles(
            b: hb.batch.Batch,
            job_name: str = None,
            vcf: hb.ResourceGroup = None,
            ref_fasta: hb.ResourceGroup = None,
            allele_mode: str = "top",
            output_name: str = None,
            out_dir: str = None,
            ncpu: int = 8,
            memory: str = 'standard',
            storage: int = None,
            img: str = 'docker.io/lindonkambule/gwaspy_phase_impute:latest',
    ) -> Job:
        j = b.new_job(name=f'Fix alleles: {job_name}')

        j.image(img)
        j.memory(memory)
        j.cpu(ncpu)
        j.storage(f'{storage}Gi')

        j.declare_resource_group(
            fixed_file={
                'bcf': '{root}.bcf',
                'bcf.csi': '{root}.bcf.csi'
            }
        )

        j.command(
            f"""
            bcftools +fixref {vcf['vcf']} -Ob -o {j.fixed_file['bcf']} -- -f {ref_fasta['ref_fasta']} -m {allele_mode}
            """
        )

        b.write_output(j.stats,
                       f'{out_dir}/check_alleles/{output_name}.alleles.fixed')

        return j

    ref_fasta_in = batch.read_input_group(**{'ref_fasta': reference_path,
                                             'ref_fasta_index': f'{reference_path}.fai'})
    ref_size = round(size(reference_path))

    if "CNUMBER" in input_path:  # input VCF is already split by chromosome
        for i in range(1, 23):
            vcf_path = input_path.replace('CNUMBER', str(i))
            input_idx = f'{vcf_path}.tbi' if hfs.exists(f'{vcf_path}.tbi') else f'{vcf_path}.csi'

            if not hfs.exists(input_idx):
                raise SystemExit('Input file needs to be indexed (.tbi or .csi). Found none, exiting')

            chrom_vcf = batch.read_input_group(**{'vcf': vcf_path,
                                                  'index': input_idx})
            vcf_size = round(size(vcf_path))
            disk_size = int(round(5.0 + vcf_size + ref_size))

            if step == "check":
                get_stats(
                    b=batch,
                    job_name=vcf_path,
                    vcf=chrom_vcf,
                    ref_fasta=ref_fasta_in,
                    output_name=output_filename,
                    out_dir=output_path,
                    storage=disk_size
                )
            else:
                fix_alleles(
                    b=batch,
                    job_name=vcf_path,
                    vcf=chrom_vcf,
                    ref_fasta=ref_fasta_in,
                    allele_mode=fix_mode,
                    output_name=output_filename,
                    out_dir=output_path,
                    storage=disk_size
                )

    else: # one input file with all the chromosomes
        vcf_path = input_path
        input_idx = f'{vcf_path}.tbi' if hfs.exists(f'{vcf_path}.tbi') else f'{vcf_path}.csi'

        if not hfs.exists(input_idx):
            raise SystemExit('Input file needs to be indexed (.tbi or .csi). Found none, exiting')

        chrom_vcf = batch.read_input_group(**{'vcf': input_path,
                                              'index': input_idx})

        vcf_size = round(size(vcf_path))
        disk_size = int(round(5.0 + vcf_size + ref_size))

        if step == "check":
            get_stats(
                b=batch,
                job_name=vcf_path,
                vcf=chrom_vcf,
                ref_fasta=ref_fasta_in,
                output_name=output_filename,
                out_dir=output_path,
                storage=disk_size
            )
        else:
            fix_alleles(
                b=batch,
                job_name=vcf_path,
                vcf=chrom_vcf,
                ref_fasta=ref_fasta_in,
                allele_mode=fix_mode,
                output_name=output_filename,
                out_dir=output_path,
                storage=disk_size
            )

    batch.run()
