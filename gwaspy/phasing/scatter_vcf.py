__author__ = 'Lindo Nkambule'

import pandas as pd, numpy as np
from gwaspy.utils.get_file_size import bytes_to_gb
from gwaspy.phasing.get_filebase import get_vcf_filebase
import hailtop.batch as hb
import pathlib
import ntpath
from typing import Union


def create_windows_bed(reference: str = 'GRCh38',
                       max_win_size_cm: float = 10.0,
                       overlap_size_cm: float = 2.0,
                       vcf_filebase: str = None,
                       out_dir: str = None):
    print('creating BED file with overlapping windows to be used in splitting the input VCF')
    maps_path = 'https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables'
    if reference == 'GRCh38':
        chrom_lens = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259,
                      'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 'chr9': 248956422, 'chr10': 242193529,
                      'chr11': 198295559, 'chr12': 190214555, 'chr13': 181538259, 'chr14': 170805979,
                      'chr15': 159345973, 'chr16': 145138636, 'chr17': 248956422, 'chr18': 242193529,
                      'chr19': 198295559, 'chr20': 190214555, 'chr21': 181538259, 'chr22': 170805979, 'chrX': 159345973}

        df_map = pd.read_csv(f'{maps_path}/genetic_map_hg38_withX.txt.gz', delim_whitespace=True,
                             compression='gzip', header=0, names=['CHR', 'POS', 'RATE', 'CM'])
    else:
        chrom_lens = {'1': 249250621, '2': 243199373, '3': 198022430, '4': 191154276, '5': 180915260, '6': 171115067,
                      '7': 159138663, '8': 146364022, '9': 141213431, '10': 135534747, '11': 135006516, '12': 133851895,
                      '13': 115169878, '14': 107349540, '15': 102531392, '16': 90354753, '17': 81195210, '18': 78077248,
                      '19': 59128983, '20': 63025520, '21': 48129895, '22': 51304566, '23': 155270560}

        df_map = pd.read_csv(f'{maps_path}/genetic_map_hg19_withX.txt.gz', delim_whitespace=True,
                             compression='gzip', header=0, names=['CHR', 'POS', 'RATE', 'CM'])
    df_out = {}

    for chrom, df_group in df_map.groupby('CHR'):
        fai_chr = str(chrom) if str(chrom) in chrom_lens else 'chr' + str(chrom) if 'chr' + str(
            chrom) in chrom_lens else 'X' if 'X' in chrom_lens else 'chrX' if 'chrX' in chrom_lens else None
        if fai_chr:
            chr_cm_len = max(df_group['CM'])
            n_win = np.ceil((chr_cm_len - overlap_size_cm) / (max_win_size_cm - overlap_size_cm))
            win_size = (chr_cm_len - overlap_size_cm) / n_win + overlap_size_cm
            cm_begs = (win_size - overlap_size_cm) * np.arange(1, n_win)
            cm_ends = (win_size - overlap_size_cm) * np.arange(1, n_win) + overlap_size_cm
            pos_begs = np.concatenate(
                ([1], np.interp(cm_begs, df_group['CM'], df_group['POS'], period=np.inf).astype(int)))
            pos_ends = np.concatenate(
                (np.interp(cm_ends, df_group['CM'], df_group['POS'], period=np.inf).astype(int), [chrom_lens[fai_chr]]))
            df_out[fai_chr] = pd.DataFrame.from_dict({'CHR': fai_chr, 'BEG': pos_begs, 'END': pos_ends})

    df = pd.concat([df_out[fai_chr] for fai_chr in chrom_lens.keys()])

    df[['CHR', 'BEG', 'END']].to_csv(f'{out_dir}/GWASpy/{vcf_filebase}/Phasing/refscatter.bed', sep='\t', header=False,
                                     index=False)


def vcf_scatter(b: hb.batch.Batch,
                vcf_file: str = None,
                intervals_bed: str = None,
                docker_img: str = 'docker.io/lindonkambule/gwaspy:v1',
                memory: float = 26,
                out_dir: str = None):

    global cmd

    vcf_size = bytes_to_gb(vcf_file)
    disk_size = round(10.0 + 3.0 * vcf_size)
    cpu = 2 * round(memory / 13) if memory > 6.5 else 1
    threads = cpu - 1

    vcf_filename_no_ext = get_vcf_filebase(vcf_file)

    vcf = b.read_input(vcf_file)
    bed = b.read_input(intervals_bed)

    scatter = b.new_job(name=f'scatter-{vcf_filename_no_ext}')
    scatter.cpu(cpu)
    scatter.memory(f'{memory}Gi')
    scatter.storage(f'{disk_size}Gi')
    scatter.image(docker_img)

    # work out if file is BCF or VCF
    vcf_basename = ntpath.basename(vcf_file)
    if '.vcf' in pathlib.Path(vcf_basename).suffixes:
        cmd = f'''
            set -euo pipefail
            mkdir vcfs
            awk -F"\t" '{{print $1":"$2"-"$3"\t"NR-1}}' {bed} > regions.lines
            bcftools query --list-samples "{vcf}" | tee "{vcf_filename_no_ext}.sample_id.lines" | wc -l > n_smpls.int
            bcftools view --threads {threads} --output-type u "{vcf}" | \
            bcftools annotate \
                --no-version \
                --force \
                --output-type u \
                --remove ID,QUAL,INFO,^FMT/GT \
                --threads {threads} | \
            bcftools +scatter \
                --no-version \
                --output-type b \
                --output vcfs \
                --threads {threads} \
                --scatter-file regions.lines \
                --prefix "{vcf_filename_no_ext}."
        '''
    elif '.bcf' in pathlib.Path(vcf_basename).suffixes:
        cmd = f'''
            set -euo pipefail
            mkdir vcfs
            awk -F"\t" '{{print $1":"$2"-"$3"\t"NR-1}}' {bed} > regions.lines
            bcftools query --list-samples "{vcf}" | tee "{vcf_filename_no_ext}.sample_id.lines" | wc -l > n_smpls.int
            bcftools annotate \
                --no-version \
                --force \
                --output-type u \
                --remove ID,QUAL,INFO,^FMT/GT \
                --threads {threads} \
                "{vcf}" | \
            bcftools +scatter \
                --no-version \
                --output-type b \
                --output vcfs \
                --threads {threads} \
                --scatter-file regions.lines \
                --prefix "{vcf_filename_no_ext}."
        '''

    else:
        raise SystemExit('Unsupported file format. Choices are VCF OR BCF')

    scatter.command(cmd)

    scatter.command(f'mv vcfs {scatter.vcfs}')
    scatter.command(f'mv regions.lines {scatter.regions}')
    b.write_output(scatter.vcfs, f'{out_dir}/GWASpy/{vcf_filename_no_ext}/Phasing/scatter_vcfs')
    b.write_output(scatter.regions, f'{out_dir}/GWASpy/{vcf_filename_no_ext}/Phasing/regions.lines')


def run_scatter(backend: Union[hb.ServiceBackend, hb.LocalBackend] = None,
                input_vcf: str = None,
                reference: str = 'GRCh38',
                max_win_size_cm: float = 10.0,
                overlap_size_cm: float = 2.0,
                scatter_memory: int = 26,
                out_dir: str = None):

    print(f'\n1. SCATTER {input_vcf}\n')
    vcf_filebase = get_vcf_filebase(input_vcf)

    create_windows_bed(reference=reference, max_win_size_cm=max_win_size_cm, out_dir=out_dir, vcf_filebase=vcf_filebase,
                       overlap_size_cm=overlap_size_cm)

    scatter_b = hb.Batch(backend=backend, name=f'scatter-{vcf_filebase}')

    vcf_scatter(
        b=scatter_b,
        vcf_file=input_vcf,
        intervals_bed=f'{out_dir}/GWASpy/{vcf_filebase}/Phasing/refscatter.bed',
        memory=scatter_memory,
        out_dir=out_dir
    )

    scatter_b.run()


