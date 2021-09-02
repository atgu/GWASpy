import pandas as pd, numpy as np
from gwaspy.utils.get_file_size import bytes_to_gb
from gwaspy.phasing.get_filebase import get_vcf_filebase
import hailtop.batch as hb


def create_windows_bed(reference: str = 'GRCh38',
                       max_win_size_cm: float = 10.0,
                       overlap_size_cm: float = 2.0,
                       out_dir: str = None):
    print('creating BED file with overlapping windows to be used in splitting the input VCF')
    if reference == 'GRCh38':
        chrom_lens = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259,
                      'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 'chr9': 248956422, 'chr10': 242193529,
                      'chr11': 198295559, 'chr12': 190214555, 'chr13': 181538259, 'chr14': 170805979,
                      'chr15': 159345973, 'chr16': 145138636, 'chr17': 248956422, 'chr18': 242193529,
                      'chr19': 198295559, 'chr20': 190214555, 'chr21': 181538259, 'chr22': 170805979, 'chrX': 159345973}

        df_map = pd.read_csv('gs://african-seq-data/GWASpy/maps/genetic_map_hg38_withX.txt.gz', delim_whitespace=True,
                             compression='gzip', header=0, names=['CHR', 'POS', 'RATE', 'CM'])
    else:
        chrom_lens = {'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276, 'chr5': 180915260,
                      'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022, 'chr9': 141213431, 'chr10': 135534747,
                      'chr11': 135006516, 'chr12': 133851895, 'chr13': 115169878, 'chr14': 107349540,
                      'chr15': 102531392, 'chr16': 90354753, 'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983,
                      'chr20': 63025520, 'chr21': 48129895, 'chr22': 51304566, 'chrX': 155270560}

        df_map = pd.read_csv('gs://african-seq-data/GWASpy/maps/genetic_map_hg19_withX.txt.gz', delim_whitespace=True,
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

    df[['CHR', 'BEG', 'END']].to_csv(f'{out_dir}/GWASpy/Phasing/gwaspy.refscatter.bed', sep='\t', header=False,
                                     index=False)


def vcf_scatter(b: hb.batch.Batch,
                vcf_file: str = None,
                intervals_bed: str = None,
                docker_img: str = 'docker.io/lindonkambule/gwaspy:v1',
                memory: float = 26,
                out_dir: str = None):

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

    cmd = f'''
        set -euo pipefail
        mkdir vcfs
        echo *
        awk -F"\t" '{{print $1":"$2"-"$3"\t"NR-1}}' {bed} > regions.lines
        bcftools query --list-samples "{vcf}" | tee "{vcf_filename_no_ext}.sample_id.lines" | wc -l > n_smpls.int
        bcftools annotate \
            --no-version \
            --force \
            --output-type z \
            --remove ID,QUAL,INFO,^FMT/GT \
            --threads {threads} \
            "{vcf}" | \
        bcftools +scatter \
            --no-version \
            --output-type z \
            --output vcfs \
            --threads {threads} \
            --scatter-file regions.lines \
            --prefix "{vcf_filename_no_ext}."
        cut -f2 regions.lines | sed 's/^/vcfs\/{vcf_filename_no_ext}./;s/$/.vcf.gz/'
        echo vcfs/*
        '''

    scatter.command(cmd)

    scatter.command(f'mv vcfs {scatter.ofile}')
    b.write_output(scatter.ofile, f'{out_dir}/GWASpy/Phasing/{vcf_filename_no_ext}/scatter_vcfs')

