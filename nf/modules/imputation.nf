#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process IMPUTE5 {
    cpus 8
    memory { 16.GB * task.attempt }
    container 'docker.io/lindonkambule/gwaspy_phase_impute:latest'
    tag "impute: ${irg}"
    publishDir "${out_directory}", overwrite: true, mode:'copy', pattern: '*.bcf*'

    input:
        tuple val(chrom), file(input), file(input_idx), file(ref), file(ref_idx), file(ref_bin), file(ref_fam), file(map_file), val(srg), val(irg), val(chk), val(out_directory)

    output:
    tuple val(chrom), val(chk), path("${chk}imputed.chr${chrom}.bcf"), path("${chk}imputed.chr${chrom}.bcf.csi")

    // IMPUTE5 automatically indexes output file
    script:
    """
    impute5_v1.2.0_static \
        --h ${ref} \
        --g ${input} \
        --m ${map_file} \
        --r ${irg} \
        --buffer-region ${srg} \
        --o ${chk}imputed.chr${chrom}.bcf
    """
}
