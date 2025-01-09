#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process PHASE_ARRAY {
    cpus 8
    memory { 32.GB * task.attempt }
    container 'docker.io/lindonkambule/gwaspy_phase_impute:latest'
    tag "phase_array: chr${chrom}"
    publishDir "${out_directory}", overwrite: true, mode:'copy', pattern: '*.{bcf,bcf.csi,log}'

    input:
        tuple val(chrom), file(input), file(input_idx), file(ref), file(ref_idx), file(ref_bin), file(ref_fam), file(map_file)
        val(output_filename)
        val(out_directory)

    output:
        tuple val(chrom), path("${output_filename}_chr${chrom}_b0_v2.b38.sorted.phased.bcf"), path("${output_filename}_chr${chrom}_b0_v2.b38.sorted.phased.bcf.csi"), path("${output_filename}_chr${chrom}_b0_v2.b38.sorted.phased.log")

    script:
    """
    phase_common_static \
        --input ${input} \
        --reference ${ref} \
        --map ${map_file} \
        --output ${output_filename}_chr${chrom}_b0_v2.b38.sorted.phased.bcf \
        --thread ${task.cpus} \
        --log ${output_filename}_chr${chrom}_b0_v2.b38.sorted.phased.log \
        --region chr${chrom}
    """
}


process PHASE_COMMON {
    cpus 8
    memory { 16.GB * task.attempt }
    container 'docker.io/lindonkambule/gwaspy_phase_impute:latest'
    tag "phase_common: ${region}"
    publishDir "${out_directory}", overwrite: true, mode:'copy', pattern: '*.log'

    input:
        tuple val(chrom), file(input), file(input_idx), file(ref), file(ref_idx), file(ref_bin), file(ref_fam), file(map_file), val(region), val(chk)
        val(maf)
        val(out_directory)

    output:
        tuple val(chrom), val(chk), path("${chk}phased_common.chr${chrom}.bcf"), path("${chk}phased_common.chr${chrom}.bcf.csi"), path("${chk}phased_common.chr${chrom}.log")

    script:
    """
    phase_common_static \
        --input ${input} \
        --reference ${ref} \
        --map ${map_file} \
        --output ${chk}phased_common.chr${chrom}.bcf \
        --thread ${task.cpus} \
        --log ${chk}phased_common.chr${chrom}.log \
        --filter-maf ${maf} \
        --region ${region}
    """
}


process LIGATE_COMMON {
    cpus 4
    memory { 4.GB * task.attempt }
    container 'docker.io/lindonkambule/gwaspy_phase_impute:latest'
    tag "ligate_common: chr${chrom}"
    publishDir "${out_directory}", overwrite: true, mode:'copy', pattern: '*shapeit5_common.bcf*'

    input:
        tuple val(chrom), path(inputs), path(inputs_index_files)
        val(output_filename)
        val(out_directory)

    output:
        tuple val(chrom), path("${output_filename}_chr${chrom}.shapeit5_common.bcf"), path("${output_filename}_chr${chrom}.shapeit5_common.bcf.csi")

    script:
    // chunks are ordered but in one line, add | tr " " "\n" so common_chunks_list_ligate.txt has one file per line
    """
    echo ${inputs} | tr " " "\n" > common_chunks_list_ligate.txt
    ligate_static \
        --input common_chunks_list_ligate.txt \
        --output ${output_filename}_chr${chrom}.shapeit5_common.bcf \
        --thread ${task.cpus} \
        --log ${output_filename}_chr${chrom}.shapeit5_common.log \
        --index
    """
}


process PHASE_RARE {
    cpus 4
    memory { 16.GB * task.attempt }
    container 'docker.io/lindonkambule/gwaspy_phase_impute:latest'
    tag "phase_rare: ${srg}"
    publishDir "${out_directory}", overwrite: true, mode:'copy', pattern: '*.log'

    input:
        // tuple( chrom, input_vcf, input_idx, map_file, srg, irg, idx, scaffold, scaffold_idx )
        tuple val(chrom), file(input), file(input_idx), file(map_file), val(srg), val(irg), val(chk), file(scaffold), file(scaffold_idx)
        val(out_directory)

    output:
        tuple val(chrom), val(chk), path("${chk}phased_rare.chr${chrom}.bcf"), path("${chk}phased_rare.chr${chrom}.bcf.csi"), path("${chk}phased_rare.chr${chrom}.log")

    script:
    // phase_rare_static does not have --reference option (because we're using a scaffold?)
    """
    phase_rare_static \
        --input ${input} \
        --input-region ${irg} \
        --scaffold ${scaffold} \
        --scaffold-region ${srg} \
        --map ${map_file} \
        --output ${chk}phased_rare.chr${chrom}.bcf \
        --thread ${task.cpus} \
        --log ${chk}phased_rare.chr${chrom}.log
    """
}


// can be used to concatenate chunks after phasing rare variants or after imputing with IMPUTE5
process CONCATENATE_CHUNKS {
    cpus 4
    memory { 8.GB * task.attempt }
    container 'docker.io/lindonkambule/gwaspy_phase_impute:latest'
    tag "concatenate_chunks: chr${chrom}"
    publishDir "${out_directory}", overwrite: true, mode:'copy', pattern: '*.bcf*'

    input:
        tuple val(chrom), path(inputs), path(inputs_idx), val(out_filename)
        val(out_directory)

    output:
        tuple val(chrom), path("${out_filename}"), path("${out_filename}.csi")

    script:
    // chunks are ordered but in one line, add | tr " " "\n" so common_chunks_list_ligate.txt has one file per line
    """
    echo ${inputs} | tr " " "\n" > list_concatenate.txt
    bcftools concat -n -f list_concatenate.txt -o ${out_filename}
    bcftools index ${out_filename}
    """
}
