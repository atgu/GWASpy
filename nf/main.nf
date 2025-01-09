#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PHASE_ARRAY; PHASE_COMMON; LIGATE_COMMON; PHASE_RARE; CONCATENATE_CHUNKS as CONCATENATE_PHASED_CHUNKS } from './modules/phasing'
include { IMPUTE5 } from './modules/imputation'
include { CONCATENATE_CHUNKS as CONCATENATE_IMPUTED_CHUNKS } from './modules/phasing'


process CONVERT_REF {
    cpus 8
    memory { 4.GB * ta    container 'docker.io/lindonkambule/gwaspy_phase_impute:latest'
sk.attempt }
    tag "convert reference to XCF: chr${chrom}"

    input:
    tuple val(chrom), file(ref), file(ref_idx)

    output:
    tuple val(chrom), path("ref.chr${chrom}_xcf.bcf"), path("ref.chr${chrom}_xcf.bcf.csi"), path("ref.chr${chrom}_xcf.bin"), path("ref.chr${chrom}_xcf.fam")

    script:
    """
    xcftools_static view -i ${ref} -o ref.chr${chrom}_xcf.bcf -O sh -r chr${chrom} -T${task.cpus} -m 0.03125
    """
}


// This function os for extracting just the chromosome number from a filename
// Currently, it can only handle chromosomes 1-22. Working on handling strings that end with X
def getChromosome = { str ->
    def matcher = str =~ /\d{1,2}$/
    if (matcher.find()) {
        return matcher.group(0)
    } else {
        return null
    }
}

// ./nextflow run main.nf -c nextflow.config -profile gbatch -params-file params.json


workflow {
    // 1: INPUT FILES
    // 1.1: input VCF split by chromosome
    def input_chr_files = params.input.replaceFirst(/CNUMBER/, "*")
    Channel
        .fromFilePairs("${input_chr_files}.{bcf,bcf.csi,vcf,vcf.tbi}", size: 2)
        .map { chr, input_file -> tuple( getChromosome(chr), file(input_file[0]), file(input_file[1]) ) }
        .set { input_files }
    input_files.view()
    // 1.2: reference VCF split by chromosome
    def reference_chr_files = params.ref.replaceFirst(/CNUMBER/, "*")
    if ( params.ref_format == "vcf" ) {
        Channel
            .fromFilePairs("${reference_chr_files}.{bcf,bcf.csi}", size: 2)
            .map { chr, ref_file -> tuple( getChromosome(chr), file(ref_file[0]), file(ref_file[1]) ) }
            .filter { it[0] as Integer > 19 } // FOR TESTING USING FEWER CHROMOSOMES
            .set { ref_files_vcf }

        // ref_files_vcf.view()
        CONVERT_REF(ref_files_vcf)
        ref_converted_xcf = CONVERT_REF.out
        ref_converted_xcf
            .map {chr, bcf, csi, bin, fam -> tuple( chr, bcf, csi, bin, fam )}
            .set { ref_files }
    }
    else {
        Channel
            .fromFilePairs("${reference_chr_files}.{bcf,bcf.csi,bin,fam}", size: 4)
            .map { chr, ref_file -> tuple( getChromosome(chr), file(ref_file[0]), file(ref_file[1]), file(ref_file[2]), file(ref_file[3])) }
            .set { ref_files }
    }

    // ref_files.view()
    // 1.3: map file split by chromosome
    Channel
        .fromFilePairs("gs://hgdp-1kg/phasing/maps/b38/chr*.b38.gmap.gz", size: 1)
        .map { chr, map_file -> tuple( chr.replaceFirst(/chr/,""), file(map_file[0]) ) }
        .set { genetic_map_files }

    // 1.4: common variants chunks files split by chromosome
    def common_chunks_files = params.common_chunks.replaceFirst(/CNUMBER/, "*")
    Channel
        .fromPath("${common_chunks_files}", checkIfExists: true)
        .splitCsv(sep:'\t', header:['index', 'chrom', 'irg', 'org'])
        // ligate_static requires an ordered list of the phased common variants, hence we're keeping the index
        .map{ row -> tuple( row.chrom.replaceFirst(/chr/,""), row.irg, row.index) }
        .set { common_variants_chunks } // chrom, irg, index

    // 1.5: rare variants chunks files split by chromosome
    // srg = --scaffold-region
    // irg =  --input-region
    def rare_chunks_files = params.rare_chunks.replaceFirst(/CNUMBER/, "*")
    Channel
        .fromPath("${rare_chunks_files}", checkIfExists: true)
        .splitCsv(sep:'\t', header:['index', 'chrom', 'srg', 'irg', 'col5', 'col6', 'col7', 'col8'])
        .map { row -> tuple( row.chrom.replaceFirst(/chr/,""), row.srg, row.irg, row.index) }
        .set { rare_variants_chunks } // chrom, irg, org, index

    // 2: PHASING ARRAY OR SEQUENCING DATA
    if ( params.data_type == "array" ) {
        // 2.1: Phase SNP array data
        input_files // chr, input_vcf, input_idx
            .combine(ref_files, by: 0) // chr, ref_file, ref_idx, ref_bin, ref_fam
            .combine(genetic_map_files, by: 0) // chr, map_file
            .map { chrom, input_vcf, input_idx, ref_vcf, ref_idx, ref_bin, ref_fam, map_file ->
                    tuple( chrom, input_vcf, input_idx, ref_vcf, ref_idx, ref_bin, ref_fam, map_file ) }
            .set { phase_array_inputs }

        PHASE_ARRAY(phase_array_inputs, params.output_filename, "${params.output_path}/shapeit5/phase_array")

        // collect files for imputation
        phased_data = PHASE_ARRAY.out
        phased_data
            .map { chr, phased_vcf, phased_idx, phased_log -> tuple( chr, phased_vcf, phased_idx )}
            .set { fully_phased_data }

    } else {
        // 2.2 Phase WGS/WES data
        // 2.2.1: phase common chunks
        // combine the input, reference, map, and chunk files by chromosome+region
        input_files // chr, input_vcf, input_idx
            .combine(ref_files, by: 0) // chr, ref_file, ref_idx, ref_bin, ref_fam
            .combine(genetic_map_files, by: 0) // chr, map_file
            .combine(common_variants_chunks, by: 0) // chr, region, chunk_idx
            .map { chrom, input_vcf, input_index, ref_vcf, ref_idx, ref_bin, ref_fam, map_file, region, chunk_idx ->
                    tuple( chrom, input_vcf, input_index, ref_vcf, ref_idx, ref_bin, ref_fam, map_file, region, chunk_idx ) }
            .set { phase_common_inputs }
        PHASE_COMMON(phase_common_inputs, params.maf, "${params.output_path}/shapeit5/phase_common/logs")

        // 2.2.2 ligate phased common chunks
        // first make sure the phased files are grouped by chromosome and sorted by index
        phased_common_chunks_data = PHASE_COMMON.out
        phased_common_chunks_data
            .map {  chr, chunk_idx, vcf, vcf_idx, vcf_log -> [chr, [chunk_idx, vcf, vcf_idx] ] } //ligate_static requires indices
            .groupTuple(by: 0)
            .map { chr, data ->
                def sortedData = data.sort { it[0] as Integer } // leaving out as Integer will sort lexicographically instead of natural
                def sortedIndices = sortedData.collect { it[0]}
                def sortedFiles = sortedData.collect { it[1] }
                def sortedIdx = sortedData.collect { it[2] }

                return tuple(chr, sortedFiles, sortedIdx)
            }
            .set { ligate_common_inputs }

        LIGATE_COMMON(ligate_common_inputs, params.output_filename, "${params.output_path}/shapeit5/phase_common")
        scaffold_data = LIGATE_COMMON.out
        scaffold_data
            .map { chr, scaffold, scaffold_idx -> tuple(chr, scaffold, scaffold_idx)}
            .set { scaffold_data_inputs }

        // 2.2.3 phase rare chunks, using scaffolds from 2.2.2
        input_files // chr, input_vcf, input_idx
            .combine(genetic_map_files, by: 0) // chr, map_file
            .combine(rare_variants_chunks, by: 0) // chr, srg, irg, index
            .combine(scaffold_data_inputs, by: 0) // chr, scaffold, scaffold_idx
            .map { chrom, input_vcf, input_idx, map_file, srg, irg, chunk_idx, scaffold, scaffold_idx ->
                    tuple( chrom, input_vcf, input_idx, map_file, srg, irg, chunk_idx, scaffold, scaffold_idx ) }
            .set { phase_rare_inputs }

        // exit code 139 when running Docker desktop: https://www.nextflow.io/blog/2021/setup-nextflow-on-windows.html
        PHASE_RARE(phase_rare_inputs, "${params.output_path}/shapeit5/phase_rare/logs")

         // 2.2.4 concatenate phased rare chunks from 2.2.3
        phased_rare_chunks_data = PHASE_RARE.out
        phased_rare_chunks_data
            .map {  chr, chunk_idx, vcf, vcf_idx, vcf_log -> [chr, [chunk_idx, vcf, vcf_idx] ] } //ligate_static requires indices
            .groupTuple(by: 0)
            .map { chr, data ->
                def sortedData = data.sort { it[0] as Integer } // leaving out as Integer will sort lexicographically instead of natural
                def sortedIndices = sortedData.collect { it[0]}
                def sortedFiles = sortedData.collect { it[1] }
                def sortedIdx = sortedData.collect { it[2] }

                return tuple(chr, sortedFiles, sortedIdx, "${params.output_filename}_chr${chr}.full.shapeit5_rare.bcf")
            }
            .set { concatenate_rare_inputs }

        CONCATENATE_PHASED_CHUNKS(concatenate_rare_inputs, "${params.output_path}/shapeit5/phase_rare")

        // collect files for imputation
        phased_data = CONCATENATE_PHASED_CHUNKS.out
        phased_data
            .map { chr, phased_file, phased_idx -> tuple( chr, phased_file, phased_idx )}
            .set { fully_phased_data }
    }

    if ( params.impute ) {
        // 3: GENOTYPE IMPUTATION
        // 3.1 Impute genotypes in chunks
        fully_phased_data // chr, phased_chr_file, phased_chr_file_idx
            .combine(ref_files, by: 0) // chr, ref_file, ref_idx, ref_bin, ref_fam
            .combine(genetic_map_files, by: 0) // chr, map_file
            .combine(rare_variants_chunks, by: 0) // chr, srg, irg, index
            .map { chrom, input_vcf, input_idx, ref_vcf, ref_idx, ref_bin, ref_fam, map_file, srg, irg, chunk_idx ->
                    tuple( chrom, input_vcf, input_idx, ref_vcf, ref_idx, ref_bin, ref_fam, map_file, srg, irg, chunk_idx, "${params.output_path}/impute5/imputed_chunks/chr${chrom}") }
            .set { impute_data_inputs }

        IMPUTE5(impute_data_inputs)

         // 3.2 Concatenate imputed chunks from 4.1
        imputed_chunks_data = IMPUTE5.out
        imputed_chunks_data
            .map {  chr, chunk_idx, vcf, vcf_idx -> [chr, [chunk_idx, vcf, vcf_idx] ] } // ligate_static requires indices
            .groupTuple(by: 0)
            .map { chr, data ->
                def sortedData = data.sort { it[0] as Integer }
                def sortedIndices = sortedData.collect { it[0]}
                def sortedFiles = sortedData.collect { it[1] }
                def sortedIdx = sortedData.collect { it[2] }

                return tuple(chr, sortedFiles, sortedIdx, "${params.output_filename}_chr${chr}.full.shapeit5_phased.impute5_imputed.bcf")
            }
            .set { concatenate_imputed_inputs }

        CONCATENATE_IMPUTED_CHUNKS(concatenate_imputed_inputs, "${params.output_path}/impute5/imputed_merged")

    }
}
