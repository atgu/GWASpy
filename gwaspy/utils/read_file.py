import hail as hl
from gwaspy.utils.sample_annotations import add_sample_annotations


def read_plink(dirname: str, basename: str) -> hl.MatrixTable:

    in_mt: hl.MatrixTable = hl.import_plink(bed=dirname + basename + '.bed',
                                            bim=dirname + basename + '.bim',
                                            fam=dirname + basename + '.fam',
                                            block_size=16)

    return in_mt


def read_vcf(dirname: str, basename: str) -> hl.MatrixTable:
    hl._set_flags(no_whole_stage_codegen='1')
    vcf_file = '{}{}.vcf.gz'.format(dirname, basename)
    hl.import_vcf(vcf_file, force_bgz=True, block_size=16).write('{}GWASpy.preimpQC.mt'.format(dirname), overwrite=True)
    in_mt = hl.read_matrix_table('{}GWASpy.preimpQC.mt'.format(dirname))

    # Unlike array data, a VCF might have multi-allelic sites
    # split multi-allelic sites into bi-allelic
    print("Checking for multi-allelic sites")
    pre_filt_multi_n = in_mt.count_rows()
    bi = in_mt.filter_rows(hl.len(in_mt.alleles) == 2)
    bi = bi.annotate_rows(a_index=hl.missing(hl.tint))  # when we update Hail version, use hl.missing instead of hl.null
    bi = bi.annotate_rows(was_split=False)

    multi = in_mt.filter_rows(hl.len(in_mt.alleles) > 2)
    split = hl.split_multi_hts(multi)

    in_mt = split.union_rows(bi)
    pos_filt_multi_n = in_mt.count_rows()
    print("Number of multi-allelic SNPs in VCF file: {}".format(pos_filt_multi_n-pre_filt_multi_n))

    return in_mt


def read_mt(dirname: str, basename: str) -> hl.MatrixTable:
    print(dirname + basename + ".mt")
    in_mt: hl.MatrixTable = hl.read_matrix_table(dirname + basename + ".mt")

    return in_mt


def read_infile(
        input_type: str = None,
        dirname: str = None, basename: str = None,
        **kwargs):

    global mt

    # vcf = kwargs.get('vcf')
    annotations = kwargs.get('annotations')

    if input_type == 'plink':
        mt = read_plink(dirname, basename)

    elif input_type == 'vcf':
        mt = read_vcf(dirname, basename)

    else:
        mt = read_mt(dirname, basename)

    if annotations:
        mt = add_sample_annotations(mt, annotations)

    return mt
