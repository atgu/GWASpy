import hail as hl


def read_plink(dirname: str, basename: str) -> hl.MatrixTable:
    global in_mt
    if not hl.hadoop_exists(dirname + basename + '.mt'):
        in_mt: hl.MatrixTable = hl.import_plink(bed=dirname + basename + '.bed',
                                                bim=dirname + basename + '.bim',
                                                fam=dirname + basename + '.fam',
                                                block_size=16)
    return in_mt


def read_vcf(dirname: str, basename: str, annotations: str) -> hl.MatrixTable:
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

    # use annotations file to annotate VCF
    ann = hl.import_table(annotations, impute=False,
                          types={'Sample': hl.tstr, 'Sex': hl.tstr, 'Pheno': hl.tstr}).key_by('Sample')
    in_mt = in_mt.annotate_cols(annotations=ann[in_mt.s])
    # this will not work for unreported sex but will work for missing values
    in_mt = in_mt.annotate_cols(is_female=hl.if_else(((in_mt.annotations.Sex == 'F') |
                                                      (in_mt.annotations.Sex == str(2)) |
                                                      (in_mt.annotations.Sex == 'True') |
                                                      (in_mt.annotations.Sex == 'Female')),
                                                     True, False))

    # table_cols = dict(ann.row)  # get column fields
    if 'Pheno' in ann.row:
        in_mt = in_mt.annotate_cols(is_case=hl.if_else(((in_mt.annotations.Pheno == str(2)) |
                                                        (in_mt.annotations.Pheno == 'True') |
                                                        (in_mt.annotations.Pheno == 'Case')),
                                                       True, False))
    # add a check to make sure file is formatted as we'd expect else quit and throw error
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
        mt = read_vcf(dirname, basename, annotations)

    else:
        mt = read_mt(dirname, basename)

    return mt
