import hail as hl


def read_plink(dirname: str, basename: str) -> hl.MatrixTable:
    if not hl.hadoop_exists(dirname + basename + '.mt'):
        in_mt: hl.MatrixTable = hl.import_plink(bed=dirname + basename + '.bed',
                                                bim=dirname + basename + '.bim',
                                                fam=dirname + basename + '.fam',
                                                n_partitions=100)
        in_mt.write(dirname + basename + '.mt')
    return hl.read_matrix_table(dirname + basename + '.mt')


def read_vcf(dirname: str, vcf: str, annotations: str) -> hl.MatrixTable:
    hl.import_vcf(vcf).write('{}preimpQC.mt'.format(dirname), overwrite=True)
    in_mt = hl.read_matrix_table('{}preimpQC.mt'.format(dirname))
    ann = hl.import_table(annotations, impute=True).key_by('Sample')
    in_mt = in_mt.annotate_cols(annotations=ann[in_mt.s])
    # need to change reported sex to True/False, can update how this is done later, ideally don't want to hardcode
    # this will not work for unreported sex but will work for missing values
    in_mt = in_mt.annotate_cols(annotations=in_mt.annotate(Sex=hl.if_else((in_mt.annotations.Sex == 'F' |
                                                                           in_mt.annotations.Sex == 2 |
                                                                           in_mt.annotations.Sex == 'Female'),
                                                                          True, False)))
    # add a check to make sure file is formatted as we'd expect else quit and throw error
    return in_mt


def read_mt(dirname: str, basename: str) -> hl.MatrixTable:
    in_mt: hl.MatrixTable = hl.read_matrix_table(dirname + basename + ".mt")

    return in_mt


def read_infile(
        input_type: str = None,
        dirname: str = None, basename: str = None):

    global mt

    if input_type == 'plink':
        mt = read_plink(dirname, basename)

    elif input_type == 'vcf':
        print("VCF Support comming")
        # mt = read_vcf(dirname, vcf, annotations)

    else:
        mt = read_mt(dirname, basename)

    return mt
