import hail as hl


def read_plink(dirname: str, basename: str) -> hl.MatrixTable:
    if not hl.hadoop_exists(dirname + basename + '.mt'):
        mt = hl.import_plink(bed=dirname + basename + '.bed',
                             bim=dirname + basename + '.bim',
                             fam=dirname + basename + '.fam',
                             n_partitions=24)
        mt.write(dirname + basename + '.mt')
    return hl.read_matrix_table(dirname + basename + '.mt')


def read_vcf(dirname: str, vcf: str, annotations: str) -> hl.MatrixTable:
    hl.import_vcf(vcf).write('{}preimpQC.mt'.format(dirname), overwrite=True)
    mt = hl.read_matrix_table('{}preimpQC.mt'.format(dirname))
    ann = hl.import_table(annotations, impute=True).key_by('Sample')
    mt = mt.annotate_cols(annotations=ann[mt.s])
    # need to change reported sex to True/False, can update how this is done later, ideally don't want to hardcode
    # this will not work for unreported sex but will work for missing values
    mt = mt.annotate_cols(annotations=mt.annotate(Sex=hl.if_else((mt.annotations.Sex == 'F' | mt.annotations.Sex == 2 |
                                                                  mt.annotations.Sex == 'Female'), True, False)))
    # add a check to make sure file is formatted as we'd expect else quit and throw error
    return mt


def read_mt(dirname: str, basename: str) -> hl.MatrixTable:
    mt: hl.MatrixTable = hl.read_matrix_table(dirname + basename + ".mt")

    return mt
