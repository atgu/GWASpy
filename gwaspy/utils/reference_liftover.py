import hail as hl
from gwaspy.utils.read_file import read_infile


def liftover_to_grch38(
        input_type: str = None,
        dirname: str = None,
        basename: str = None):

    lifted_over = f'{dirname}{basename}.liftover.grch38.mt'
    print('\nLifting over to GRCh38')
    mt = read_infile(input_type=input_type, dirname=dirname, basename=basename)

    rg37 = hl.get_reference('GRCh37')
    rg38 = hl.get_reference('GRCh38')
    rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38)

    mt = mt.annotate_rows(new_locus=hl.liftover(mt.locus, 'GRCh38', include_strand=True), old_locus=mt.locus)
    mt = mt.filter_rows(hl.is_defined(mt.new_locus) & ~mt.new_locus.is_negative_strand)

    mt = mt.key_rows_by(locus=mt.new_locus.result, alleles=mt.alleles)

    print(f'\nWriting out data lifted-over to GRCh38 to: {lifted_over}')
    mt.write(lifted_over)

    return hl.read_matrix_table(lifted_over)
