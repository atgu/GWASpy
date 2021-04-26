import hail as hl


def export_qced_file(mt: hl.MatrixTable, out_dir: str, basename: str, export_type='hail'):
    outname = basename + '_qced'

    if export_type == 'hail':
        mt.write('{}{}.mt'.format(out_dir, outname), overwrite=True)

    elif export_type == 'plink':
        hl.export_plink(mt, '{}{}'.format(out_dir, outname))

    else:
        hl.export_vcf(mt, '{}{}.vcf.bgz'.format(out_dir, outname))
