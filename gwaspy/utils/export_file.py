import hail as hl


def export_qced_file(mt: hl.MatrixTable, out_dir: str, basename: str, export_type='hail'):
    outname = basename + '_qced'

    if export_type == 'hail':
        mt.write('{}GWASpy/Preimp_QC/{}.mt'.format(out_dir, outname), overwrite=True)

    elif export_type == 'plink':
        hl.export_plink(dataset=mt, output='{}GWASpy/Preimp_QC/{}'.format(out_dir, outname), fam_id=mt.fam_id,
                        ind_id=mt.s, pat_id=mt.pat_id, mat_id=mt.mat_id, is_female=mt.is_female, pheno=mt.is_case,
                        varid=mt.rsid)

    else:
        hl.export_vcf(mt, '{}GWASpy/Preimp_QC/{}.vcf.bgz'.format(out_dir, outname))
