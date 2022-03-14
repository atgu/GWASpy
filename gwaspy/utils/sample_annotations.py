__author__ = 'Lindo Nkambule'

import hail as hl
import sys


def add_sample_annotations(mt: hl.MatrixTable, annotations: str) -> hl.MatrixTable:
    # use annotations file to annotate VCF
    ann = hl.import_table(annotations, impute=False,
                          types={'Sample': hl.tstr, 'Sex': hl.tstr, 'Pheno': hl.tstr}).key_by('Sample')
    ann_cols = dict(ann.row)

    mt = mt.annotate_cols(annotations=ann[mt.s])

    if 'is_female' not in mt.col:
        if 'Sex' in ann_cols:
            mt = mt.annotate_cols(is_female=hl.if_else(((mt.annotations.Sex == 'F') |
                                                        (mt.annotations.Sex == str(2)) |
                                                        (mt.annotations.Sex == 'True') |
                                                        (mt.annotations.Sex == 'Female')),
                                                       True, False))
        else:
            print('Sex column is missing from annotations file. Please add it and run GWASpy again')
            sys.exit(2)

    if 'is_case' not in mt.col:
        if 'Pheno' in ann_cols:
            mt = mt.annotate_cols(is_case=hl.if_else(((mt.annotations.Pheno == str(2)) |
                                                      (mt.annotations.Pheno == 'True') |
                                                      (mt.annotations.Pheno == 'Case')),
                                                     True, False))

    return mt

