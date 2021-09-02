__author__ = 'Lindo Nkambule'

import ntpath


def get_vcf_filebase(file: str = None):

    vcf_name = ntpath.basename(file)
    if vcf_name.endswith('.gz'):
        file_no_ext = vcf_name[:-7]
    elif vcf_name.endswith('.bgz'):
        file_no_ext = vcf_name[:-8]
    else:
        file_no_ext = vcf_name[:-4]

    return file_no_ext
