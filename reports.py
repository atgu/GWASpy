__author__ = 'Lindo Nkambule'

from pylatex import Document, Section, Subsection, Command, Center, Tabular
from pylatex.utils import NoEscape, bold


class MyDocument(Document):
    def __init__(self, basename, pre_qc_conts, post_qc_conts, count_results):
        super().__init__()

        self.preamble.append(Command('title', 'QC Report of {}'.format(basename)))
        self.preamble.append(Command('author', 'Anonymous author'))
        self.preamble.append(Command('date', NoEscape(r'\today')))
        self.append(NoEscape(r'\maketitle'))

        self._pre_counts = pre_qc_conts
        self._post_counts = post_qc_conts
        self._conts = count_results

    def general_info(self):
        pre_pheno_counts = str(self._pre_counts['is_case_counts']['case']) + ',' + \
                           str(self._pre_counts['is_case_counts']['control']) + ',' + \
                           str(self._pre_counts['is_case_counts']['unknown'])

        pre_sex_counts = str(self._pre_counts['is_female_counts']['male']) + ',' + \
                         str(self._pre_counts['is_female_counts']['female']) + ',' + \
                         str(self._pre_counts['is_female_counts']['unknown'])

        post_pheno_counts = str(self._post_counts['is_case_counts']['case']) + ',' + \
                            str(self._post_counts['is_case_counts']['control']) + ',' + \
                            str(self._post_counts['is_case_counts']['unknown'])

        post_sex_counts = str(self._post_counts['is_female_counts']['male']) + ',' + \
                          str(self._post_counts['is_female_counts']['female']) + ',' + \
                          str(self._post_counts['is_female_counts']['unknown'])

        pheno_diffs = str(self._pre_counts['is_case_counts']['case']-self._post_counts['is_case_counts']['case']) + ','\
                      + str(self._pre_counts['is_case_counts']['control']-self._post_counts['is_case_counts']['control']) +\
                      ',' + str(self._pre_counts['is_case_counts']['unknown']-self._post_counts['is_case_counts']['unknown'])
        sex_diffs = str(self._pre_counts['is_female_counts']['male']-self._post_counts['is_female_counts']['male']) + ','\
                    + str(self._pre_counts['is_female_counts']['female']-self._post_counts['is_female_counts']['female']) +\
                    ',' + str(self._pre_counts['is_female_counts']['unknown']-self._post_counts['is_female_counts']['unknown'])

        with self.create(Section('Flag')):
            self.append('A flags table should be here (STILL WORKING ON IT)')

        with self.create(Section('General Info')):
            with self.create(Subsection('Size of sample')):
                with self.create(Center()) as centered:
                    with centered.create(Tabular('|c|c|c|c|')) as table:
                        table.add_hline()
                        table.add_row((bold('Test'), bold('pre QC'), bold('post QC'), bold('exlcusion-N')))
                        table.add_hline()
                        table.add_row(('Cases,Controls,Missing', pre_pheno_counts,
                                       post_pheno_counts, pheno_diffs))
                        table.add_hline()
                        table.add_row(('Males,Females,Unspec', pre_sex_counts,
                                       post_sex_counts, sex_diffs))
                        table.add_hline()
                        table.add_row(('SNPs', self._pre_counts['n_variants'],
                                       self._post_counts['n_variants'],
                                       self._pre_counts['n_variants'] - self._post_counts['n_variants']))
                        table.add_hline()

            with self.create(Subsection('Exclusion overview')):
                with self.create(Center()) as centered:
                    with centered.create(Tabular('|l|l|')) as table:
                        table.add_hline()
                        table.add_row((bold('Filter'), bold('N')))
                        table.add_hline()
                        table.add_row(('SNPs: call rate < 0.950 (pre - filter)', self._conts['pre_geno'][True]))
                        table.add_row(('IDs: call rate (cases/controls) < 0.980', self._conts['mind'][True]))
                        table.add_row(('IDs: FHET outside +- 0.20 (cases/controls)', self._conts['fstat'][True]))
                        table.add_row(('IDs: Sex violations -excluded- (N-tested)', self._conts['sex_violations'][True]))
                        table.add_row(
                            ('IIDs: Sex warnings (undefined phenotype / ambiguous genotypes)',
                             0))
                        table.add_row(('SNPs: call rate < 0.980', self._conts['geno'][True]))
                        table.add_row(('SNPs: missing diference > 0.020', self._conts['cr_diff'][True]))
                        table.add_row(('SNPs: without valid association p-value (invariant)',
                                       self._conts['monomorphic_var'][True]))
                        table.add_row(('SNPs: HWE-controls < -6', self._conts['hwe_con'][True]))
                        table.add_row(('SNPs: HWE-cases < -10', self._conts['hwe_cas'][True]))
                        table.add_hline()


