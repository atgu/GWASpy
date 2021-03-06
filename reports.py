__author__ = 'Lindo Nkambule'

from pylatex import Document, Section, Subsection, Command, Center, Tabular, NewPage, Figure, SubFigure
from pylatex.utils import NoEscape, bold


class MyDocument(Document):
    def __init__(self, basename):
        super().__init__()

        self.preamble.append(Command('title', 'QC Report of {}'.format(basename)))
        self.preamble.append(Command('author', 'Anonymous author'))
        self.preamble.append(Command('date', NoEscape(r'\today')))
        self.append(NoEscape(r'\maketitle'))

    def general_info(self, pre_qc_conts, post_qc_conts, count_results, pre_filter, id_cr, fhet_thresh, var_cr,
                     miss_diff, hwe_con, hwe_cas):
        pre_pheno_counts = str(pre_qc_conts['is_case_counts']['case']) + ', ' + \
                           str(pre_qc_conts['is_case_counts']['control']) + ', ' + \
                           str(pre_qc_conts['is_case_counts']['unknown'])

        pre_sex_counts = str(pre_qc_conts['is_female_counts']['male']) + ', ' + \
                         str(pre_qc_conts['is_female_counts']['female']) + ', ' + \
                         str(pre_qc_conts['is_female_counts']['unknown'])

        post_pheno_counts = str(post_qc_conts['is_case_counts']['case']) + ', ' + \
                            str(post_qc_conts['is_case_counts']['control']) + ', ' + \
                            str(post_qc_conts['is_case_counts']['unknown'])

        post_sex_counts = str(post_qc_conts['is_female_counts']['male']) + ', ' + \
                          str(post_qc_conts['is_female_counts']['female']) + ', ' + \
                          str(post_qc_conts['is_female_counts']['unknown'])

        pheno_diffs = str(pre_qc_conts['is_case_counts']['case']-post_qc_conts['is_case_counts']['case']) + ', '\
                      + str(pre_qc_conts['is_case_counts']['control']-post_qc_conts['is_case_counts']['control']) +\
                      ', ' + str(pre_qc_conts['is_case_counts']['unknown']-post_qc_conts['is_case_counts']['unknown'])
        sex_diffs = str(pre_qc_conts['is_female_counts']['male']-post_qc_conts['is_female_counts']['male']) + ', '\
                    + str(pre_qc_conts['is_female_counts']['female']-post_qc_conts['is_female_counts']['female']) +\
                    ', ' + str(pre_qc_conts['is_female_counts']['unknown']-post_qc_conts['is_female_counts']['unknown'])

        with self.create(Section('Flag')):
            self.append('A flags table should be here (STILL WORKING ON IT)')

        with self.create(Section('General Info')):
            with self.create(Subsection('Size of sample')):
                with self.create(Center()) as centered:
                    with centered.create(Tabular('|c|c|c|c|')) as table:
                        table.add_hline()
                        table.add_row((bold('Test'), bold('pre QC'), bold('post QC'), bold('exlcusion-N')))
                        table.add_hline()
                        table.add_row(('Cases, Controls, Missing', pre_pheno_counts,
                                       post_pheno_counts, pheno_diffs))
                        table.add_hline()
                        table.add_row(('Males, Females, Unspec', pre_sex_counts,
                                       post_sex_counts, sex_diffs))
                        table.add_hline()
                        table.add_row(('SNPs', pre_qc_conts['n_variants'],
                                       post_qc_conts['n_variants'],
                                       pre_qc_conts['n_variants'] - post_qc_conts['n_variants']))
                        table.add_hline()

            with self.create(Subsection('Exclusion overview')):
                with self.create(Center()) as centered:
                    with centered.create(Tabular('|l|l|')) as table:
                        table.add_hline()
                        table.add_row((bold('Filter'), bold('N')))
                        table.add_hline()
                        table.add_row(('SNPs: call rate < {} (pre - filter)'.format(pre_filter),
                                       count_results['pre_geno'][True]))
                        table.add_row(('IDs: call rate (cases/controls) < {}'.format(id_cr),
                                       count_results['mind'][True]))
                        table.add_row(('IDs: FHET outside +- {} (cases/controls)'.format(fhet_thresh),
                                       count_results['fstat'][True]))
                        table.add_row(('IDs: Sex violations -excluded- (N-tested)', count_results['sex_violations'][True]))
                        table.add_row(
                            ('IDs: Sex warnings (undefined phenotype / ambiguous genotypes)',
                             count_results['sex_warnings'][True]))
                        table.add_row(('SNPs: call rate < {}'.format(var_cr), count_results['geno'][True]))
                        table.add_row(('SNPs: missing diference > {}'.format(miss_diff), count_results['cr_diff'][True]))
                        table.add_row(('SNPs: without valid association p-value (invariant)',
                                       count_results['monomorphic_var'][True]))
                        table.add_row(('SNPs: HWE-controls < {}'.format(hwe_con), count_results['hwe_con'][True]))
                        table.add_row(('SNPs: HWE-cases < {}'.format(hwe_cas), count_results['hwe_cas'][True]))
                        table.add_hline()

    def manhattan_sec(self, qq_pre_path, qq_pos_path, man_pre_path, man_pos_path, table_results):
        self.append(NewPage())

        with self.create(Section('Manhattan')):
            with self.create(Subsection('Basic stats')):
                with self.create(Center()) as centered:
                    with centered.create(Tabular('|c|c|c|')) as table:
                        table.add_hline()
                        table.add_row((bold('Description'), bold('Pre-QC'), bold('Post-QC')))
                        table.add_hline()
                        table.add_row(('Number of GWAS hits', table_results[0], table_results[1]))
                        table.add_hline()
                        table.add_row(('Lambda GC', table_results[2], table_results[3]))
                        table.add_hline()
                        # calculate lambda 1000
                        table.add_row(('Lambda 1000', 0, 0))
                        table.add_hline()

            self.append(NewPage())

            with self.create(Subsection('Manhattan Plot - pre-QC')):
                with self.create(Figure(position='h!')) as pre_man_images:
                    self.append(Command('centering'))
                    with self.create(
                            SubFigure(position='c', width=NoEscape(r'1\linewidth'))) as man_pre_images:
                        man_pre_images.add_image(man_pre_path, width=NoEscape(r'1\linewidth'))
                with self.create(Figure(position='h!')) as pre_man_images:
                    self.append(Command('centering'))
                    with self.create(
                            SubFigure(position='c', width=NoEscape(r'0.85\linewidth'))) as qq_pre_images:
                        qq_pre_images.add_image(qq_pre_path, width=NoEscape(r'0.85\linewidth'))

            self.append(NewPage())

            with self.create(Subsection('Manhattan Plot - post-QC')):
                with self.create(Figure(position='h!')) as pos_man_images:
                    self.append(Command('centering'))
                    with self.create(
                            SubFigure(position='c', width=NoEscape(r'1\linewidth'))) as man_pos_images:
                        man_pos_images.add_image(man_pos_path, width=NoEscape(r'1\linewidth'))
                with self.create(Figure(position='h!')) as pos_man_images:
                    self.append(Command('centering'))
                    with self.create(
                            SubFigure(position='c', width=NoEscape(r'0.85\linewidth'))) as qq_pos_images:
                        qq_pos_images.add_image(qq_pos_path, width=NoEscape(r'0.85\linewidth'))

    def individual_char(self, id_con_pre_path, id_cas_pre_path, id_con_pos_path, id_cas_pos_path, fstat_fig_path):
        self.append(NewPage())

        with self.create(Section('Per Individual Characteristics Analysis')):
            with self.create(Subsection('Missing Rates - pre-QC')):
                with self.create(Figure(position='h!')) as pre_idcr_images:
                    self.append(Command('centering'))
                    with self.create(
                            SubFigure(position='c', width=NoEscape(r'1\linewidth'))) as idcr_con_pre_images:
                        idcr_con_pre_images.add_image(id_con_pre_path, width=NoEscape(r'1\linewidth'))
                with self.create(Figure(position='h!')) as pre_idcr_images:
                    self.append(Command('centering'))
                    with self.create(
                            SubFigure(position='c', width=NoEscape(r'1\linewidth'))) as idcr_cas_pre_images:
                        idcr_cas_pre_images.add_image(id_cas_pre_path, width=NoEscape(r'1\linewidth'))

            self.append(NewPage())

            with self.create(Subsection('Missing Rates - post-QC')):
                with self.create(Figure(position='h!')) as pos_idcr_images:
                    self.append(Command('centering'))
                    with self.create(
                            SubFigure(position='c', width=NoEscape(r'1\linewidth'))) as idcr_con_pos_images:
                        idcr_con_pos_images.add_image(id_con_pos_path, width=NoEscape(r'1\linewidth'))
                with self.create(Figure(position='h!')) as pos_idcr_images:
                    self.append(Command('centering'))
                    with self.create(
                            SubFigure(position='c', width=NoEscape(r'1\linewidth'))) as idcr_cas_pos_images:
                        idcr_cas_pos_images.add_image(id_cas_pos_path, width=NoEscape(r'1\linewidth'))

            self.append(NewPage())

            with self.create(Subsection('Fstat - Sex Violations')):
                with self.create(Figure(position='h!')) as fstat_image:
                    self.append(Command('centering'))
                    with self.create(
                            SubFigure(position='c', width=NoEscape(r'1\linewidth'))) as fstat_fig:
                        fstat_fig.add_image(fstat_fig_path, width=NoEscape(r'1\linewidth'))

