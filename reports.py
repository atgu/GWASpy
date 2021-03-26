__author__ = 'Lindo Nkambule'

from pylatex import Document, Section, Subsection, Command, Center, Tabular, NewPage, Figure, SubFigure, TextColor
from pylatex.utils import NoEscape, bold


class MyDocument(Document):
    def __init__(self, basename):
        super().__init__()

        self.preamble.append(Command('title', 'QC Report of {}'.format(basename)))
        self.preamble.append(Command('author', 'Anonymous author'))
        self.preamble.append(Command('date', NoEscape(r'\today')))
        self.append(NoEscape(r'\maketitle'))

    def flags_table(self, pre_qc_counts, pos_qc_counts, results, lambda_gc, sig_vars):
        cas_con_ratio = round(pos_qc_counts['is_case_counts']['case'] / pos_qc_counts['is_case_counts']['control'], 4)
        nids_lost = (pre_qc_counts['n_samples'] - pos_qc_counts['n_samples']) / pre_qc_counts['n_samples']
        nids_sex_check = results['sex_warnings'][True] / pre_qc_counts['n_samples']

        tbl = [['nsnps-postqc', pos_qc_counts['n_variants'], 250000, 200000, 0, 'green'],
               ['ncases-postqc', pos_qc_counts['is_case_counts']['case'], 100, 50, 0, 'green'],
               ['ncontrols-postqc', pos_qc_counts['is_case_counts']['control'], 100, 50, 0, 'green'],
               ['case-control-ratio-postqc', cas_con_ratio, 0.0625, 0.0278, 0, 'green'],
               ['nids-lost-ratio', nids_lost, 0.01, 0.1, 0, 'green'],
               ['n-nopt-postqc', pos_qc_counts['is_case_counts']['unknown'], 0, 10, 0, 'green'],
               ['nids-sexcheck-ratio', nids_sex_check, 0.005, 0.0025, 0, 'green'],
               ['lambda-postqc', lambda_gc, 1.1, 1.2, 0, 'green'],
               ['nsnps-gws', sig_vars, 0, 1, 0, 'green']]

        for i in tbl:
            if (i[0] == 'nsnps-postqc') | (i[0] == 'ncases-postqc') | (i[0] == 'ncontrols-postqc') | (
                    i[0] == 'case-control-ratio-postqc'):
                if (i[1] <= i[2]) & (i[1] >= i[3]):
                    i[4] = 1
                    i[5] = 'orange'
                elif i[1] < i[3]:
                    i[4] = 2
                    i[5] = 'red'
                else:
                    i[4] = 0
                    i[5] = 'green'

            if i[0] == 'nids-lost-ratio':
                if (i[1] > i[2]) & (i[1] <= i[3]):
                    i[4] = 1
                    i[5] = 'orange'
                elif i[1] > i[3]:
                    i[4] = 2
                    i[5] = 'red'
                else:
                    i[4] = 0
                    i[5] = 'green'

            if (i[0] == 'nids-sexcheck-ratio') | (i[0] == 'lambda-postqc'):
                if (i[1] >= i[2]) & (i[1] <= i[3]):
                    i[4] = 1
                    i[5] = 'orange'
                elif i[1] > i[3]:
                    i[4] = 2
                    i[5] = 'red'
                else:
                    i[4] = 0
                    i[5] = 'green'

            if (i[0] == 'nsnps-gws') | (i[0] == 'n-nopt-postqc'):
                if i[1] == i[2]:
                    i[4] = 0
                    i[5] = 'green'
                elif (i[1] > i[2]) & (i[1] <= i[3]):
                    i[4] = 1
                    i[5] = 'orange'
                else:
                    i[4] = 2
                    i[5] = 'red'

        with self.create(Section('Flags')):
            with self.create(Center()) as centered:
                with centered.create(Tabular('|l|c|c|c|c|c|')) as table:
                    table.add_hline()
                    table.add_row((bold('Flagname'), bold('Value'), bold('yellow-th'), bold('red-th'), bold('flag'),
                                   bold('color')))
                    table.add_hline()
                    table.add_row(tbl[0][0], tbl[0][1], tbl[0][2], tbl[0][3], tbl[0][4], TextColor(tbl[0][5],tbl[0][5]))
                    table.add_hline()
                    table.add_row(tbl[1][0], tbl[1][1], tbl[1][2], tbl[1][3], tbl[1][4], TextColor(tbl[1][5],tbl[1][5]))
                    table.add_hline()
                    table.add_row(tbl[2][0], tbl[2][1], tbl[2][2], tbl[2][3], tbl[2][4], TextColor(tbl[2][5],tbl[2][5]))
                    table.add_hline()
                    table.add_row(tbl[3][0], tbl[3][1], tbl[3][2], tbl[3][3], tbl[3][4], TextColor(tbl[3][5],tbl[3][5]))
                    table.add_hline()
                    table.add_row(tbl[4][0], tbl[4][1], tbl[4][2], tbl[4][3], tbl[4][4], TextColor(tbl[4][5],tbl[4][5]))
                    table.add_hline()
                    table.add_row(tbl[5][0], tbl[5][1], tbl[5][2], tbl[5][3], tbl[5][4], TextColor(tbl[5][5],tbl[5][5]))
                    table.add_hline()
                    table.add_row(tbl[6][0], tbl[6][1], tbl[6][2], tbl[6][3], tbl[6][4], TextColor(tbl[6][5],tbl[6][5]))
                    table.add_hline()
                    table.add_row(tbl[7][0], tbl[7][1], tbl[7][2], tbl[7][3], tbl[7][4], TextColor(tbl[7][5],tbl[7][5]))
                    table.add_hline()
                    table.add_row(tbl[8][0], tbl[8][1], tbl[8][2], tbl[8][3], tbl[8][4],
                                  TextColor(tbl[8][5], tbl[8][5]))
                    table.add_hline()

    def general_info(self, pre_qc_conts, post_qc_conts, count_results, pre_filter, id_cr, fhet_thresh, var_cr,
                     miss_diff, hwe_con, hwe_cas, data_type):
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

        self.append(NewPage())

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
                        if data_type == "Case-only":
                            table.add_row(('SNPs: HWE-cases < {}'.format(hwe_cas), count_results['hwe_cas'][True]))
                            table.add_hline()
                        if data_type == "Control-only":
                            table.add_row(('SNPs: HWE-controls < {}'.format(hwe_con), count_results['hwe_con'][True]))
                            table.add_hline()
                        if data_type == "Case-Control":
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
                        table.add_row(('Lambda 1000', table_results[4], table_results[5]))
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

    def individual_char(self, id_con_pre_path, id_cas_pre_path, fstat_fig_path, data_type):
        self.append(NewPage())

        with self.create(Section('Per Individual Characteristics Analysis')):
            with self.create(Subsection('Missing Rates - pre-QC')):
                if data_type == "Case-only":
                    with self.create(Figure(position='h!')) as pre_idcr_images:
                        self.append(Command('centering'))
                        with self.create(
                                SubFigure(position='c', width=NoEscape(r'1\linewidth'))) as idcr_cas_pre_images:
                            idcr_cas_pre_images.add_image(id_cas_pre_path, width=NoEscape(r'1\linewidth'))
                    self.append(NewPage())

                if data_type == "Control-only":
                    with self.create(Figure(position='h!')) as pre_idcr_images:
                        self.append(Command('centering'))
                        with self.create(
                                SubFigure(position='c', width=NoEscape(r'1\linewidth'))) as idcr_con_pre_images:
                            idcr_con_pre_images.add_image(id_con_pre_path, width=NoEscape(r'1\linewidth'))
                    self.append(NewPage())

                if data_type == "Case-Control":
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

            with self.create(Subsection('Fstat - Sex Violations')):
                with self.create(Figure(position='h!')) as fstat_image:
                    self.append(Command('centering'))
                    with self.create(
                            SubFigure(position='c', width=NoEscape(r'1\linewidth'))) as fstat_fig:
                        fstat_fig.add_image(fstat_fig_path, width=NoEscape(r'1\linewidth'))

    def snp_char(self, var_con_pre_path, var_cas_pre_path, data_type):
        self.append(NewPage())

        with self.create(Section('Per SNP Characteristics Analysis')):
            with self.create(Subsection('Missing Rates - pre-QC')):
                if data_type == "Case-only":
                    with self.create(Figure(position='h!')) as pre_varcr_images:
                        self.append(Command('centering'))
                        with self.create(
                                SubFigure(position='c', width=NoEscape(r'1\linewidth'))) as varcr_cas_pre_images:
                            varcr_cas_pre_images.add_image(var_cas_pre_path, width=NoEscape(r'1\linewidth'))
                    self.append(NewPage())

                if data_type == "Control-only":
                    with self.create(Figure(position='h!')) as pre_varcr_images:
                        self.append(Command('centering'))
                        with self.create(
                                SubFigure(position='c', width=NoEscape(r'1\linewidth'))) as varcr_con_pre_images:
                            varcr_con_pre_images.add_image(var_con_pre_path, width=NoEscape(r'1\linewidth'))
                    self.append(NewPage())

                if data_type == "Case-Control":
                    with self.create(Figure(position='h!')) as pre_varcr_images:
                        self.append(Command('centering'))
                        with self.create(
                                SubFigure(position='c', width=NoEscape(r'1\linewidth'))) as varcr_con_pre_images:
                            varcr_con_pre_images.add_image(var_con_pre_path, width=NoEscape(r'1\linewidth'))
                    with self.create(Figure(position='h!')) as pre_varcr_images:
                        self.append(Command('centering'))
                        with self.create(
                                SubFigure(position='c', width=NoEscape(r'1\linewidth'))) as varcr_cas_pre_images:
                            varcr_cas_pre_images.add_image(var_cas_pre_path, width=NoEscape(r'1\linewidth'))
                    self.append(NewPage())

