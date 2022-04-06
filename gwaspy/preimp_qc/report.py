__author__ = 'Lindo Nkambule'

from pylatex import Document, Section, Subsection, Command, Center, Tabular, NewPage, Figure, SubFigure, TextColor
from pylatex.utils import NoEscape, bold


class MyDocument(Document):
    def __init__(self, basename):
        super().__init__()

        self.preamble.append(Command('title', 'QC Report of {}'.format(basename)))
        # self.preamble.append(Command('author', 'Anonymous author'))
        self.preamble.append(Command('date', NoEscape(r'\today')))
        self.append(NoEscape(r'\maketitle'))

    def flags_table(self, pre_qc_counts=None, pos_qc_counts=None, results=None, lambda_gc=None, sig_vars=None):
        nids_lost = (pre_qc_counts['n_samples'] - pos_qc_counts['n_samples']) / pre_qc_counts['n_samples']
        nids_lost = round(nids_lost, 4)
        nids_sex_check = results['sex_warnings'][True] / pre_qc_counts['n_samples']
        nids_sex_check = round(nids_sex_check, 4)

        if 'is_case_counts' in pre_qc_counts.keys() | pos_qc_counts.keys():
            cas_con_ratio = round(pos_qc_counts['is_case_counts']['case'] / pos_qc_counts['is_case_counts']['control'], 4)

            tbl = [['Number of SNPs Post-QC', pos_qc_counts['n_variants'], 250000, 200000, 0, 'green'],
                   ['Number of Cases Post-QC', pos_qc_counts['is_case_counts']['case'], 100, 50, 0, 'green'],
                   ['Number of Controls Post-QC', pos_qc_counts['is_case_counts']['control'], 100, 50, 0, 'green'],
                   ['Case-Control ratio Post-QC', cas_con_ratio, 0.0625, 0.0278, 0, 'green'],
                   ['Number of IDs lost ratio', nids_lost, 0.01, 0.1, 0, 'green'],
                   ['Number of IDs with no Phenotype Post-QC', pos_qc_counts['is_case_counts']['unknown'], 0, 10, 0, 'green'],
                   ['Ratio of IDs that failed sex checks', nids_sex_check, 0.005, 0.025, 0, 'green'],
                   ['Lambda GC Post-QC', lambda_gc, 1.1, 1.2, 0, 'green'],
                   ['Number of Significant GWAS hits Post-QC', sig_vars, 0, 1, 0, 'green']]

        else:
            tbl = [['Number of SNPs Post-QC', pos_qc_counts['n_variants'], 250000, 200000, 0, 'green'],
                   ['Number of IDs lost ratio', nids_lost, 0.01, 0.1, 0, 'green'],
                   ['Ratio of IDs that failed sex checks', nids_sex_check, 0.005, 0.025, 0, 'green']]

        for i in tbl:
            if (i[0] == 'Number of SNPs Post-QC') | (i[0] == 'Number of Cases Post-QC') |\
                    (i[0] == 'Number of Controls Post-QC') | (i[0] == 'Case-Control ratio Post-QC'):
                if (i[1] <= i[2]) & (i[1] >= i[3]):
                    i[4] = 1
                    i[5] = 'orange'
                elif i[1] < i[3]:
                    i[4] = 2
                    i[5] = 'red'
                else:
                    i[4] = 0
                    i[5] = 'green'

            if i[0] == 'Number of IDs lost ratio':
                if (i[1] > i[2]) & (i[1] <= i[3]):
                    i[4] = 1
                    i[5] = 'orange'
                elif i[1] > i[3]:
                    i[4] = 2
                    i[5] = 'red'
                else:
                    i[4] = 0
                    i[5] = 'green'

            if (i[0] == 'Ratio of IDs that failed sex checks') | (i[0] == 'Lambda GC Post-QC'):
                if (i[1] >= i[2]) & (i[1] <= i[3]):
                    i[4] = 1
                    i[5] = 'orange'
                elif i[1] > i[3]:
                    i[4] = 2
                    i[5] = 'red'
                else:
                    i[4] = 0
                    i[5] = 'green'

            if (i[0] == 'Number of Significant GWAS hits Post-QC') | (i[0] == 'Number of IDs with no Phenotype Post-QC'):
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
                    for i in range(0, len(tbl)):
                        table.add_row(tbl[i][0], tbl[i][1], tbl[i][2], tbl[i][3], tbl[i][4],
                                      TextColor(tbl[i][5], tbl[i][5]))
                        table.add_hline()

    def general_info(self, pre_qc_conts, post_qc_conts, count_results, pre_filter, id_cr, fhet_thresh, var_cr,
                     miss_diff, hwe_con, hwe_cas, hwe_all, data_type):
        global pre_pheno_counts, post_pheno_counts, pheno_diffs
        if 'is_case_counts' in pre_qc_conts.keys() | post_qc_conts.keys():
            pre_pheno_counts = str(pre_qc_conts['is_case_counts']['case']) + ', ' + \
                               str(pre_qc_conts['is_case_counts']['control']) + ', ' + \
                               str(pre_qc_conts['is_case_counts']['unknown'])

            post_pheno_counts = str(post_qc_conts['is_case_counts']['case']) + ', ' + \
                                str(post_qc_conts['is_case_counts']['control']) + ', ' + \
                                str(post_qc_conts['is_case_counts']['unknown'])

            pheno_diffs = str(pre_qc_conts['is_case_counts']['case'] - post_qc_conts['is_case_counts']['case']) + ', ' \
                          + str(pre_qc_conts['is_case_counts']['control'] - post_qc_conts['is_case_counts']['control'])\
                          + ', ' +\
                          str(pre_qc_conts['is_case_counts']['unknown'] - post_qc_conts['is_case_counts']['unknown'])

        pre_sex_counts = str(pre_qc_conts['is_female_counts']['male']) + ', ' + \
                         str(pre_qc_conts['is_female_counts']['female']) + ', ' + \
                         str(pre_qc_conts['is_female_counts']['unknown'])

        post_sex_counts = str(post_qc_conts['is_female_counts']['male']) + ', ' + \
                          str(post_qc_conts['is_female_counts']['female']) + ', ' + \
                          str(post_qc_conts['is_female_counts']['unknown'])

        sex_diffs = str(pre_qc_conts['is_female_counts']['male']-post_qc_conts['is_female_counts']['male']) + ', '\
                    + str(pre_qc_conts['is_female_counts']['female']-post_qc_conts['is_female_counts']['female']) +\
                    ', ' + str(pre_qc_conts['is_female_counts']['unknown']-post_qc_conts['is_female_counts']['unknown'])

        with self.create(Section('General Info')):
            # with self.create(Subsection('Size of sample')):
            with self.create(Center()) as centered:
                with centered.create(Tabular('|c|c|c|c|')) as table:
                    table.add_hline()
                    table.add_row((bold('Test'), bold('pre QC'), bold('post QC'), bold('exlcusion-N')))
                    table.add_hline()
                    if 'is_case_counts' in pre_qc_conts.keys() | post_qc_conts.keys():
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

            # with self.create(Subsection('Exclusion overview')):
            with self.create(Center()) as centered:
                with centered.create(Tabular('|l|l|')) as table:
                    table.add_hline()
                    table.add_row((bold('Filter'), bold('N ')))
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
                    if data_type != "no-pheno":
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
                    if data_type == "no-pheno":
                        table.add_row(('SNPs: HWE < {}'.format(hwe_all), count_results['hwe_all'][True]))
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

    def individual_char(self, id_con_pre_path, id_cas_pre_path, id_all_path, fstat_fig_path, data_type):
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

                if data_type == "no-pheno":
                    with self.create(Figure(position='h!')) as pre_idcr_images:
                        self.append(Command('centering'))
                        with self.create(
                                SubFigure(position='c', width=NoEscape(r'1\linewidth'))) as idcr_all_pre_images:
                            idcr_all_pre_images.add_image(id_all_path, width=NoEscape(r'1\linewidth'))
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

    def snp_char(self, var_con_pre_path, var_cas_pre_path, var_all_path, data_type):
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

                if data_type == "no-pheno":
                    with self.create(Figure(position='h!')) as pre_varcr_images:
                        self.append(Command('centering'))
                        with self.create(
                                SubFigure(position='c', width=NoEscape(r'1\linewidth'))) as varcr_all_pre_images:
                            varcr_all_pre_images.add_image(var_all_path, width=NoEscape(r'1\linewidth'))
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

