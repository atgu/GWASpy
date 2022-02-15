import hail as hl
from hail import Table
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plt_hist(expression: hl.Expression, bins: int = 50, range: list = None, threshold: float = None,
             title: str = None, x_label: str = None, y_label: str = None, log: bool = False, figsize: tuple = (12, 8)):
    exprs = expression.collect()
    df = pd.DataFrame({'col': exprs})

    title = f'{title}' if title else ''
    fig = plt.figure(figsize=figsize)
    if log is True:
        df = df[df['col'] != 0]
        plt.hist(np.log10(df.col), edgecolor='black', density=False, bins=bins, color='tab:blue')
    else:
        plt.hist(df['col'], edgecolor='black', density=False, bins=bins, color='tab:blue')
    if threshold:
        plt.axvline(x=threshold, color='red', linestyle='--')
    if range is not None:
        plt.xlim(xmin=range[0], xmax=range[1])
    plt.title(title, fontsize=20)
    plt.ylabel(y_label if y_label else 'Frequency', fontsize=15)
    plt.xlabel(x_label if x_label else '', fontsize=15)
    plt.close()

    return fig


def fstat_plot(df_female, df_male, f_stat_x: float = 0.4, f_stat_y: float = 0.8, figsize: tuple = (12, 8)):
    fig, axs = plt.subplots(2, figsize=figsize)
    axs[0].hist(df_female['filters'], bins=30, histtype='bar', alpha=0.8, fill=True, color='tab:blue', edgecolor="k")
    axs[0].axvline(x=f_stat_y, color='red', linestyle='--')
    axs[0].set_title('Female Fstat', fontsize=20)
    axs[0].set_xlabel(xlabel='Fstat', fontsize=15)
    axs[0].set_ylabel(ylabel='Frequency', fontsize=15)
    axs[0].tick_params(axis='both', which='major', labelsize=12)
    axs[1].hist(df_male['filters'], bins=40, histtype='bar', alpha=0.8, fill=True, color='tab:blue', edgecolor="k")
    axs[1].axvline(x=f_stat_x, color='red', linestyle='--')
    axs[1].set_title('Male Fstat', fontsize=20)
    axs[1].set_xlabel(xlabel='Fstat', fontsize=15)
    axs[1].set_ylabel(ylabel='Frequency', fontsize=15)
    axs[1].tick_params(axis='both', which='major', labelsize=12)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.close()

    return fig


def qqplot(pvals, title: str = None, figsize: tuple = (10, 10)):
    source = pvals._indices.source
    if isinstance(source, Table):
        ht = source.select(p_value=pvals)
    else:
        ht = source.select_rows(p_value=pvals).rows()

    ht = ht.key_by().select('p_value').key_by('p_value').persist()
    lambda_gc = hl.lambda_gc(ht['p_value'])
    n = ht.count()
    ht = ht.annotate(
        observed_p=-hl.log10(ht['p_value']),
        expected_p=-hl.log10((hl.scan.count() + 1) / n),
        p_val=ht['p_value']
    ).persist()

    p_val_pd = ht.to_pandas()
    p_val_pd['observed_p'].values[p_val_pd['observed_p'] > 10] = 10
    mini = min(p_val_pd['expected_p'].max(), p_val_pd['observed_p'].max())
    maxi = max(p_val_pd['expected_p'].max(), p_val_pd['observed_p'].max())

    title = f'{title}' if title else 'QQ Plot'

    fig = plt.figure(figsize=figsize)
    plt.scatter(p_val_pd['expected_p'], p_val_pd['observed_p'], c='black', s=0.5)
    plt.plot((0, mini), (0, mini), 'red')
    plt.xlim([0, maxi + 0.5])
    plt.ylim([0, maxi + 0.5])
    plt.title(title, fontsize=20)
    plt.ylabel('Observed -log10(' + r'$p$' + ')', fontsize=15)
    plt.xlabel('Expected -log10(' + r'$p$' + ')', fontsize=15)
    plt.close()

    return fig, round(lambda_gc, 3)


def manhattan_plot(pvals, significance_threshold: float = -np.log10(5E-08), title: str = None,
                   figsize: tuple = (17, 11), annotate_sig: bool = False):
    source = pvals._indices.source

    if isinstance(source, Table):
        ht = source.select(p_value=pvals)
    else:
        ht = source.select_rows(p_value=pvals).rows()

    data = ht.to_pandas()
    data = data.drop('alleles', 1)  # remove the 'allele' column
    data['locus'] = data['locus'].astype(str)
    data[['CHROM', 'POS']] = data.locus.str.split(":", expand=True)
    data.columns = ['locus', 'p', 'chromosome', 'position'] # rename columns
    data['position'] = data['position'].astype(int)
    data['chromosome'].replace({"X": 23, "Y": 24, "MT": 25}, inplace=True)
    data['chromosome'] = data['chromosome'].astype(int)
    data.dropna(subset=['p'], inplace=True)  # drop NAs as log10(val) won't work

    data['-log10(p_value)'] = -np.log10(data['p'])  # compute log10(pvals)
    data['-log10(p_value)'].values[data['-log10(p_value)'] > 199] = 199
    data['chromosome'] = data['chromosome'].astype('category')
    data['ind'] = range(len(data))
    data_grouped = data.groupby('chromosome')

    title = f'{title}' if title else 'Manhattan Plot'

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot()

    x_labels = []
    x_labels_pos = []

    colors = ['tab:orange', 'tab:blue']

    for num, (name, group) in enumerate(data_grouped):
        group.plot(kind='scatter', x='ind', y='-log10(p_value)', marker='o', color=colors[int(name) % len(colors)],
                   ax=ax, s=10)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0]) / 2))

    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    # ax.set_xlim([0, len(data)])
    ax.margins(0.05)
    # ax.set_ylim([0, data['-log10(p_value)'].max() + 1])
    ax.set_xlabel('Chromosome', fontsize=15)
    ax.set_ylabel('-log10(p_value)', fontsize=15)
    plt.title(title, fontsize=20)
    plt.axhline(y=significance_threshold, color='red', linestyle='--', linewidth=2)
    plt.xticks(fontsize=10, rotation=90)
    plt.yticks(fontsize=10)
    if annotate_sig is True:
        for index, row in data.iterrows():
            if row['-log10(p_value)'] >= significance_threshold:
                ax.annotate('{}:{}'.format(row['chromosome'], row['position']),
                            xy=(index, row['-log10(p_value)'] + 0.1),
                            bbox=dict(boxstyle="round", fc="0.8"))
    plt.close()

    return fig
