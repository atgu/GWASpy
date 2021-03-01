import hail as hl
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
    plt.title(title)
    plt.ylabel(y_label if y_label else 'Frequency')
    plt.xlabel(x_label if x_label else '')
    plt.close()

    return fig


def fstat_plot(imputed_sex_ht: hl.Table, f_stat_x: float = 0.4, f_stat_y: float = 0.8, figsize: tuple = (12, 8)):
    fstat_df = imputed_sex_ht.to_pandas()

    fstat_df['is_female'] = fstat_df['is_female'].astype(str)
    fstat_df['is_female'] = fstat_df['is_female'].replace(['True', 'False', 'None'], ['female', 'male', 'unspecified'])

    fig = plt.figure(figsize=figsize)
    plt.hist(fstat_df['f_stat'],
             bins=40,
             histtype='bar',
             alpha=0.8,
             fill=True,
             color='tab:blue',
             edgecolor="k")
    plt.axvline(x=f_stat_y, color='red', linestyle='--')
    plt.axvline(x=f_stat_x, color='red', linestyle='--')
    plt.close()

    return fig
