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
