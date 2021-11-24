import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def make_nucleotide_content_df(nucleotide_content_per_base):
    position_number_list = [i + 1 for i in range(len(nucleotide_content_per_base['A']))]
    nucleotide_content_per_base_df = pd.DataFrame({"Position": position_number_list,
                                                   "A": nucleotide_content_per_base['A'],
                                                   "C": nucleotide_content_per_base['C'],
                                                   "G": nucleotide_content_per_base['G'],
                                                   "T": nucleotide_content_per_base['T']})
    return nucleotide_content_per_base_df


def make_N_content_df(nucleotide_content_per_base):
    position_number_list = [i + 1 for i in range(len(nucleotide_content_per_base['A']))]
    N_content_per_base_df = pd.DataFrame({"Position": position_number_list,
                                          "N": nucleotide_content_per_base['N']})
    return N_content_per_base_df


def make_GC_content_df(gc_dic_rounded, gc_theor_dic_rounded):
    GC_content_rounded_df = pd.DataFrame({"GC content": gc_dic_rounded.keys(),
                                          "Number": gc_dic_rounded.values()})
    GC_content_theoretical_rounded_df = pd.DataFrame({"GC content": gc_theor_dic_rounded.keys(),
                                                      "Number": gc_theor_dic_rounded.values()})

    GC_content_rounded_df_sorted = GC_content_rounded_df.sort_values(by=['GC content'])
    GC_content_theoretical_rounded_df_sorted = GC_content_theoretical_rounded_df.sort_values(by=['GC content'])

    GC_content_df = pd.DataFrame({"GC": list(GC_content_rounded_df_sorted['GC content']),
                                  "GC count per read": list(GC_content_rounded_df_sorted['Number']),
                                  "Theoretical Distribution": list(GC_content_theoretical_rounded_df_sorted['Number'])})
    return GC_content_df


def make_nucleotide_content_plot(nucleotide_content_per_base_df):
    sns.set_style('whitegrid')
    fig, ax = plt.subplots()
    fig.set_size_inches(11.7, 8.27)
    ax.set(xlim=(0, len(nucleotide_content_per_base_df) + 1), ylim=(0, 100))
    palette = ['green', 'blue', 'black', 'red']
    ax.xaxis.set_major_locator(ticker.MultipleLocator(8))
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_major_locator(plt.MaxNLocator(10))
    nucleotide_content_per_base_df_melted = nucleotide_content_per_base_df.melt('Position', var_name='Nucleotide',
                                                                                value_name='vals')
    nucleotide_plot = sns.lineplot(x="Position", y="vals", hue='Nucleotide',
                                   data=nucleotide_content_per_base_df_melted, ax=ax, palette=palette,
                                   linewidth=1, dashes=False)
    nucleotide_plot.set(title='Sequence content across all bases',
                        xlabel='Position in read (bp)', ylabel='Nucleotide content (%)')
    plt.savefig('./nucleotide_content.png', format='png', dpi=300)
#    plt.show()


def make_N_content_plot(N_content_per_base_df):
    sns.set_style('whitegrid')

    fig, ax = plt.subplots()
    fig.set_size_inches(11.7, 8.27)
    ax.set(xlim=(0, len(N_content_per_base_df) + 1), ylim=(0, 100))
    palette = ['red']
    ax.xaxis.set_major_locator(ticker.MultipleLocator(8))
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_major_locator(plt.MaxNLocator(10))
    N_content_per_base_df_melted = N_content_per_base_df.melt('Position', var_name='Nucleotide', value_name='vals')
    N_plot = sns.lineplot(x="Position", y="vals", hue='Nucleotide',
                          data=N_content_per_base_df_melted, ax=ax, palette=palette,
                          linewidth=1, dashes=False)
    N_plot.set(title='N content across all bases', xlabel='Position in read (bp)', ylabel='Nucleotide content (%)')
    plt.savefig('./N_content.png', format='png', dpi=300)
#    plt.show()


def make_GC_content_plot(GC_content_df):
    sns.set_style('whitegrid')
    fig, ax = plt.subplots()
    fig.set_size_inches(11.7, 8.27)
    ax.set(xlim=(1, 100))
    palette = ['red', 'blue']
    ax.xaxis.set_major_locator(ticker.MultipleLocator(4))
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    GC_content_df_melted = GC_content_df.melt('GC', var_name='GC content', value_name='vals')
    GC_content_plot = sns.lineplot(x="GC", y="vals", hue='GC content',
                                   data=GC_content_df_melted, ax=ax, palette=palette,
                                   linewidth=1, dashes=False)
    GC_content_plot.set(title='GC distribution over all sequences', xlabel='Mean GC content (%)',
                        ylabel='Number of reads')
    plt.savefig('./GC_content.png', format='png', dpi=300)
#    plt.show()


def make_nucleotide_and_gc_plots(nucleotide_content_per_base, gc_dic_rounded, gc_theor_dic_rounded):
    nucleotide_content_per_base_df = make_nucleotide_content_df(nucleotide_content_per_base)

    N_content_per_base_df = make_N_content_df(nucleotide_content_per_base)

    GC_content_df = make_GC_content_df(gc_dic_rounded, gc_theor_dic_rounded)

    make_nucleotide_content_plot(nucleotide_content_per_base_df)

    make_N_content_plot(N_content_per_base_df)

    make_GC_content_plot(GC_content_df)
