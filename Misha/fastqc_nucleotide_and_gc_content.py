import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import statistics
import time

start = time.time()


def read_file_and_generate_dictionary_with_fastq_data(input_fastq):
    with open(input_fastq) as r:
        fastq_file = r.readlines()

        sequences_list = []

        for i in range(1, len(fastq_file), 4):
            sequences_list.append(fastq_file[i][:-1])

    return sequences_list


def calculate_gc(sequences_list):
    gc_content_list = []
    for sequence in sequences_list:
        gc_content_list.append(sequence.lower().count('g') + sequence.lower().count('c') * 100 / len(sequence))
    GC_mean = round(statistics.mean(gc_content_list), 2)
    return gc_content_list, GC_mean


def get_nucleotide_per_position(sequences_list):
    nucleotide_per_position_dic = {}
    for sequence in sequences_list:
        for position_number in range(len(sequence)):
            if position_number in nucleotide_per_position_dic:
                nucleotide_per_position_dic[position_number].append(sequence[position_number])
            else:
                nucleotide_per_position_dic[position_number] = [sequence[position_number]]
    return nucleotide_per_position_dic


def get_nucleotide_content(nucleotide_per_position_dic):
    nucleotide_content_per_base = {}
    for position_number in nucleotide_per_position_dic:
        for var in {'A', 'C', 'G', 'T', 'N'}:
            if var in nucleotide_content_per_base:
                nucleotide_content_per_base[var].append(nucleotide_per_position_dic[position_number].count(var) /
                                                        len(nucleotide_per_position_dic[position_number]) * 100)
            else:
                nucleotide_content_per_base[var] = [nucleotide_per_position_dic[position_number].count(var) /
                                                    len(nucleotide_per_position_dic[position_number]) * 100]

    return nucleotide_content_per_base


def make_nucleotide_content_df(nucleotide_content_per_base):
    position_number_list = [i + 1 for i in range(len(nucleotide_content_per_base['A']))]
    nucleotide_content_per_base_df = pd.DataFrame({"Position": position_number_list,
                                                   "A": nucleotide_content_per_base['A'],
                                                   "C": nucleotide_content_per_base['C'],
                                                   "G": nucleotide_content_per_base['G'],
                                                   "T": nucleotide_content_per_base['T']})
    return nucleotide_content_per_base_df


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


def make_N_content_df(nucleotide_content_per_base):
    position_number_list = [i + 1 for i in range(len(nucleotide_content_per_base['A']))]
    N_content_per_base_df = pd.DataFrame({"Position": position_number_list,
                                          "N": nucleotide_content_per_base['N']})
    return N_content_per_base_df


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


def get_GC_theoretical_distribution(gc_content_list):
    GC_mode = statistics.mode(gc_content_list)
    GC_sd = statistics.stdev(gc_content_list)
    number_of_sequences = len(gc_content_list)
    GC_theoretical_distribution = np.random.normal(GC_mode, GC_sd, number_of_sequences)
    return GC_theoretical_distribution


def make_GC_content_integer_rounding(gc_content_list):
    gc_dic_rounded = {}
    for i in range(len(gc_content_list)):
        rounded = round(gc_content_list[i])
        if rounded in gc_dic_rounded:
            gc_dic_rounded[rounded] += 1
        else:
            gc_dic_rounded[rounded] = 1
    return gc_dic_rounded


def make_GC_content_and_GC_theoretical_have_equal_length(gc_theor_dic_rounded, gc_dic_rounded):
    for gc in gc_dic_rounded:
        if gc not in gc_theor_dic_rounded:
            gc_theor_dic_rounded[gc] = 0

    for gc in gc_theor_dic_rounded:
        if gc not in gc_dic_rounded:
            gc_dic_rounded[gc] = 0
    return gc_theor_dic_rounded, gc_dic_rounded


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


def make_nucleotide_and_gc_plots(sequences_list, gc_content_list):
    nucleotide_content_per_base = get_nucleotide_content(get_nucleotide_per_position(sequences_list))

    nucleotide_content_per_base_df = make_nucleotide_content_df(nucleotide_content_per_base)

    make_nucleotide_content_plot(nucleotide_content_per_base_df)

    N_content_per_base_df = make_N_content_df(nucleotide_content_per_base)

    make_N_content_plot(N_content_per_base_df)

    GC_theoretical_distribution = get_GC_theoretical_distribution(gc_content_list)

    gc_theor_dic_rounded = make_GC_content_integer_rounding(GC_theoretical_distribution)

    gc_dic_rounded = make_GC_content_integer_rounding(gc_content_list)

    gc_theor_dic_rounded, gc_dic_rounded = make_GC_content_and_GC_theoretical_have_equal_length(gc_theor_dic_rounded,
                                                                                                gc_dic_rounded)

    GC_content_df = make_GC_content_df(gc_dic_rounded, gc_theor_dic_rounded)

    make_GC_content_plot(GC_content_df)


if __name__ == "__main__":
    input_fastq = '/mnt/d/IB2021_2022/command_line/Project_1/raw_data/amp_res_1.fastq'
    # input_fastq = '/mnt/d/IB2021_2022/command_line/Project_2/SRR1705851.fastq'
    sequences_list = read_file_and_generate_dictionary_with_fastq_data(input_fastq)
    gc_content_list, GC_mean = calculate_gc(sequences_list)
    make_nucleotide_and_gc_plots(sequences_list, gc_content_list)
    print(f'Performing time = {time.time() - start} seconds')
