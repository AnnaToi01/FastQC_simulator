import matplotlib.pyplot as plt
import math
import numpy as np
from pathlib import Path
import seaborn as sns
from PIL import Image
import imgkit
import matplotlib.ticker as mticker
import plotly.graph_objects as go
import pandas as pd


# Anna
def plot_basic_statistics(basicstats, output_dir):
    """
    Styles and saves basic statistics dataframe as .png
    :param basicstats: pd Dataframe, basic statistics
    :param output_dir: directory where picture will be stored
    :return: None, save basic statistics table as .png
    """
    basicstats_styled = basicstats.style.set_table_styles(
        [{
            'selector': 'th.col_heading',
            'props': [('border', 'solid white'), ('background-color', '#000080'), ('text-align', 'center'),
                      ("font-family", "sans-serif"), ("color", "white")]
        },
            {'selector': 'td',

             'props': [('border', 'solid white'), ('background-color', '#EEEEEE'), ('text-align', 'left'),
                       ("font-family", "monospace"), ("padding", "0.4em"), ("color", "#000000")]
             }
        ]).hide_index()

    html = basicstats_styled.render()
    imgkit.from_string(html, Path(output_dir, 'basic_statistics.png'))
    return


# Ivan
def draw_boxplot_with_qscores(df, phred, encoding, stats_for_boxplots, means, output_dir):
    """
    Draw boxplots with quality scores per position
    :param df: data frame with unique quality scores per position and their counts
    :param phred: int, type of phred quality (33 or 64)
    :param encoding: type of Illumina
    :param stats_for_boxplots: dict with statistics of quality distribution per position
    :param means: numpy.array with quality means per position
    :param output_dir: directory where picture will be stored
    """
    positions = [col_name for col_name in df.columns if col_name not in ['ascii_symb', 'quality']]
    y_max = max([stats['whishi'] for stats in stats_for_boxplots])
    _, ax = plt.subplots(figsize=(14, 10))
    bplot = ax.bxp(stats_for_boxplots, showfliers=False, patch_artist=True, widths=0.7, zorder=2)
    for box, med in zip(bplot['boxes'], bplot['medians']):
        box.set_facecolor('yellow')
        box.set_linewidth(1)
        med.set_color('red')
        med.set_linewidth(2)
    ax.set_xlabel('Position in read (bp)', fontsize=14)
    ax.set_title(f'Quality scores across all bases ({encoding})', fontsize=14)
    ax.xaxis.grid(True, zorder=1)
    length = len(positions)
    if length <= 150:
        xlabels = np.array(positions[::2])
        xticks = np.arange(1, length + 1, 2)
    else:
        xlabels = np.array(positions[::3])
        xticks = np.arange(1, length + 1, 3)
    ylabels = yticks = np.arange(0, 42, 2)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, rotation=45)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.set_xlim([0, length + 1])
    ax.set_ylim([0, y_max + 0.5])
    ax.plot(np.arange(1, len(means) + 1), means, color='blue', zorder=3)
    ax.axhspan(28, 42, facecolor='green', alpha=0.2)
    ax.axhspan(20, 28, facecolor='orange', alpha=0.2)
    ax.axhspan(0, 20, facecolor='red', alpha=0.2)
    plt.savefig(Path(output_dir, 'per_base_quality.png'))


# Ivan
def draw_per_sequence_quality_scores(read_quality_dict, output_dir):
    """
    Draw Quality score distribution over all sequences
    :param read_quality_dict: dict, where keys - quality scores per read, values - their counts
    :param output_dir: directory where picture will be stored
    """
    q_scores_unique = list(read_quality_dict.keys())
    counts = list(read_quality_dict.values())
    min_quality = min(q_scores_unique)
    max_quality = max(q_scores_unique)
    plt.figure(figsize=(14, 10))
    sns.lineplot(q_scores_unique, counts, color='red', label='Average Quality per read')
    plt.grid(axis="y")
    plt.xlabel('Mean Sequence Quality (Phred Score)', fontsize=14)
    plt.title('Quality score distribution over all sequences', fontsize=14)
    plt.xlim([min_quality - 1, max_quality + 2])
    plt.ylim([0, None])
    xticks = np.arange(min_quality, max_quality + 2, 2)
    plt.xticks(xticks)
    for x in xticks[::2]:
        plt.axvspan(x + 1, x + 3, color='black', alpha=0.1, zorder=1)
    leg = plt.legend(loc="upper right", prop={'size': 14}, handlelength=0, borderaxespad=0)
    for line, text in zip(leg.get_lines(), leg.get_texts()):
        text.set_color(line.get_color())
    plt.savefig(Path(output_dir, 'per_sequence_quality.png'))


# Anna
def plot_tile_sequence(tile_quality, max_length, output_dir):
    """
    Plots per tile sequence quality
    :param tile_quality: dic, tile_number: array(average quality per position)
    :param output_dir: directory where picture will be stored
    :return: None, saves the generated figure and shows it
    """
    if tile_quality is None:
        return

    def makeColors():
        """
        Generates color palette
        :return:
            colors: list, color palette
        """
        min_num = 0 - 50 ** 0.5
        max_num = (99 - 50) ** 0.5
        colors = [0] * 100
        for c in range(100):
            actualC = c - 50
            if actualC < 0:
                actualC = 0 - actualC
            corrected = actualC ** 0.5
            if c < 50 and corrected > 0:
                corrected = 0 - corrected
            r, g, b = getRGB(corrected, min_num, max_num)
            colors[c] = [r, g, b]
        return colors

    def getColor(value, min_num, max_num):
        """
        Gets the color of the specific tile
        :param value: average quality value
        :param min_num: int, 0
        :param max_num: int, 10
        :return:
                color: list, of the corresponding tile
        """
        colors = makeColors()
        percentage = int(100 * (value - min_num) / (max_num - min_num))
        if percentage > 100:
            percentage = 100
        if percentage < 1:
            percentage = 1

        return colors[percentage - 1]

    def getRGB(value, min_num, max_num):
        """
        Generates tuple of RGV
        :param value: color value
        :param min_num: minimal value, 0 - 50 ** 0.5
        :param max_num: maximal value, (99 - 50) ** 0.5
        :return: tuple, RGB: r - red, g - green, b - blue
        """
        diff = max_num - min_num
        if value < (min_num + diff * 0.25):
            red = 0
            blue = 200
            green = int(200 * ((value - min_num) / (diff * 0.25)))
        elif value < (min_num + diff * 0.5):
            red = 0
            green = 200
            blue = int(200 - (200 * ((value - (min_num + (diff * 0.25))) / (diff * 0.25))))
        elif value < (min_num + diff * 0.75):
            green = 200
            blue = 0
            red = int(200 * ((value - (min_num + (diff * 0.5))) / (diff * 0.25)))
        else:
            red = 200
            blue = 0
            green = int(200 - (200 * ((value - (min_num + (diff * 0.75))) / (diff * 0.25))))
        return red, green, blue

    _ = plt.figure(figsize=(14, 10))
    ax = plt.subplot(111)

    for i, nmtile in enumerate(tile_quality.keys()):
        # length = np.count_nonzero(~np.isnan(tile_quality[nmtile]))
        for j in range(max_length):
            color = getColor(-tile_quality[nmtile][j], 0, 10)
            im = Image.new('RGB', (1, 1), tuple(color))
            plt.imshow(im, aspect='equal', extent=np.array([j, j + 1, i, i + 1]), interpolation="nearest")

    plt.xlim(0, max_length)
    plt.ylim(0, len(tile_quality.keys()))
    plt.xticks(np.arange(1, max_length + 1, 5))
    labels = list(tile_quality.keys())[::2]
    ax.set_yticks(np.arange(len(tile_quality.keys()))[::2])
    ax.set_yticklabels(labels)
    plt.xlabel('Position in read (bp)', fontsize=14)
    plt.title('Quality per tile', fontsize=14)
    plt.savefig(Path(output_dir, 'per_tile_sequence_quality.png'))


# Anna
def plot_length_distribution(length_distribution_dic, output_dir):
    """
    Plot Sequence length distribution
    :param length_distribution_dic: dictionary, read_length:count
    :param output_dir: directory where picture will be stored
    :return: plot
    """
    x = [int(key) for key in length_distribution_dic.keys()]
    min_x = min(x) - 1
    max_x = max(x) + 1
    x.append(min_x)
    x.append(max_x)
    y = [int(value) for value in length_distribution_dic.values()]
    y.append(0)
    y.append(0)
    max_y = max(y)
    x, y = zip(*sorted(zip(x, y)))
    _, ax = plt.subplots(figsize=(14, 10))
    plt.plot(x, y, label="Sequence length", color='r')
    plt.grid(axis="x")
    if max_x - min_x < 20:
        step = 1
    else:
        step = math.ceil((max_x - min_x) / 10)
    xticks = np.arange(min_x, max_x + step, step=step)
    for x0, x1 in zip(xticks[::2], xticks[1::2]):
        plt.axvspan(x0 + step / 2, x1 + step / 2, color='black', alpha=0.1, zorder=0)
    plt.xticks(xticks)
    difference_yticks = max_y // 10
    digits_zero = 1
    while difference_yticks // 10 > 1:
        difference_yticks = difference_yticks // 10
        digits_zero *= 10
    step_y = math.floor(max_y / difference_yticks / digits_zero) * digits_zero
    plt.yticks(np.arange(0, math.ceil(max(y) / 1000) * 1000 + step_y, step=step_y))
    ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())

    ax.ticklabel_format(style='plain', axis='y')
    leg = plt.legend(handlelength=0, loc="upper right", borderaxespad=0, prop={'size': 14})
    for line, text in zip(leg.get_lines(), leg.get_texts()):
        text.set_color(line.get_color())
    plt.title("Distribution of sequence lengths over all sequences", fontsize=14)
    plt.xlabel("Sequence length", fontsize=14)
    plt.savefig(Path(output_dir, 'sequence_length_distribution.png'))
    # plt.show()


# Misha
def make_nucleotide_content_df(nucleotide_content_per_base):
    """
    Make Nucleotide content pandas DataFrame for "A", "C", "G" and "T"
    :param nucleotide_content_per_base: dic, Nucleotide: Nucleotide content per position in read (%)
    :return:
            nucleotide_content_per_base_df: pandas DataFrame, Position in read (bp), Nucleotide content per position (%)
    """
    position_number_list =\
        [i + 1 for i in range(len(nucleotide_content_per_base[list(nucleotide_content_per_base.keys())[0]][:-1]))]
    nucleotide_content_per_base_df = pd.DataFrame({"Position": position_number_list,
                                                   "A": nucleotide_content_per_base['A'][:-1],
                                                   "C": nucleotide_content_per_base['C'][:-1],
                                                   "G": nucleotide_content_per_base['G'][:-1],
                                                   "T": nucleotide_content_per_base['T'][:-1]})
    return nucleotide_content_per_base_df


# Misha
def make_N_content_df(nucleotide_content_per_base):
    """
    Make N content pandas DataFrame
    :param nucleotide_content_per_base: dic, Nucleotide: Nucleotide content per position in read (%)
    :return:
            N_content_per_base_df: pandas DataFrame, Position in read (bp), N content per position (%)
    """
    position_number_list =\
        [i + 1 for i in range(len(nucleotide_content_per_base[list(nucleotide_content_per_base.keys())[0]][:-1]))]
    N_content_per_base_df = pd.DataFrame({"Position": position_number_list,
                                          "N": nucleotide_content_per_base['N'][:-1]})
    return N_content_per_base_df


# Misha
def make_GC_content_df(gc_dic_rounded, gc_theor_dic_rounded):
    """
    Make GC content pandas DataFrame
    :param gc_dic_rounded: dic, GC content (%):	count in Experiment
    :param gc_theor_dic_rounded: dic, GC content (%): count in Theoretical Distribution
    :return:
            GC_content_df: pandas DataFrame, GC content (%), count in Experiment, count in Theoretical Distribution
    """
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


# Misha
def make_nucleotide_content_plot(nucleotide_content_per_base_df, output_dir):
    """
    Plot Sequence content across all bases
    :param nucleotide_content_per_base_df: pandas DataFrame, Position in read (bp), Nucleotide content per position (%)
    :param output_dir: str, path to output directory
    :return: plot
    """
    sns.set_style('whitegrid')
    fig, ax = plt.subplots(figsize=(14, 10))
    # fig.set_size_inches(11.7, 8.27)
    ax.set(xlim=(0, len(nucleotide_content_per_base_df) + 1), ylim=(0, 100))
    palette = ['green', 'blue', 'black', 'red']
    ax.xaxis.set_major_locator(mticker.MultipleLocator(8))
    ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
    plt.grid(axis="x")
    ax.yaxis.set_major_locator(plt.MaxNLocator(10))
    nucleotide_content_per_base_df_melted = nucleotide_content_per_base_df.melt('Position', var_name='Nucleotide',
                                                                                value_name='vals')
    nucleotide_plot = sns.lineplot(x="Position", y="vals", hue='Nucleotide',
                                   data=nucleotide_content_per_base_df_melted, ax=ax, palette=palette,
                                   linewidth=1, dashes=False)

    ax.set_title(label='Sequence content across all bases', fontsize=14)
    ax.set_xlabel('Position in read (bp)', fontsize=14)
    ax.set_ylabel('Nucleotide content (%)', fontsize=14)
    step = 4
    xticks = np.arange(0, len(nucleotide_content_per_base_df) + step, step=step)
    for x0, x1 in zip(xticks[::2], xticks[1::2]):
        plt.axvspan(x0 + step / 2, x1 + step / 2, color='black', alpha=0.1, zorder=0)
    plt.xticks(xticks)
    leg = plt.legend(labels=["% " + i for i in nucleotide_content_per_base_df.keys()[1:]], handlelength=0,
                     loc="upper right", borderaxespad=0, prop={'size': 8})
    for line, text in zip(leg.get_lines(), leg.get_texts()):
        text.set_color(line.get_color())
    plt.savefig(Path(output_dir, 'nucleotide_content.png'), format='png')
    #  plt.show()


# Misha
def make_N_content_plot(N_content_per_base_df, output_dir):
    """
    Plot N content across all bases
    :param N_content_per_base_df: pandas DataFrame, Position in read (bp), N content per position (%)
    :param output_dir: str, path to output directory
    :return: plot
    """
    sns.set_style('whitegrid')

    fig, ax = plt.subplots(figsize=(14, 10))
    # fig.set_size_inches(11.7, 8.27)
    ax.set(xlim=(0, len(N_content_per_base_df) + 1), ylim=(0, 100))
    palette = ['red']
    ax.xaxis.set_major_locator(mticker.MultipleLocator(8))
    ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
    plt.grid(axis="x")
    ax.yaxis.set_major_locator(plt.MaxNLocator(10))
    N_content_per_base_df_melted = N_content_per_base_df.melt('Position', var_name='Nucleotide', value_name='vals')
    N_plot = sns.lineplot(x="Position", y="vals", hue='Nucleotide',
                          data=N_content_per_base_df_melted, ax=ax, palette=palette,
                          linewidth=1, dashes=False)
    ax.set_title(label='N content across all bases', fontsize=14)
    ax.set_xlabel('Position in read (bp)', fontsize=14)
    ax.set_ylabel('Nucleotide content (%)', fontsize=14)
    step = 4
    xticks = np.arange(0, len(N_content_per_base_df) + 1, step=step)
    for x0, x1 in zip(xticks[::2], xticks[1::2]):
        plt.axvspan(x0 + step / 2, x1 + step / 2, color='black', alpha=0.1, zorder=0)
    plt.xticks(xticks)
    leg = plt.legend(labels=["% N"], handlelength=0, loc="upper right", borderaxespad=0, prop={'size': 8})
    for line, text in zip(leg.get_lines(), leg.get_texts()):
        text.set_color(line.get_color())
    plt.savefig(Path(output_dir, 'N_content.png'), format='png')
    #  plt.show()


# Misha
def make_GC_content_plot(GC_content_df, output_dir):
    """
    Plot GC distribution over all sequences
    :param GC_content_df: pandas DataFrame, GC content (%),	count in Experiment, count in Theoretical Distribution
    :param output_dir: str, path to output directory
    :return: plot
    """
    sns.set_style('whitegrid')
    fig, ax = plt.subplots(figsize=(14, 10))
    # fig.set_size_inches(11.7, 8.27)
    ax.set(xlim=(1, 100))
    palette = ['red', 'blue']
    ax.xaxis.set_major_locator(mticker.MultipleLocator(4))
    ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
    plt.grid(axis="x")
    GC_content_df_melted = GC_content_df.melt('GC', var_name='GC content', value_name='vals')
    GC_content_plot = sns.lineplot(x="GC", y="vals", hue='GC content',
                                   data=GC_content_df_melted, ax=ax, palette=palette,
                                   linewidth=1, dashes=False)
    ax.set_title(label='GC distribution over all sequences', fontsize=14)
    ax.set_xlabel('Mean GC content (%)', fontsize=14)
    ax.set_ylabel('Number of reads', fontsize=14)
    step = 4
    xticks = np.arange(0, 100, step=step)
    for x0, x1 in zip(xticks[::2], xticks[1::2]):
        plt.axvspan(x0 + step / 2, x1 + step / 2, color='black', alpha=0.1, zorder=0)
    plt.xticks(xticks)
    leg = plt.legend(handlelength=0, loc="upper right", borderaxespad=0, prop={'size': 8})
    for line, text in zip(leg.get_lines(), leg.get_texts()):
        text.set_color(line.get_color())
    plt.savefig(Path(output_dir, 'GC_content.png'), format='png')
    #  plt.show()


# Misha
def make_nucleotide_and_gc_plots(nucleotide_content_per_base, gc_dic_rounded, gc_theor_dic_rounded, output_dir):
    """
    Calls and executes the functions that draw plots: Sequence content across all bases, N content across all bases,
    GC distribution over all sequences and generates the pandas DataFrames necessary for this
    :param nucleotide_content_per_base: dic, Nucleotide: Nucleotide content per position in read (%)
    :param gc_dic_rounded: dic, GC content (%):	count in Experiment
    :param gc_theor_dic_rounded: dic, GC content (%): count in Theoretical Distribution
    :param output_dir: str, path to output directory
    :return: plots
    """
    nucleotide_content_per_base_df = make_nucleotide_content_df(nucleotide_content_per_base)

    N_content_per_base_df = make_N_content_df(nucleotide_content_per_base)

    GC_content_df = make_GC_content_df(gc_dic_rounded, gc_theor_dic_rounded)

    make_nucleotide_content_plot(nucleotide_content_per_base_df, output_dir)

    make_N_content_plot(N_content_per_base_df, output_dir)

    make_GC_content_plot(GC_content_df, output_dir)


# Anton
def plot_duplications_and_distinct_duplications(frequency_procent, dist_frequency_procent, perc_if_dupl, output_dir):
    """
    Plots Sequence Duplication Levels. Take arguments from duplicates function from calculations module
    :param frequency_procent: list of values,each value is a procent of sequences with duplication
    level less than value in x defined below
    :param dist_frequency_procent:
    :param frequency_procent: list of values,each value is a procent of distinct sequences with duplication
    level less than value in x defined velow
    :param perc_if_dupl: total unique reads/total reads
    :param output_dir: directory where picture will be stored
    :return: None, generates Sequence Duplication Levels
    """

    x = ["1", "2", "3", "4", "5", "6", "7", "8", "9", ">10", ">50", ">100", ">500", ">1k", ">5k", ">10k"]
    y = frequency_procent[1:]
    z = dist_frequency_procent[1:]

    plt.figure(figsize=(14, 10))
    plt.tick_params(bottom=False, left=False)

    plt.plot(x, z, color="red", label="some")
    plt.plot(x, y, color="blue")

    plt.title(f"Percent of seqs remaining if duplicated {perc_if_dupl}%", fontsize=14)
    plt.xlabel("Sequence duplication level", fontsize=14)

    plt.ylim(0, 1)
    plt.yticks(ticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1],
               labels=[0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    plt.xticks(ticks=x)

    for i in range(0, 16, 2):
        plt.axvspan(i + 0.5, i + 1.5, color="#e6e6e6", zorder=0)
    plt.grid(axis="x")
    leg = plt.legend(labels=["% Total sequences", "% Deduplicated sequences"], handlelength=0, loc="upper right", borderaxespad=0, prop={'size': 14})
    for line, text in zip(leg.get_lines(), leg.get_texts()):
        text.set_color(line.get_color())

    plt.savefig(Path(output_dir, 'sequence_duplication_levels.png'))
    return None


# Anton
def plot_table_of_overrepresent(read_sequence, read_count, read_percentage, read_source, output_dir):
    """
    Plots a table of a overrepresented sequences. Takes parameters from overrepresented_sequences
    function from calculations module.
    :param read_sequence: sequence of the read that has presentation more than 0.1% of total reads
    :param read_count: count of this sequence
    :param read_percentage: percentage of count
    :param read_source: source of sequence. Default: No hit
    :param output_dir: directory where picture will be stored
    :return: None, generates table of overrepresented sequences
    """
    if len(read_sequence) > 0:
        fig = go.Figure(data=[go.Table(
            header=dict(values=["Sequence", "Count", "Percentage", "Possible Source"],
                        fill_color="#000080",
                        font=dict(color='white', family="Times New Roman", size=16)),
            cells=dict(values=[read_sequence, read_count, read_percentage, read_source],
                       align="left",
                       font=dict(color='black', family="Times New Roman", size=12)),
            columnwidth=[0.3, 0.04, 0.09, 0.08])
        ])
        fig.update_layout(width=890, height=235 + 20 * len(read_sequence))
        fig.write_image(Path(output_dir, 'overrepresented.png'))
        return None
    plt.figure(figsize=(4, 1))
    plt.text(-0.1, 0.93, "Overrepresented sequences", color="#830006", fontsize=15)
    plt.text(-0.15, 0.75, "No overrepresented sequences", fontsize=8)
    plt.axis('off')
    plt.savefig(Path(output_dir, 'overrepresented.png'), dpi=200)
    return None


# Anna
def plot_adapters_content(adapters_count, plot_edge, output_dir):
    """
    Plot Adapter content
    :param adapters_count: dic, adapter_name:count each position
    :return: plot for adapter count
    """
    max_x = plot_edge
    x = np.arange(1, max_x + 1)
    _, ax = plt.subplots(figsize=(14, 10))
    for key in adapters_count.keys():
        plt.plot(x, adapters_count[key][:max_x], label=key)
    plt.xlim(1, max_x + 1)
    plt.ylim(0, 100)
    plt.grid(axis="x")
    step = 5
    xticks = np.arange(1, max_x + 1, step=step)
    for x0, x1 in zip(xticks[::2], xticks[1::2]):
        plt.axvspan(x0 + step / 2, x1 + step / 2, color='black', alpha=0.1, zorder=0)
    plt.xticks(xticks)
    plt.yticks(np.arange(0, 101, step=10))
    leg = plt.legend(handlelength=0, loc="upper right", borderaxespad=0, prop={'size': 14})
    for line, text in zip(leg.get_lines(), leg.get_texts()):
        text.set_color(line.get_color())
    plt.title("% Adapter", fontsize=14)
    plt.xlabel("Position in read (bp)", fontsize=14)
    plt.savefig(Path(output_dir, 'adapter_content.png'))
    # plt.show()
