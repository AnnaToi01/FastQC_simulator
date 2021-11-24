import numpy as np
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt


def draw_boxplot_with_qscores(df, phred, encoding, stats_for_boxplots, means, output_dir):
    '''
    Draw boxplots with quality scores per position
    :param df: data frame with unique quality scores per position and their counts
    :param phred: int, type of phred quality (33 or 64)
    :param encoding: type of Illumina
    :param stats_for_boxplots: dict with statistics of quality distribution per position
    :param means: numpy.array with quality means per position
    :param output_dir: directory where picture will be stored
    '''
    positions = [col_name for col_name in df.columns if col_name not in ['ascii_symb', 'quality']]
    y_max = max([stats['whishi'] for stats in stats_for_boxplots])
    _, ax = plt.subplots(figsize=(25, 15))
    bplot = ax.bxp(stats_for_boxplots, showfliers=False, patch_artist=True, widths=0.7, zorder=2)
    for box, med in zip(bplot['boxes'], bplot['medians']):
        box.set_facecolor('yellow')
        box.set_linewidth(1)
        med.set_color('red')
        med.set_linewidth(2)
    ax.set_xlabel('Position in read (bp)', fontsize=18)
    ax.set_title(f'Quality scores across all bases ({encoding})', fontsize=20)
    ax.xaxis.grid(True, zorder=1)
    length = len(positions)
    if length <= 150:
        xlabels = np.array(positions[::2])
        xticks = np.arange(1, length+1, 2)
    else:
        xlabels = np.array(positions[::3])
        xticks = np.arange(1, length+1, 3)
    ylabels = yticks = np.arange(0, 42, 2)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, rotation=45)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_xlim([0, length+1])
    ax.set_ylim([0, y_max+0.5])
    ax.plot(np.arange(1, len(means)+1), means, color='blue', zorder=3)
    ax.axhspan(28, 42, facecolor='green', alpha=0.2)
    ax.axhspan(20, 28, facecolor='orange', alpha=0.2)
    ax.axhspan(0, 20, facecolor='red', alpha=0.2)
    plt.savefig(Path(output_dir, 'per_base_quality.png'))


def draw_per_sequence_quality_scores(read_quality_dict, output_dir):
    '''
    Draw Quality score distribution over all sequences
    :param read_quality_dict: dict, where keys - quality scores per read, values - their counts
    :param output_dir: directory where picture will be stored
    '''
    q_scores_unique = list(read_quality_dict.keys())
    counts = list(read_quality_dict.values())
    min_quality = min(q_scores_unique)
    max_quality = max(q_scores_unique)
    plt.figure(figsize=(14, 10))
    sns.lineplot(q_scores_unique, counts, color='red', label='Average Quality per read')
    plt.grid(axis="y")
    plt.xlabel('Mean Sequence Quality (Phred Score)', fontsize=14)
    plt.title('Quality score distribution over all sequences', fontsize=14)
    plt.xlim([min_quality-1, max_quality+2])
    plt.ylim([0, None])
    xticks = np.arange(min_quality, max_quality+2, 2)
    plt.xticks(xticks)
    for x in xticks[::2]:
        plt.axvspan(x + 1, x + 3, color='black', alpha=0.1, zorder=1)
    leg = plt.legend(loc="upper right", prop={'size': 14}, handlelength=0, borderaxespad=0)
    for line, text in zip(leg.get_lines(), leg.get_texts()):
        text.set_color(line.get_color())
    plt.savefig(Path(output_dir, 'per_sequence_quality.png'))
