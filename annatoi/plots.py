import matplotlib.pyplot as plt
import math
import numpy as np
import dataframe_image as dfi


def plot_basic_statistics(basicstats):
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

    dfi.export(basicstats_styled, 'basicstats.png')


def plot_length_distribution(length_distribution_dic):
    '''
    Plot Sequence length distribution
    :param length_distribution_dic: dictionary, read_length:count
    :return: plot
    '''
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
    plt.plot(x, y, label="Sequence length", color='r')
    plt.grid(axis="y")
    if max_x - min_x < 20:
        step = 1
    else:
        step = math.ceil((max_x - min_x) / 10)
    xticks = np.arange(min_x, max_x+step, step=step)
    for x0, x1 in zip(xticks[::2], xticks[1::2]):
        plt.axvspan(x0 + step / 2, x1 + step / 2, color='black', alpha=0.1, zorder=0)
    plt.xticks(xticks)
    difference_yticks = max_y // 10
    digits_zero = 1
    while difference_yticks // 10 > 1:
        difference_yticks = difference_yticks // 10
        digits_zero *= 10
    step_y = math.floor(max_y / difference_yticks / digits_zero) * digits_zero
    plt.yticks(np.arange(0, math.ceil(max(y)/1000)*1000+step_y, step=step_y))
    plt.ticklabel_format(style="plain", useOffset=False)
    leg = plt.legend(handlelength=0, loc="upper right", borderaxespad=0, prop={'size': 6.8})
    for line, text in zip(leg.get_lines(), leg.get_texts()):
        text.set_color(line.get_color())
    plt.title("Distribution of sequence lengths over all sequences")
    plt.xlabel("Sequence length")
    plt.savefig("sequence_length_distribution.png")
    plt.show()


def plot_adapters_content(adapters_count, plot_edge):
    '''
    Plot Adapter content
    :param adapters_count: dic, adapter_name:count each position
    :return: plot for adapter count
    '''
    max_x = plot_edge + 2
    x = np.arange(1, max_x + 1)
    for key in adapters_count.keys():
        plt.plot(x, adapters_count[key][:max_x], label=key)

    plt.xlim(1, max_x + 1)
    plt.ylim(0, 100)
    plt.grid(axis="y")
    step = 5
    xticks = np.arange(1, max_x + 1, step=step)
    for x0, x1 in zip(xticks[::2], xticks[1::2]):
        plt.axvspan(x0 + step / 2, x1 + step / 2, color='black', alpha=0.1, zorder=0)
    plt.xticks(xticks)
    plt.yticks(np.arange(0, 101, step=10))
    leg = plt.legend(handlelength=0, loc="upper right", borderaxespad=0, prop={'size': 6.8})
    for line, text in zip(leg.get_lines(), leg.get_texts()):
        text.set_color(line.get_color())
    plt.title("% Adapter")
    plt.xlabel("Position in read (bp)")
    plt.savefig("adapter_content.png")
    plt.show()
