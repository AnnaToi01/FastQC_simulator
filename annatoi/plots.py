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
    # plt.show()


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
    # plt.show()


def plot_tile_sequence(tile_quality, max_length):
    """
    Plots per tile sequence quality
    :param tile_quality: dic, tile_number: array(average quality per position)
    :return: None, saves the generated figure and shows it
    """

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


    fig = plt.figure(figsize=(20, 20))
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
    plt.savefig("tile_sequence.png")

