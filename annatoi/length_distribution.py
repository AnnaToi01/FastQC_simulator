import numpy as np
import matplotlib.pyplot as plt
import math


def calculate_read_length(path_to_file):
    '''
    Calculate read length distribution
    :param path_to_file: str, path to fastq file
    :return: dictionary, read_length:count
    '''
    length_dic = {}
    with open(path_to_file, 'r') as fastq_file:
        for line_num, sequence in enumerate(fastq_file):
            if line_num % 4 == 1:
                length = len(sequence.strip())
                length_dic[length] = length_dic.get(length, 0) + 1
    global max_length
    max_length = max(length_dic, key=int)
    return length_dic, max_length


def plot_length_distribution(length_distribution_dic):
    '''
    Plot Sequence length distribution
    :param length_distribution_dic: dictionary, read_length:count
    :return: plot
    '''
    x = [int(key) for key in length_distribution_dic.keys()]
    min_x = min(x)
    max_x = max(x)
    x.append(min_x - 1)
    x.append(max_x + 1)
    min_x -= 1
    max_x += 1
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
    diffference_yticks = max_y // 10
    digits_zero = 1
    while diffference_yticks // 10 > 1:
        diffference_yticks = diffference_yticks // 10
        digits_zero *= 10
    step_y = math.floor(max_y / diffference_yticks / digits_zero) * digits_zero
    plt.yticks(np.arange(0, math.ceil(max(y)/1000)*1000+step_y, step=step_y))
    plt.ticklabel_format(style="plain", useOffset=False)
    leg = plt.legend(handlelength=0, loc="upper right", borderaxespad=0, prop={'size': 6.8})
    for line, text in zip(leg.get_lines(), leg.get_texts()):
        text.set_color(line.get_color())
    plt.title("Distribution of sequence lengths over all sequences")
    plt.xlabel("Sequence length")
    plt.savefig("sequence_length_distribution.png")
    plt.show()


if __name__ == "__main__":
    path_to_file = "amp_res_1P_30.fastq"
    length_distribution_dic = calculate_read_length(path_to_file)
    plot_length_distribution(length_distribution_dic)
