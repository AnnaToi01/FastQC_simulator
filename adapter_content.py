import numpy as np
import matplotlib.pyplot as plt
from length_distribution import calculate_read_length


def get_adapters(path_to_file, max_length):
    '''
        Get adapter names  and their sequences
    :param path_to_file: str, path to file with adapter
    :param max_length: int, maximum read length
    :return: adapters - dic, adapter_name:sequence
             adapters_count - dic, adapter_name:count each position(0)
    '''

    adapters = {}
    adapters_count = {}

    with open("adapters.txt", "r") as f:
        for line in f:
            k, v = line.strip().split("\t")
            adapters[k] = v
            adapters_count[k] = np.zeros(max_length)

    return adapters, adapters_count


def adapter_content(path_to_file, adapters, adapters_count):
    '''
    Gets the cumulative adapter frequency for each adapter and position in read
    :param path_to_file: str, path to fastq
    :param adapters: dic, adapter_name:sequence
    :param adapters_count: dic, adapter_name:count each position(0)
    :return: adapters_count,
             max(plot_edge) + 2 - int, xlim
    '''
    plot_edge = []
    with open(path_to_file, "r") as fastq_file:
        for key in adapters.keys():
            n = 0
            for line_num, sequence in enumerate(fastq_file):
                if line_num % 4 == 1:
                    n += 1
                    ind = sequence.strip().find(adapters[key])
                    if ind != -1:
                        adapters_count[key][ind] += 1
            plot_edge.append(np.argmax(adapters_count[key]))
            adapters_count[key] = np.cumsum(adapters_count[key]) / n * 100
            fastq_file.seek(0)
    return adapters_count, max(plot_edge)


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


if __name__ == "__main__":
    path_to_file = "amp_res_1P_30.fastq"
    _, max_length = calculate_read_length(path_to_file)
    adapters, adapter_count = get_adapters(path_to_file, max_length)
    adapter_count, plot_edge = adapter_content(path_to_file, adapters, adapter_count)
    plot_adapters_content(adapter_count, plot_edge)



