from pathlib import Path
import numpy as np
import pandas as pd


class QualityError(Exception):
    """
    Exception is raised when the quality score is undefined.
    """

    def __init__(self, lowest_char, message):
        self.lowest_char = lowest_char
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f'{self.message}  (Yours was "{self.lowest_char}" with value "{ord(self.lowest_char)}"'


def calculate_read_length(path_to_file):
    '''
    Calculate read length distribution
    :param path_to_file: str, path to fastq file
    :return:
            length_dic: dic, read_length:count
            max_length: int, maximum read length
            min_length: int, minimum read length
            total_reads: int, total count of reads
    '''
    length_dic = {}
    with open(path_to_file, 'r') as fastq_file:
        total_reads = 0
        for line_num, sequence in enumerate(fastq_file):
            if line_num % 4 == 1:
                total_reads += 1
                length = len(sequence.strip())
                length_dic[length] = length_dic.get(length, 0) + 1
    max_length = max(length_dic, key=int)
    min_length = min(length_dic, key=int)
    return length_dic, max_length, min_length, total_reads


def get_adapters(path_to_file, max_length):
    '''
        Get adapter names  and their sequences
    :param path_to_file: str, path to file with adapter
    :param max_length: int, maximum read length
    :return: adapters: dic, adapter_name:sequence
             adapters_count: dic, adapter_name:count each position(0)
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
    :return: adapters_count: dic, adapter_name:count for each position
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


def basic_statistics(path_to_file, lowest_char, mean_gc_content, min_length, max_length, total_reads):
    """
    Generates the r
    :param path_to_file: str, path to fastq
    :param lowest_char: str, character with lowest ASCII order in quality
    :param mean_gc_content: int, mean GC content as percentage
    :param min_length: int, minimal read length
    :param max_length: int, maximal read length
    :param total_reads: int, total count of reads
    :return:
            basicstats: pd Dataframe, basic statistics graph
            offset: int, offset for conversion of ASCII to Q score and vice versa
            encoding: str, Illumina type

    """
    file_name = Path(path_to_file).name
    sanger_offset = 33
    illumina_1_3_offset = 64

    file_type = "Conventional base calls"

    ord_lowest_char = ord(lowest_char)

    if ord_lowest_char < 33:
        raise QualityError(lowest_char, message="No known encodings with chars < 33")
    elif ord_lowest_char < 64:
        encoding = "Sanger / Illumina 1.9"
        offset = sanger_offset
    elif ord_lowest_char == illumina_1_3_offset + 1:
        encoding = "Illumina 1.3"
        offset = illumina_1_3_offset
    elif ord_lowest_char <= 126:
        encoding = "Illumina 1.5"
        offset = illumina_1_3_offset
    else:
        raise QualityError(lowest_char, message="No known encodings with chars > 126")

    if max_length == min_length:
        sequence_length = str(max_length)
    else:
        sequence_length = "-".join([str(min_length), str(max_length)])

    poor_quality = 0

    basicstats = pd.DataFrame(
        data = [
            ["Filename", file_name],
            ["File type", file_type],
            ["Encoding", encoding],
            ["Total Sequences", total_reads],
            ["Sequences flagged as poor quality", poor_quality],
            ["Sequence length", sequence_length],
            ["%GC", mean_gc_content]
        ],
        columns=["Measure", "Value"]
    )

    return basicstats, offset, encoding


def tile_sequence(path_to_file, max_length):
    """
    Generates dictionary for plotting of tile sequence quality
    :param path_to_file: path to file
    :param max_length: maximum read length
    :return:
            tile_quality: dic, tile_number: array(average quality per position)
    """
    tile_quality = {}
    tile_count = {}

    def quality_per_position(qual_sequence, max_length):
        """
        Returns array with the ASCII order of symbols in a sequence
        :param qual_sequence: sequence of letters
        :param max_length: maximum read length
        :return:
                qual_array: np.array, array of ASCII order of symbol in each position
        """
        qual_array = np.zeros(max_length)
        for i, symbol in enumerate(qual_sequence):
            qual_array[i] = ord(symbol)
        return qual_array

    with open(path_to_file, "r") as fastq_file:
        for line_num, qual_sequence in enumerate(fastq_file):
            if line_num % 4 == 0:
                try:
                    tile_num = qual_sequence.split(":")[4]  # get tile number in sequence identity
                except IndexError:
                    return  # Sometimes no tile numbers included, tile sequence not generated then

            elif line_num % 4 == 3:
                qual_sequence = qual_sequence.strip()  # get quality sequence
                len_qual = len(qual_sequence)  # length of quality sequence
                # For each tile number generate array with counts of positions included in counting
                tile_count[tile_num] = tile_count.get(tile_num, np.zeros(max_length)) + \
                                       np.array([1] * (len_qual) + [0] * (max_length - len_qual))
                # For each tile number generate array with sum of qualities per position
                tile_quality[tile_num] = tile_quality.get(tile_num, np.zeros(max_length)) + \
                                         quality_per_position(qual_sequence, max_length)
    # Ignore some errors
    # np.seterr(invalid='ignore')
    # Normalizing the quality per count
    lowest_char = "#"
    min_num = ord(lowest_char)
    for tile_num in tile_count.keys():
        tile_count_array = np.where(tile_count[tile_num] == 0, 1, tile_count[tile_num])
        tile_quality_array = np.where(tile_quality[tile_num] == 0, min_num, tile_quality[tile_num])

        tile_quality[tile_num] = np.divide(tile_quality_array, tile_count_array)

    # Normalizing quality per position in whole tile - create array where quality per position is summed up
    av_qual_position = np.zeros(max_length)
    for tile_num in tile_quality.keys():
        for i in range(max_length):
            av_qual_position[i] += tile_quality[tile_num][i]
    # Normalize the generated array by the number of tiles -> average quality per position
    for i in range(max_length):
        av_qual_position[i] /= len(tile_quality.keys())
    # Comparing the generated tile_quality array to the means per position by subtracting
    for tile_num in tile_quality.keys():
        for i in range(max_length):
            tile_quality[tile_num][i] -= av_qual_position[i]

    return tile_quality