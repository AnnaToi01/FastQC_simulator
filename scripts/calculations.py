import numpy as np
import pandas as pd
import statistics
from itertools import islice
from collections import Counter
import os
from pathlib import Path
import subprocess
import logging


# Ivan
def SetLogger(logger_name):
    """
    Create custom logger and set its configuration
    :param logger_name: name of created logger
    """
    logger = logging.getLogger('argparse')  # Create a custom logger
    logger.setLevel(logging.INFO)
    c_handler = logging.StreamHandler()  # Create handlers
    c_handler.setLevel(logging.INFO)
    c_format = logging.Formatter('%(levelname)s: %(message)s')  # Create formatters and add them to handlers
    c_handler.setFormatter(c_format)
    logger.addHandler(c_handler)  # Add handlers to the logger


# Anna
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


# Anna
# Misha
def calculate_read_length(path_to_file):
    """
    Calculate read length distribution
    :param path_to_file: str, path to fastq file
    :return:
            length_dic: dic, read_length:count
            max_length: int, maximum read length
            min_length: int, minimum read length
            total_reads: int, total count of reads
            gc_content_list: list, GC content per read (%)
            mean_gc_content: float, average GC content
            nucleotide_per_position_dic: dic, Position in read (bp): Nucleotide in this position per read
    """
    length_dic = {}

    gc_content_list = []
    nucleotide_per_position_dic = {}

    with open(path_to_file, 'r') as fastq_file:
        total_reads = 0
        for line_num, sequence in enumerate(fastq_file):
            if line_num % 4 == 1:
                total_reads += 1
                length = len(sequence.strip())
                length_dic[length] = length_dic.get(length, 0) + 1

                gc_content_list.append((sequence.lower().count('g') + sequence.lower().count('c')) * 100 /
                                       len(sequence))

                for position_number in range(len(sequence)):
                    if position_number in nucleotide_per_position_dic:
                        nucleotide_per_position_dic[position_number].append(sequence[position_number])
                    else:
                        nucleotide_per_position_dic[position_number] = [sequence[position_number]]

        mean_gc_content = round(statistics.mean(gc_content_list), 2)

    max_length = max(length_dic, key=int)
    min_length = min(length_dic, key=int)
    return (length_dic, max_length, min_length,
            total_reads, gc_content_list,
            mean_gc_content, nucleotide_per_position_dic)


# Ivan
def calc_quality_scores_per_pos(path_to_file, max_length):
    """
    Create dictionary of dictionaries , where keys - nucleotide positions in reads.
    Values - dictionaries: keys - unique quality scores in ascii format, values - their counts.
    1) Function create tmp file from fastq where only each fourth line with qualities is contained.
    2) Then this tmp file is read by chunks and matrix 10000*n (n - max read length) is created.
    3) For each matrix column unique ascii symbols and their counts are calculated.
    :path_to_file: str, path to fastq file
    :param max_length: int, max read length in fastq file
    :return: 1) df_quality - data frame with quality counts per position
             2) phred - int, phred quality
             3) ascii_symb_min - str, ascii symbol corresponding to minimal quality
    """
    file_name = Path(path_to_file).name.split('.')[0] + '_quality_tmp.txt'
    path_to_tmp_file = Path(path_to_file).parent / file_name
    subprocess.run("awk 'NR%4==0' {} > {}".format(path_to_file, path_to_tmp_file), shell=True)
    quality_counts_per_pos = [{}] * max_length
    with open(path_to_tmp_file, 'r') as quality_file:
        while True:
            chunk = list(islice(quality_file, 10000))
            if not chunk:
                break
            chunk_processed = pd.DataFrame([list(seq.strip()) for seq in chunk])
            quality_counts_per_pos_chunk = list(chunk_processed.apply(lambda x: dict(Counter(x)), axis=0))
            quality_counts_per_pos = [dict(Counter(x) + Counter(y)) for x, y in zip(quality_counts_per_pos,
                                                                                    quality_counts_per_pos_chunk)]
    os.remove(path_to_tmp_file)
    quality_dict = {str(i + 1): dct for i, dct in enumerate(quality_counts_per_pos)}
    df_quality, phred, ascii_symb_min = convert_quality_per_pos_dict_to_df(quality_dict=quality_dict)
    if max_length >= 110:
        df_quality = compress_quality_per_pos_df(df_quality=df_quality)
    return df_quality, phred, ascii_symb_min


# Ivan
def convert_quality_per_pos_dict_to_df(quality_dict):
    """
    Convert dict with unique quality scores per position and their counts to data frame
    :param quality_dict: dict with quality values
    :return: 1) df_quality - data frame with quality counts per position
             2) phred - int, phred quality
             3) ascii_symb_min - str, ascii symbol corresponding to minimal quality
    """
    df_quality = pd.DataFrame(quality_dict)
    df_quality = df_quality.reset_index().rename(columns={'index': 'ascii_symb'})
    is_not_none = [x is not None for x in df_quality.ascii_symb]
    df_quality = df_quality[is_not_none]
    min_score = df_quality.ascii_symb.apply(lambda x: ord(x)).min()
    ascii_symb_min = chr(min_score)
    if min_score >= 64:
        phred = 64
    else:
        phred = 33
    df_quality['quality'] = df_quality.ascii_symb.apply(lambda x: ord(x) - phred)
    df_quality = df_quality.fillna(0)
    return df_quality, phred, ascii_symb_min


# Ivan
def compress_quality_per_pos_df(df_quality):
    """
    Reduce data frame with quality scores for each position in read by collapsing neighbouring
    columns into one. Apply for fastq file with large reads.
    :param df_quality: data frame with quality values and their counts for each position in read
    :return df_quality_red - collapsed data frame
    """
    df_quality_red = df_quality[['ascii_symb', 'quality']].copy()
    positions = [col_name for col_name in df_quality.columns if col_name not in ['ascii_symb', 'quality']]
    if len(positions) % 2 == 0:
        for i in range(1, len(positions) + 1, 2):
            cols = [str(i), str(i + 1)]
            df_quality_red[f'{i}-{i + 1}'] = df_quality[cols].apply(lambda x: x[0] + x[1], axis=1)
    else:
        for i in range(1, len(positions) - 2, 2):
            cols = [str(i), str(i + 1)]
            df_quality_red[f'{i}-{i + 1}'] = df_quality[cols].apply(lambda x: x[0] + x[1], axis=1)
        cols = [str(i + 2), str(i + 3), str(i + 4)]
        df_quality_red[f'{i + 2}-{i + 4}'] = df_quality[cols].apply(lambda x: x[0] + x[1] + x[2], axis=1)
    return df_quality_red


# Ivan
def calc_stats_for_boxplot(df):
    """
    Calculate necessary statistics for drawing boxplots
    :param df: data frame with quality scores per position and their counts
    :return: 1) stats_for_boxplots - dictionary with statistics of quality distribution per position
             2) means - np.array with quality means per position
    """
    positions = [col_name for col_name in df.columns if col_name not in ['ascii_symb', 'quality']]
    stats_for_boxplots = []
    means = np.zeros(len(positions))
    for i, pos in enumerate(positions):
        # obtain spreaded list
        df_subset = df[[pos, 'quality']]
        df_subset = df_subset[df_subset[pos] != 0]
        df_subset[pos] = df_subset[pos].astype(int)
        quality_list = [[q_score] * count for q_score, count in zip(df_subset['quality'], df_subset[pos])]
        quality_list = sum(quality_list, [])
        # calculate mean quality
        mean_quality = np.mean(quality_list)
        means[i] = mean_quality
        # calculate statistics for boxplots
        stats = np.quantile(quality_list, [0.25, 0.5, 0.75, 0.1, 0.9])
        names = ['q1', 'med', 'q3', 'whislo', 'whishi']
        stats_dict = {name: st for name, st in zip(names, stats)}
        stats_for_boxplots.append(stats_dict)
    return stats_for_boxplots, means


# Ivan
def calc_quality_per_read(ascii_seq, max_length, phred):
    """
    Caclulate quality for each position in a read
    :param ascii_seq: str, fourth line from fastq file with quality values
    :param max_length: int, max read length in fastq file
    :phred: int, type of phred quality (33 or 64)
    :return: array of shape (1, max_length) filled with phred quality scores
    """
    quality_array = np.zeros(max_length)
    for i, symbol in enumerate(ascii_seq):
        quality_array[i] = ord(symbol) - phred
    return quality_array


# Ivan
def calc_quality_scores_per_read(path_to_file, total_reads, max_length, phred):
    """
    Calculate mean quality scores per read and their counts
    :param path_to_file: str, path to fastq file
    :param total_reads: int, total number of reads in fastq file
    :param max_length: int, max read length in fastq file
    :phred: int, type of phred quality (33 or 64)
    :return: read_quality_dict - dict, where keys - quality scores per read, values - their counts
    """
    mean_quality_per_read = np.zeros(total_reads)
    with open(path_to_file, "r") as fastq_file:
        for i, line in enumerate(fastq_file):
            if i % 4 == 3:
                ascii_seq = line.strip()
                read_length = len(ascii_seq)
                quality_values = calc_quality_per_read(ascii_seq=ascii_seq, max_length=max_length, phred=phred)
                mean_read_quality = np.mean(quality_values[:read_length])
                mean_quality_per_read[i // 4] = mean_read_quality
    max_quality = max(mean_quality_per_read)
    thres = max(np.floor(mean_quality_per_read))
    mean_quality_per_read = [round(val) if val < thres else val for val in mean_quality_per_read]
    mean_quality_per_read = [np.floor(val) if thres < val < max_quality else val for val in mean_quality_per_read]
    read_quality_dict = Counter(mean_quality_per_read)
    return read_quality_dict


# Misha
def get_nucleotide_content(nucleotide_per_position_dic):
    """
    Make Nucleotide content per position in read (%) dictionary
    :param nucleotide_per_position_dic: dic, Position in read (bp): Nucleotide in this position per read
    :return:
            nucleotide_content_per_base: dic, Nucleotide: Nucleotide content per position in read (%)
    """
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


# Misha
def get_GC_theoretical_distribution(gc_content_list):
    """
    Make GC content Theoretical Distribution based on the sequencing data
    :param gc_content_list: list, GC content per read (%)
    :return:
            GC_theoretical_distribution: list, GC content per read (%) obtained from Theoretical Distribution
    """
    GC_mode = statistics.mode(gc_content_list)
    GC_sd = statistics.stdev(gc_content_list)
    number_of_sequences = len(gc_content_list)
    GC_theoretical_distribution = np.random.normal(GC_mode, GC_sd, number_of_sequences)
    return GC_theoretical_distribution


# Misha
def make_GC_content_integer_rounding(gc_content_list):
    """
    Make GC content integer rounding for further plotting
    :param gc_content_list: list, GC content per read (%)
    :return:
            gc_dic_rounded: dic, GC content (%): count
    """
    gc_dic_rounded = {}
    for i in range(len(gc_content_list)):
        rounded = round(gc_content_list[i])
        if rounded in gc_dic_rounded:
            gc_dic_rounded[rounded] += 1
        else:
            gc_dic_rounded[rounded] = 1
    return gc_dic_rounded


# Misha
def make_GC_content_and_GC_theoretical_have_equal_length(gc_theor_dic_rounded, gc_dic_rounded):
    """
    Make GC content in Experiment and Theoretical Distribution dictionaries have the same GC content range of values
    and the same length for further merging in pandas DataFrame
    :param gc_dic_rounded: dic, GC content (%):	count in Experiment
    :param gc_theor_dic_rounded: dic, GC content (%): count in Theoretical Distribution
    :return:
            gc_dic_rounded: dic, GC content (%): count in Experiment
            gc_theor_dic_rounded: dic, GC content (%): count in Theoretical Distribution
    """
    for gc in gc_dic_rounded:
        if gc not in gc_theor_dic_rounded:
            gc_theor_dic_rounded[gc] = 0

    for gc in gc_theor_dic_rounded:
        if gc not in gc_dic_rounded:
            gc_dic_rounded[gc] = 0
    return gc_theor_dic_rounded, gc_dic_rounded


# Misha
def make_nucleotide_and_gc_calculations(gc_content_list, nucleotide_per_position_dic):
    """
    Calls and executes the functions that perform data preprocessing for further drawing following plots:
    Sequence content across all bases, N content across all bases, GC distribution over all sequences
    :param gc_content_list, list, GC content per read (%)
    :param nucleotide_per_position_dic: dic, Position in read (bp): Nucleotide in this position per read
    :return:
            nucleotide_content_per_base: dic, Nucleotide: Nucleotide content per position in read (%)
            gc_dic_rounded: dic, GC content (%): count in Experiment
            gc_theor_dic_rounded: dic, GC content (%): count in Theoretical Distribution
    """
    nucleotide_content_per_base = get_nucleotide_content(nucleotide_per_position_dic)

    GC_theoretical_distribution = get_GC_theoretical_distribution(gc_content_list)

    gc_theor_dic_rounded = make_GC_content_integer_rounding(GC_theoretical_distribution)

    gc_dic_rounded = make_GC_content_integer_rounding(gc_content_list)

    gc_theor_dic_rounded, gc_dic_rounded = make_GC_content_and_GC_theoretical_have_equal_length(gc_theor_dic_rounded,
                                                                                                gc_dic_rounded)
    return nucleotide_content_per_base, gc_dic_rounded, gc_theor_dic_rounded


# Anton
def preprocess(path_to_file, mode="fast"):
    """
    Preproceess isolate sequences from fast file, trim sequences with length >75 bp to 50bp,
    count reads, then sort counted reads to 2 files: counted_reads.txt and counted_reads_reverse.txt
    :param path_to_file: path to input file
    :param mode: basic mode uses command line utilities for core calculations, fast mode uses
    :return: returns None, genereates 2 files: counted_reads.txt and counted_reads_reverse.txt
    """

    with open(f"{path_to_file}") as reads, open("trimmed.txt", mode="w") as trim:
        count = 0
        for line in reads:
            count += 1
            if count % 4 == 2:
                if len(line) > 75:
                    trim.write(line[0:51] + '\n')
                else:
                    trim.write(line)

    if mode == "basic":

        cmd_1 = "cat trimmed.txt | sort | uniq -c | sort -k1 -r -n\
        | tee counted_reads_reverse.txt | sort -k1 -n | awk '{print $1}' > counted_reads.txt"
        subprocess.run(cmd_1, shell=True, stdout=subprocess.PIPE)

    elif mode == "fast":

        with open("trimmed.txt") as base, open("counted_reads.txt", mode="w") as result_1, open(
                "counted_reads_reverse.txt", mode="w") as result_2:
            countings = {}
            for line in base:
                if line in countings:
                    countings[line] += 1
                else:
                    countings[line] = 1
            countings_list = list(countings.items())
            countings_list.sort(key=lambda x: x[1], reverse=True)
            for i in countings_list:
                result_2.write(str(i[1]) + " ")
                result_2.write(i[0])
            for i in range(len(countings_list) - 1, -1, -1):
                result_1.write(str(countings_list[i][1]) + "\n")

    os.remove("trimmed.txt")
    return


# Anton
def duplicates(total_reads):
    """
    Calculates duplications, distinct duplications and percent if duplicated (material
    for Sequence Duplication Levels graph)
    NOTE! This function may be used only after preprocess function as it uses counted_reads.txt file
    :param total_reads: total number of sequences in fastq file taken by preprocess
    :return: frequency_procent, dist_frequency_procent, perc_if_dupl
    """

    cmd = "wc -l counted_reads.txt | awk '{print $1}'"
    cmd_result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    total_lines = int(cmd_result.stdout.decode())
    perc_if_dupl = round((total_lines / total_reads) * 100, 2)

    duplication_level = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 50, 100, 500, 1000, 5000, 10000, 100000]
    frequency = [0] * 18
    dist_frequency = [0] * 18
    with open("counted_reads.txt") as reads:
        read_count = 0
        for i in range(1, 17):
            if read_count >= total_lines:
                break
            line = reads.readline()
            read_count += 1
            d_level = int(line)
            if duplication_level[i - 1] < d_level <= duplication_level[i]:
                count = 0
                dist_count = 0
                while duplication_level[i - 1] < d_level <= duplication_level[i]:
                    count += d_level
                    dist_count += 1
                    if read_count < total_lines:
                        line = reads.readline()
                        read_count += 1
                        d_level = int(line)
                    else:
                        break
                frequency[i] += count
                dist_frequency[i] += dist_count
            else:
                frequency[i + 1] += d_level
                dist_frequency[i + 1] += 1

    frequency = frequency[0:10] + frequency[11:]
    frequency_procent = [i / total_reads for i in frequency]
    dist_frequency = dist_frequency[0:10] + dist_frequency[11:]
    dist_frequency_procent = [i / total_lines for i in dist_frequency]
    os.remove("counted_reads.txt")
    return (frequency_procent, dist_frequency_procent, perc_if_dupl)


# Anton
def overrepresented_sequences(total_reads):
    """
    Search for overrepresented sequences (sequences with count above 0,1% of total reads in fastq
    file given in preprocess)
    NOTE! This function may be used only after preprocess function as it uses
    counted_reads_reverse.txt file
    :param total_reads: total_reads: total number of sequences in fastq file taken by preprocess
    :return: read_sequence, read_count, read_percentage,read_source
    """
    threshold = round(total_reads * 0.001)
    read_sequence = []
    read_count = []
    read_percentage = []
    with open("counted_reads_reverse.txt") as count:
        for line in count:
            line = line.split()
            if int(line[0]) > threshold:
                read_sequence.append(line[1][:50])
                read_count.append(line[0])
                read_percentage.append(int(line[0]) / total_reads * 100)
            else:
                break
    read_source = ["No hit"] * len(read_sequence)
    os.remove("counted_reads_reverse.txt")
    return (read_sequence, read_count, read_percentage, read_source)


# Anna
def get_adapters(path_to_file, max_length):
    """
    Get adapter names  and their sequences
    :param path_to_file: str, path to file with adapter
    :param max_length: int, maximum read length
    :return: adapters: dic, adapter_name:sequence
             adapters_count: dic, adapter_name:count each position(0)
    """
    adapters = {}
    adapters_count = {}

    with open(path_to_file, "r") as f:
        for line in f:
            k, v = line.strip().split("\t")
            adapters[k] = v
            adapters_count[k] = np.zeros(max_length)

    return adapters, adapters_count


def adapter_content(path_to_file, adapters, adapters_count, max_length):
    """
    Gets the cumulative adapter frequency for each adapter and position in read
    :param path_to_file: str, path to fastq
    :param adapters: dic, adapter_name:sequence
    :param adapters_count: dic, adapter_name:count each position(0)
    :param max_length: maximum read length
    :return: adapters_count,
             max(plot_edge) + 2 - int, xlim
    """
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
    if all(key == 0 for key in plot_edge):
        plot_edge = max_length
    else:
        plot_edge = max(plot_edge)
    return adapters_count, plot_edge


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
        data=[
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


def tile_sequence(path_to_file, max_length, lowest_char):
    """
    Generates dictionary for plotting of tile sequence quality
    :param path_to_file: str, path to file
    :param max_length: int, maximum read length
    :param lowest_char: str, character signifying lowest character quality
    :return:
            tile_quality: dic, tile_number: array(average quality per position)
    """
    tile_quality = {}
    tile_count = {}

    def quality_per_position(qual_sequence, max_length):
        """
        Returns array with the ASCII order of symbols in a sequence
        :param qual_sequence: str, sequence of letters
        :param max_length: int, maximum read length
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
                count_array = np.array([1] * (len_qual) + [0] * (max_length - len_qual))
                # For each tile number generate array with counts of positions included in counting
                tile_count[tile_num] = tile_count.get(tile_num, np.zeros(max_length)) + count_array
                qual_array = quality_per_position(qual_sequence, max_length)
                # For each tile number generate array with sum of qualities per position
                tile_quality[tile_num] = tile_quality.get(tile_num, np.zeros(max_length)) + qual_array

    # Ignore some errors
    # np.seterr(invalid='ignore')
    # Normalizing the quality per count
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
