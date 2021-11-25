import numpy as np
import statistics


def calculate_read_length(path_to_file):
    '''
    Calculate read length distribution
    :param path_to_file: str, path to fastq file
    :return:
            length_dic: dic, read_length:count
            max_length: int, maximum read length
            min_length: int, minimum read length
            total_reads: int, total count of reads
            gc_content_list: list, GC content per read (%)
            GC_mean: float, average GC content
            nucleotide_per_position_dic: dic, Position in read (bp): Nucleotide in this position per read 
    '''
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

        GC_mean = round(statistics.mean(gc_content_list), 2)

    max_length = max(length_dic, key=int)
    min_length = min(length_dic, key=int)
    return length_dic, max_length, min_length, total_reads, gc_content_list, GC_mean, nucleotide_per_position_dic


def get_nucleotide_content(nucleotide_per_position_dic):
    '''
    Make Nucleotide content per position in read (%) dictionary
    :param nucleotide_per_position_dic: dic, Position in read (bp): Nucleotide in this position per read
    :return:
            nucleotide_content_per_base: dic, Nucleotide: Nucleotide content per position in read (%)
    '''
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


def get_GC_theoretical_distribution(gc_content_list):
    '''
    Make GC content Theoretical Distribution based on the sequencing data
    :param gc_content_list: list, GC content per read (%)
    :return:
            GC_theoretical_distribution: list, GC content per read (%) obtained from Theoretical Distribution
    '''
    GC_mode = statistics.mode(gc_content_list)
    GC_sd = statistics.stdev(gc_content_list)
    number_of_sequences = len(gc_content_list)
    GC_theoretical_distribution = np.random.normal(GC_mode, GC_sd, number_of_sequences)
    return GC_theoretical_distribution


def make_GC_content_integer_rounding(gc_content_list):
    '''
    Make GC content integer rounding for further plotting
    :param gc_content_list: list, GC content per read (%)
    :return:
            gc_dic_rounded: dic, GC content (%): count
    '''
    gc_dic_rounded = {}
    for i in range(len(gc_content_list)):
        rounded = round(gc_content_list[i])
        if rounded in gc_dic_rounded:
            gc_dic_rounded[rounded] += 1
        else:
            gc_dic_rounded[rounded] = 1
    return gc_dic_rounded


def make_GC_content_and_GC_theoretical_have_equal_length(gc_theor_dic_rounded, gc_dic_rounded):
    '''
    Make GC content in Experiment and Theoretical Distribution dictionaries have the same GC content range of values
    and the same length for further merging in pandas DataFrame
    :param gc_dic_rounded: dic, GC content (%):	count in Experiment
    :param gc_theor_dic_rounded: dic, GC content (%): count in Theoretical Distribution
    :return:
            gc_dic_rounded: dic, GC content (%): count in Experiment
            gc_theor_dic_rounded: dic, GC content (%): count in Theoretical Distribution
    '''
    for gc in gc_dic_rounded:
        if gc not in gc_theor_dic_rounded:
            gc_theor_dic_rounded[gc] = 0

    for gc in gc_theor_dic_rounded:
        if gc not in gc_dic_rounded:
            gc_dic_rounded[gc] = 0
    return gc_theor_dic_rounded, gc_dic_rounded


def make_nucleotide_and_gc_calculations(gc_content_list, nucleotide_per_position_dic):
    '''
    Calls and executes the functions that perform data preprocessing for further drawing following plots: 
    Sequence content across all bases, N content across all bases, GC distribution over all sequences
    :param gc_content_list, list, GC content per read (%)
    :param nucleotide_per_position_dic: dic, Position in read (bp): Nucleotide in this position per read 
    :return:
            nucleotide_content_per_base: dic, Nucleotide: Nucleotide content per position in read (%)
            gc_dic_rounded: dic, GC content (%): count in Experiment
            gc_theor_dic_rounded: dic, GC content (%): count in Theoretical Distribution
    '''
    nucleotide_content_per_base = get_nucleotide_content(nucleotide_per_position_dic)

    GC_theoretical_distribution = get_GC_theoretical_distribution(gc_content_list)

    gc_theor_dic_rounded = make_GC_content_integer_rounding(GC_theoretical_distribution)

    gc_dic_rounded = make_GC_content_integer_rounding(gc_content_list)

    gc_theor_dic_rounded, gc_dic_rounded = make_GC_content_and_GC_theoretical_have_equal_length(gc_theor_dic_rounded,
                                                                                                gc_dic_rounded)
    return nucleotide_content_per_base, gc_dic_rounded, gc_theor_dic_rounded
