import numpy as np
import statistics


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


def make_nucleotide_and_gc_calculations(sequences_list, gc_content_list):
    nucleotide_content_per_base = get_nucleotide_content(get_nucleotide_per_position(sequences_list))

    GC_theoretical_distribution = get_GC_theoretical_distribution(gc_content_list)

    gc_theor_dic_rounded = make_GC_content_integer_rounding(GC_theoretical_distribution)

    gc_dic_rounded = make_GC_content_integer_rounding(gc_content_list)

    gc_theor_dic_rounded, gc_dic_rounded = make_GC_content_and_GC_theoretical_have_equal_length(gc_theor_dic_rounded,
                                                                                                gc_dic_rounded)
    return nucleotide_content_per_base, gc_dic_rounded, gc_theor_dic_rounded
