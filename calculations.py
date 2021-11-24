import numpy as np
import pandas as pd
from itertools import islice
from collections import Counter
import os
from pathlib import Path
import subprocess


def calc_quality_scores_per_pos(path_to_file, max_length):
    '''
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
    '''
    file_name = Path(path_to_file).name.split('.')[0] + '_quality_tmp.txt'
    path_to_tmp_file = Path(path_to_file).parent / file_name
    subprocess.run("awk 'NR%4==0' {} > {}".format(path_to_file, path_to_tmp_file), shell=True)
    quality_counts_per_pos = [{}]*max_length
    with open(path_to_tmp_file, 'r') as quality_file:
        while True:
            chunk = list(islice(quality_file, 10000))
            if not chunk:
                break
            chunk_processed = pd.DataFrame([list(seq.strip()) for seq in chunk])
            quality_counts_per_pos_chunk = list(chunk_processed.apply(lambda x: dict(Counter(x)), axis=0))
            quality_counts_per_pos = [dict(Counter(x)+Counter(y)) for x, y in zip(quality_counts_per_pos,
                                                                                  quality_counts_per_pos_chunk)]
    os.remove(path_to_tmp_file)
    quality_dict = {str(i+1): dct for i, dct in enumerate(quality_counts_per_pos)}
    df_quality, phred, ascii_symb_min = convert_quality_per_pos_dict_to_df(quality_dict=quality_dict)
    if max_length >= 110:
        df_quality = compress_quality_per_pos_df(df_quality=df_quality)
    return df_quality, phred, ascii_symb_min


def convert_quality_per_pos_dict_to_df(quality_dict):
    '''
    Convert dict with unique quality scores per position and their counts to data frame
    :param quality_dict: dict with quality values
    :return: 1) df_quality - data frame with quality counts per position
             2) phred - int, phred quality
             3) ascii_symb_min - str, ascii symbol corresponding to minimal quality
    '''
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
    df_quality['quality'] = df_quality.ascii_symb.apply(lambda x: ord(x)-phred)
    df_quality = df_quality.fillna(0)
    return df_quality, phred, ascii_symb_min


def compress_quality_per_pos_df(df_quality):
    '''
    Reduce data frame with quality scores for each position in read by collapsing neighbouring
    columns into one. Apply for fastq file with large reads.
    :param df_quality: data frame with quality values and their counts for each position in read
    :return df_quality_red - collapsed data frame
    '''
    df_quality_red = df_quality[['ascii_symb', 'quality']].copy()
    positions = [col_name for col_name in df_quality.columns if col_name not in ['ascii_symb', 'quality']]
    if len(positions) % 2 == 0:
        for i in range(1, len(positions)+1, 2):
            cols = [str(i), str(i+1)]
            df_quality_red[f'{i}-{i+1}'] = df_quality[cols].apply(lambda x: x[0]+x[1], axis=1)
    else:
        for i in range(1, len(positions)-2, 2):
            cols = [str(i), str(i+1)]
            df_quality_red[f'{i}-{i+1}'] = df_quality[cols].apply(lambda x: x[0]+x[1], axis=1)
        cols = [str(i+2), str(i+3), str(i+4)]
        df_quality_red[f'{i+2}-{i+4}'] = df_quality[cols].apply(lambda x: x[0]+x[1]+x[2], axis=1)
    return df_quality_red


def calc_stats_for_boxplot(df):
    '''
    Calculate necessary statistics for drawing boxplots
    :param df: data frame with quality scores per position and their counts
    :return: 1) stats_for_boxplots - dictionary with statistics of quality distribution per position
             2) means - np.array with quality means per position
    '''
    positions = [col_name for col_name in df.columns if col_name not in ['ascii_symb', 'quality']]
    stats_for_boxplots = []
    means = np.zeros(len(positions))
    for i, pos in enumerate(positions):
        # obtain spreaded list
        df_subset = df[[pos, 'quality']]
        df_subset = df_subset[df_subset[pos] != 0]
        df_subset[pos] = df_subset[pos].astype(int)
        quality_list = [[q_score]*count for q_score, count in zip(df_subset['quality'], df_subset[pos])]
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


def calc_quality_per_read(ascii_seq, max_length, phred):
    '''
    Caclulate quality for each position in a read
    :param ascii_seq: str, fourth line from fastq file with quality values
    :param max_length: int, max read length in fastq file
    :phred: int, type of phred quality (33 or 64)
    :return: array of shape (1, max_length) filled with phred quality scores
    '''
    quality_array = np.zeros(max_length)
    for i, symbol in enumerate(ascii_seq):
        quality_array[i] = ord(symbol)-phred
    return quality_array


def calc_quality_scores_per_read(path_to_file, total_reads, max_length, phred):
    '''
    Calculate mean quality scores per read and their counts
    :param path_to_file: str, path to fastq file
    :param total_reads: int, total number of reads in fastq file
    :param max_length: int, max read length in fastq file
    :phred: int, type of phred quality (33 or 64)
    :return: read_quality_dict - dict, where keys - quality scores per read, values - their counts
    '''
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
