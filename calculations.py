import subprocess
import os

def count_lines_in_file(path_to_file):
    """Counts number of lines in text file

     NOTE! This function will be deleted after Anyas code deployment!!!"""

    cmd = f"wc -l {path_to_file}" + " | awk '{print $1}'"
    cmd_result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    return int(cmd_result.stdout.decode())

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
                    trim.write(line[0:51]+'\n')
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
            for i in range(len(countings_list)-1,-1,-1):
                result_1.write(str(countings_list[i][1]) + "\n")

    os.remove("trimmed.txt")
    return None


def duplicates(total_reads):
    """
    Calculates duplications, distinct duplications and percent if duplicated (material
    for Sequence Duplication Levels graph)
    NOTE!  This function may be used only after preprocess function as it uses counted_reads.txt file
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


def overrepresented_sequences(total_reads):
    """
    Search for overrepresented sequences (sequences with count above 0,1% of total reads in fastq
    file given in preprocess)
    NOTE!  This function may be used only after preprocess function as it uses
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
    return (read_sequence, read_count, read_percentage,read_source)
