import calculations as calc
import plots as plot


def read_file_and_generate_dictionary_with_fastq_data(input_fastq):
    with open(input_fastq) as r:
        fastq_file = r.readlines()

        sequences_list = []

        for i in range(1, len(fastq_file), 4):
            sequences_list.append(fastq_file[i][:-1])

    return sequences_list


if __name__ == "__main__":
    input_fastq = '/mnt/d/IB2021_2022/command_line/Project_1/raw_data/amp_res_1.fastq'
    # input_fastq = '/mnt/d/IB2021_2022/command_line/Project_2/SRR1705851.fastq'
    sequences_list = read_file_and_generate_dictionary_with_fastq_data(input_fastq)
    gc_content_list, GC_mean = calc.calculate_gc(sequences_list)
    plot.make_nucleotide_and_gc_plots(*calc.make_nucleotide_and_gc_calculations(sequences_list, gc_content_list))
