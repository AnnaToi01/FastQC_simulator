import calculations as calc
import plots as plot


if __name__ == "__main__":
    # input_fastq = '/mnt/d/IB2021_2022/command_line/Project_1/raw_data/amp_res_1.fastq'
    # input_fastq = '/mnt/d/IB2021_2022/command_line/Project_2/SRR1705851.fastq'
    path_to_file = '/mnt/d/IB2021_2022/command_line/Project_1/raw_data/amp_res_1.fastq'
    output_dir = '/mnt/d/IB2021_2022/python/FastQC_simulator/results/'

    length_dic, max_length, min_length, total_reads, gc_content_list, GC_mean, nucleotide_per_position_dic =\
        calc.calculate_read_length(path_to_file)

    plot.make_nucleotide_and_gc_plots(*calc.make_nucleotide_and_gc_calculations(gc_content_list,
                                                                                nucleotide_per_position_dic),
                                      output_dir)
