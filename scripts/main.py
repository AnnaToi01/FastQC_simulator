import calculations as calc
import plotting as plot
import argparse
from pathlib import Path
import sys
import warnings
import logging
import time
import datetime
from PIL import Image

warnings.simplefilter('ignore')

if __name__ == '__main__':
    # create argparse logger
    calc.SetLogger(logger_name='argparse')
    logger = logging.getLogger('argparse')
    # argparse main parameters
    parser = argparse.ArgumentParser(description='FastQC simulator', epilog='Enjoy the program! :)',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', type=str, help='Path to fastq file', required=True, metavar='')
    parser.add_argument('-o', '--output', type=str, help='Path to output folder for storing results',
                        required=True, metavar='')
    parser.add_argument('-a', '--adapters', type=str, help='Path to file with adapters. Default: ./adapters.txt',
                        default='./adapters.txt', metavar='')
    args = parser.parse_args()
    if not Path(args.input).exists():
        logger.warning('File {} does not exist!'.format(args.input))
        logger.info('Abort calculations. Specify correct path to file.')
        sys.exit()
    if not Path(args.adapters).exists():
        logger.warning('File {} does not exist!'.format(args.adapters))
        logger.info('Abort calculations. Specify correct path to file.')
        sys.exit()
    logger.info('Start fastq quality control...')
    start_time = time.process_time()
    Path(args.output).mkdir(parents=True, exist_ok=True)
    # insert all functions for calculations and plotting here

    # Anna
    length_dic, max_length, min_length, total_reads, gc_content_list, mean_gc_content, nucleotide_per_position_dic = \
        calc.calculate_read_length(path_to_file=args.input)

    # Ivan
    df_quality, phred, ascii_symb_min = calc.calc_quality_scores_per_pos(path_to_file=args.input, max_length=max_length)
    stats_for_boxplots, means = calc.calc_stats_for_boxplot(df=df_quality)
    read_quality_dict = calc.calc_quality_scores_per_read(path_to_file=args.input, total_reads=total_reads,
                                                          max_length=max_length, phred=phred)

    # Anton
    calc.preprocess(path_to_file=args.input)
    read_sequence, read_count, read_percentage, read_source = calc.overrepresented_sequences(
        total_reads=total_reads)  # see comment in calculations
    frequency_procent, dist_frequency_percent, perc_if_dupl = calc.duplicates(
        total_reads=total_reads)  # see comment in calculations

    # Anna
    adapters, adapters_count = calc.get_adapters(path_to_file=args.adapters, max_length=max_length)
    adapters_count, plot_edge = calc.adapter_content(path_to_file=args.input, adapters=adapters,
                                                     adapters_count=adapters_count, max_length=max_length)
    tile_quality = calc.tile_sequence(path_to_file=args.input, max_length=max_length, lowest_char=ascii_symb_min)
    basicstats, offset, encoding = calc.basic_statistics(path_to_file=args.input, lowest_char=ascii_symb_min,
                                                         mean_gc_content=mean_gc_content, min_length=min_length,
                                                         max_length=max_length, total_reads=total_reads)

    # Plotting

    # Anna
    plot.plot_basic_statistics(basicstats=basicstats, output_dir=args.output)

    # Ivan
    plot.draw_boxplot_with_qscores(df=df_quality, phred=phred, encoding=encoding, stats_for_boxplots=stats_for_boxplots,
                                   means=means, output_dir=args.output)
    # Anna
    plot.plot_tile_sequence(tile_quality=tile_quality, max_length=max_length, output_dir=args.output)
    # Ivan
    plot.draw_per_sequence_quality_scores(read_quality_dict=read_quality_dict, output_dir=args.output)
    # Misha
    plot.make_nucleotide_and_gc_plots(*calc.make_nucleotide_and_gc_calculations(gc_content_list,
                                                                                nucleotide_per_position_dic),
                                      output_dir=args.output)

    # Anton
    plot.plot_table_of_overrepresent(read_sequence=read_sequence, read_count=read_count,
                                     read_percentage=read_percentage, read_source=read_source, output_dir=args.output)
    plot.plot_duplications_and_distinct_duplications(frequency_procent=frequency_procent,
                                                     dist_frequency_procent=dist_frequency_percent,
                                                     perc_if_dupl=perc_if_dupl, output_dir=args.output)

    # Anna
    plot.plot_length_distribution(length_distribution_dic=length_dic, output_dir=args.output)
    plot.plot_adapters_content(adapters_count=adapters_count, plot_edge=plot_edge, output_dir=args.output)

    images = ["basic_statistics.png", "per_base_quality.png", "per_tile_sequence_quality.png",
              "per_sequence_quality.png", "nucleotide_content.png", "GC_content.png",
              "N_content.png", "sequence_length_distribution.png", "sequence_duplication_levels.png",
              "overrepresented.png", "adapter_content.png"
              ]
    im_list_out = []
    for image in images:
        im = Image.open(Path(args.output, image))
        im_list_out.append(im.convert("RGB"))

    im_list_out[0].save(str(Path(args.output, "result.pdf")), save_all=True, append_images=im_list_out[1:])
    end_time = time.process_time()
    total_time = str(datetime.timedelta(seconds=end_time - start_time))
    logger.info('End fastq quality_control')
    logger.info(f'Total time for FastQC analysis: {total_time}')
