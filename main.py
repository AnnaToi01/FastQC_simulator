import calculations as calc
import plots as plot

path_to_file = "filename.fastq"

calc.preprocess(path_to_file)

read_sequence, \
read_count, \
read_percentage,\
read_source = calc.overrepresented_sequences(total_reads)

plot.plot_table_of_overrepresent(read_sequence, read_count, read_percentage,read_source)

frequency_procent, \
dist_frequency_procent, \
perc_if_dupl = calc.duplicates(total_reads)

plot.plot_duplications_and_distinct_duplications(frequency_procent, dist_frequency_procent, perc_if_dupl)

