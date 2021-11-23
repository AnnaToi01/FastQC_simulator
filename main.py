import calculations as calc
import plots as plot

filename = "filename.fastq"

calc.preprocess(filename)

read_sequence, read_count, read_percentage,read_source = calc.overrepresented_sequences(filename) # see comment in calculations

plot.plot_table_of_overrepresent(read_sequence, read_count, read_percentage,read_source)

frequency_procent, dist_frequency_procent, perc_if_dupl = calc.duplicates(filename) # see comment in calculations

plot.plot_duplications_and_distinct_duplications(frequency_procent, dist_frequency_procent, perc_if_dupl)

