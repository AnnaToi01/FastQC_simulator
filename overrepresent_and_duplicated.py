import matplotlib.pyplot as plt
import subprocess
import plotly.graph_objects as go


def preprocess(filename, mode="fast"):
    with open(f"{filename}") as reads, open("for_trim.txt", mode="w") as trim:
        count = 0
        for line in reads:
            count += 1
            if count % 4 == 2:
                trim.write(line)

    with open("for_trim.txt") as trim, open("trimmed.txt", mode="w") as result:
        for line in trim:
            if len(line) > 75:
                result.write(line[0:51] + '\n')
            else:
                result.write(line)

    cmd = "rm for_trim.txt"
    subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)

    if mode == "basic":

        cmd_1 = "cat trimmed.txt | sort | uniq -c | sort -k1 -r -n\
        | tee counted_reads_reverse.txt | sort -k1 -n | awk '{print $1}' > counted_reads.txt"
        subprocess.run(cmd_1, shell=True, stdout=subprocess.PIPE)

    elif mode == "fast":

        with open("trimmed.txt") as base, open("counted_reads.txt", mode="w") as result_1, open(
                "counted_reads_reverse.txt",
                mode="w") as result_2:
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
            countings_list.sort(key=lambda x: x[1])
            for i in countings_list:
                result_1.write(str(i[1]) + "\n")

    cmd_2 = "rm trimmed.txt"
    subprocess.run(cmd_2, shell=True, stdout=subprocess.PIPE)


def count_lines_in_file(filename):
    cmd = f"wc -l {filename}" + " | awk '{print $1}'"
    cmd_result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    return int(cmd_result.stdout.decode())


def overrepresented_sequences(filename):
    total_reads = count_lines_in_file(filename) // 4
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
    if len(read_sequence) > 0:
        fig = go.Figure(data=[go.Table(
            header=dict(values=["Sequence", "Count", "Percentage", "Possible Source"],
                        fill_color="blue",
                        font=dict(color='white', family="Times New Roman", size=16)),
            cells=dict(values=[read_sequence, read_count, read_percentage, read_source],
                       align="left",
                       font=dict(color='black', family="Times New Roman", size=12)),
            columnwidth=[0.3, 0.04, 0.09, 0.08])
        ])
        fig.update_layout(width=890, height=1000)
        fig.write_image("overrepresented.png")
        return
    plt.figure(figsize=(4, 1))
    plt.text(-0.1, 0.93, "Overrepresented sequences", color="#830006", fontsize=15)
    plt.text(-0.15, 0.75, "No overrepresented sequences", fontsize=8)
    plt.axis('off')
    plt.savefig("overrepresented.png", dpi=400)
    return


def dublicates(filename):
    total_reads = count_lines_in_file(filename) // 4
    total_lines = count_lines_in_file("counted_reads.txt")
    perc_if_dupl = round((total_lines / total_reads) * 100, 2)
    duplication_level = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 50, 100, 500, 1000, 5000, 10000, 100000]
    frequency = [0] * 18
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
                while duplication_level[i - 1] < d_level <= duplication_level[i]:
                    count += d_level
                    if read_count < total_lines:
                        line = reads.readline()
                        read_count += 1
                        d_level = int(line)
                    else:
                        break
                frequency[i] += count
            else:
                frequency[i + 1] += d_level

    frequency = frequency[0:10] + frequency[11:]
    frequency_procent = [i / total_reads for i in frequency]

    dist_duplication_level = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 50, 100, 500, 1000, 5000, 10000, 100000]
    dist_frequency = [0] * 18
    with open("counted_reads.txt") as reads:
        read_count = 0
        for i in range(1, 17):
            if read_count >= total_lines:
                break
            line = reads.readline()
            read_count += 1
            d_level = int(line)
            if dist_duplication_level[i - 1] < d_level <= dist_duplication_level[i]:
                count = 0
                while dist_duplication_level[i - 1] < d_level <= dist_duplication_level[i]:
                    count += 1
                    if read_count < total_lines:
                        line = reads.readline()
                        read_count += 1
                        d_level = int(line)
                    else:
                        break
                dist_frequency[i] += count
            else:
                dist_frequency[i + 1] += 1

    dist_frequency = dist_frequency[0:10] + dist_frequency[11:]
    dist_frequency_procent = [i / total_lines for i in dist_frequency]

    x = ["1", "2", "3", "4", "5", "6", "7", "8", "9", ">10", ">50", ">100", ">500", ">1k", ">5k", ">10k"]
    y = frequency_procent[1:]
    z = dist_frequency_procent[1:]

    plt.figure(figsize=(8, 6))
    plt.tick_params(bottom=False, left=False)

    plt.plot(x, z, color="red", label="some")
    plt.plot(x, y, color="blue")

    plt.title(f"Percent of seqs remaining if duplicated {perc_if_dupl}%", fontsize=9)
    plt.xlabel("Sequence duplication level", fontsize=7)

    plt.ylim(0, 1)
    plt.yticks(ticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1],
               labels=[0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100], fontsize=7)
    plt.xticks(ticks=x, fontsize=8)

    for i in range(0, 16, 2):
        plt.axvspan(i + 0.5, i + 1.5, color="#e6e6e6")
    plt.grid(axis="y")

    textstr = "% Deduplicated sequences\n "
    props = dict(boxstyle='square', facecolor='white', edgecolor='grey', alpha=1)

    plt.text(12.542, 0.992, textstr, fontsize=7,
             verticalalignment='top', bbox=props, color="red")
    plt.text(12.542, 0.948, "% Total sequences", color="blue", fontsize=7)

    plt.savefig("doublicates.png", dpi=300)


def overseq_and_doublicates(filename, mode="fast"):
    preprocess(filename, mode)
    overrepresented_sequences(filename)
    dublicates(filename)
    cmd = "rm counted_reads.txt counted_reads_reverse.txt"
    subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)

#if __name__ == "main":
overseq_and_doublicates("reads.fastq")
