import matplotlib.pyplot as plt
import plotly.graph_objects as go


def plot_table_of_overrepresent(read_sequence, read_count, read_percentage,read_source):

    """Takes arguments from overrepresented_sequences function (calculations module). Saves png picture with
     table of the overrepresented sequences"""

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
        fig.update_layout(width=890, height=235+20*len(read_sequence))
        fig.write_image("overrepresented.png")
        return
    plt.figure(figsize=(4, 1))
    plt.text(-0.1, 0.93, "Overrepresented sequences", color="#830006", fontsize=15)
    plt.text(-0.15, 0.75, "No overrepresented sequences", fontsize=8)
    plt.axis('off')
    plt.savefig("overrepresented.png", dpi=400)
    return

def plot_duplications_and_distinct_duplications(frequency_procent, dist_frequency_procent, perc_if_dupl):
    """Takes arguments from duplicates function (calculations module). Saves png picture with duplications and
     distinct duplicatins"""

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
    return