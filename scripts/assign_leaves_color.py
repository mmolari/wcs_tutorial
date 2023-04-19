import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from Bio import Phylo


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", type=str, required=True)
    parser.add_argument("--paf", type=str, required=True)
    parser.add_argument("--color_csv", type=str, required=True)
    parser.add_argument("--fig", type=str, required=True)
    return parser.parse_args()


def isolates_from_paf(paf_file):
    """return all isolates containing a match from paf file"""
    paf = pd.read_csv(paf_file, sep="\t", header=None)
    return paf[0].unique()


def color_generator(N):
    cmap = plt.get_cmap("Dark2")
    for x in np.linspace(0, 1, N):
        yield cmap(x)


def plot_tree(tree, selected, colors, fig_name):
    """plot tree coloring only selected leaves"""

    fig, ax = plt.subplots(1, 1, figsize=(5, 12))
    Phylo.draw(
        tree,
        label_func=lambda l: l.name if l.name in selected else "",
        label_colors=colors,
        axes=ax,
        do_show=False,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(fig_name, facecolor="white", dpi=150)
    plt.close(fig)


if __name__ == "__main__":

    args = parse_args()

    # load newick tree
    tree = Phylo.read(args.tree, "newick")

    # get list of isolates containing the bla gene
    selected = isolates_from_paf(args.paf)

    # generate dictionary of colors for leaves
    cg = color_generator(len(selected))
    colors = {l.name: next(cg) for l in tree.get_terminals() if l.name in selected}

    # plot the tree, with selected leaves in corresponding colors
    plot_tree(tree, selected, colors, args.fig)

    # save colors to csv file
    selected_colors = {l: mpl.colors.to_hex(v) for l, v in colors.items()}
    cdf = pd.Series(selected_colors, name="color")
    cdf.to_csv(args.color_csv, index_label="isolate")
