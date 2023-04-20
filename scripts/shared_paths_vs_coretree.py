import argparse
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from Bio import Phylo
from itertools import product, combinations


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--shared_len_df", type=str, required=True)
    parser.add_argument("--tree", type=str, required=True)
    parser.add_argument("--leaves_colors", type=str, required=True)
    parser.add_argument("--fig_scatter", type=str, required=True)
    parser.add_argument("--fig_tree", type=str, required=True)
    return parser.parse_args()


def tree_leaves_distances(tree):
    """given a phylogenetic tree evaluates all pairwise distances for
    the leaves and returns them in a dictionary"""
    leaves = tree.get_terminals()
    d = {}
    for (i, j) in product(leaves, repeat=2):
        label = (i.name, j.name)
        d[label] = tree.distance(i, j)
    return d


def scatterplot(df, fig_savename):
    """Draw a scatter-plot of shared path length vs tree distance"""
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.scatter(df["shared_L"] / 1000, df["tree_dist"], alpha=0.2)
    ax.set_xlabel("shared path length (kbp)")
    ax.set_ylabel("tree distance")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.tight_layout()
    fig.savefig(fig_savename, dpi=300, bbox_inches="tight")
    plt.close(fig)


def __draw_ellipse(y1, y2, color, ax):
    """utility function to draw a half-ellipse between two points"""
    theta = np.linspace(0, np.pi, 30)
    dy = np.abs(y2 - y1) / 2
    y0 = np.mean([y1, y2])
    dx = dy / 2
    x = np.sin(theta) * dx
    y = np.cos(theta) * dy + y0
    ax.plot(x, y, color=color)


def tree_plot(tree, df, fig_savename):
    """Draw a plot to compare the distance on core-genome tree to the length of
    shared paths"""

    # set of bla-containing isolates
    leaves_idx = {l.name: n + 1 for n, l in enumerate(tree.get_terminals())}
    selected_leaves = list(leaves_colors.keys())

    fig, axs = plt.subplots(
        1, 2, figsize=(6, 10), sharey=True, gridspec_kw={"width_ratios": [2, 1]}
    )

    # draw tree, displaying the label of sleected isolates
    ax = axs[0]
    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        label_func=lambda x: x.name if x.name in selected_leaves else "",
        label_colors=leaves_colors,
    )
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.grid(True, axis="y", alpha=0.5)

    # draw links between isolates, color representing the shared path length
    ax = axs[1]

    cmap = plt.get_cmap("Blues")
    Lmin, Lmax = df["shared_L"].min() / 1000, df["shared_L"].max() / 1000
    norm = plt.Normalize(vmin=Lmin, vmax=Lmax)
    
    for _, x in df.sort_values("shared_L").iterrows():
        l1 = x["p1"].split("-")[0]
        l2 = x["p2"].split("-")[0]
        L = x["shared_L"] / 1000
        i1 = leaves_idx[l1]
        i2 = leaves_idx[l2]
        color = cmap(norm(L))
        __draw_ellipse(i1, i2, color, ax)

    ax.set_xticks([])
    for sp in ax.spines:
        ax.spines[sp].set_visible(False)

    # add colorbar
    mapp = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    plt.colorbar(mapp, ax=ax, label="shared path length (kbp)", shrink=0.5, pad=0.05)

    plt.tight_layout()
    plt.savefig(fig_savename, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":

    args = parse_args()

    # load core genome tree
    tree = Phylo.read(args.tree, "newick")

    # load shared path distance (select only one item per pair)
    df_l = pd.read_csv(args.shared_len_df)
    df_l = df_l[df_l["p1"] > df_l["p2"]]

    # load leaves colors
    leaves_colors = pd.read_csv(args.leaves_colors, index_col=0)["color"].to_dict()

    # evaluate tree distance
    tree_d = tree_leaves_distances(tree)

    # add tree distance to dataframe
    df_l["tree_dist"] = df_l.apply(
        lambda x: tree_d[(x["p1"].split("-")[0], x["p2"].split("-")[0])], axis=1
    )

    # perform plots
    scatterplot(df_l, args.fig_scatter)
    tree_plot(tree, df_l, args.fig_tree)
