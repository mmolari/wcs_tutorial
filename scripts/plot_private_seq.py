import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Bio import Phylo


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", type=str, help="core-genome tree")
    parser.add_argument("--dist_df", type=str, help="pairwise distance dataframe")
    parser.add_argument("--fig_scatter", type=str, help="output scatterplot")
    parser.add_argument("--fig_matrix", type=str, help="output distance matrix")
    return parser.parse_args()


def prune_tree(tree, keep_leaves):
    """Prune tree to keep only selected leaves."""

    to_prune = []
    for leaf in tree.get_terminals():
        if leaf.name not in keep_leaves:
            to_prune.append(leaf)

    for leaf in to_prune:
        tree.prune(leaf)

    return tree


def add_tree_distances(tree, df):
    """add tree distances to dataframe"""

    leaves_dict = {leaf.name: leaf for leaf in tree.get_terminals()}

    df["tree_d"] = df.apply(
        lambda x: tree.distance(leaves_dict[x["p1"]], leaves_dict[x["p2"]]),
        axis=1,
        result_type="expand",
    )

    return df


def scatterplot(df, fig_savename):
    """Produce a scatterplot of private sequence vs core genome tree distance"""

    # only select one occurrence per pair
    sdf = df[df["p1"] > df["p2"]]

    fig, ax = plt.subplots(1, 1, figsize=(4, 3))
    ax.scatter(sdf["tree_d"], sdf["private_seq"] / 1000, alpha=0.4)
    ax.set_xlabel("core genome tree distance")
    ax.set_ylabel("private sequence distance (kbp)")
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    plt.tight_layout()
    plt.savefig(fig_savename, facecolor="white", dpi=150)
    plt.close(fig)


def matrixplot(tree, M, fig_savename):
    """Plot the private sequence distance matrix next to the phylogenetic tree"""

    # list of isolates in the order in which they appear in the tree
    iso_order = [leaf.name for leaf in tree.get_terminals()]
    N = len(iso_order)

    # set the same order as the tree
    M = M.loc[iso_order, iso_order]

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))

    ax = axs[0]
    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        label_func=lambda x: x.name if x in tree.get_terminals() else "",
    )
    for k in ["top", "right", "left"]:
        ax.spines[k].set_visible(False)
    ax.set_ylabel("")
    ax.set_yticks([])

    ax = axs[1]
    g = ax.matshow(M / 1000)
    ax.set_xticks(range(N))
    ax.set_xticklabels(iso_order, rotation=90)
    ax.set_yticks(range(N))
    ax.set_yticklabels([])
    plt.colorbar(g, ax=ax, label="private seq. (kbp)", fraction=0.046, pad=0.04)

    plt.tight_layout()
    fig.subplots_adjust(wspace=-0.2, top=0.65)
    plt.savefig(fig_savename, facecolor="white", dpi=150)
    plt.close(fig)


if __name__ == "__main__":

    args = parse_args()

    # load pairwise distance dataframe
    dist_df = pd.read_csv(args.dist_df)

    # list of selected isolates
    isolates = dist_df["p1"].unique()

    # load tree and prune all non-selected sequences
    tree = Phylo.read(args.tree, "newick")
    tree = prune_tree(tree, isolates)

    # add tree distances to dataframe
    dist_df = add_tree_distances(tree, dist_df)

    # produce a scatter-plot of private sequence vs tree distance
    scatterplot(dist_df, args.fig_scatter)

    # compute pairwise distance matrix
    M = dist_df.pivot(index="p1", columns="p2", values="private_seq")

    # plot distance matrix vs tree
    matrixplot(tree, M, args.fig_matrix)
