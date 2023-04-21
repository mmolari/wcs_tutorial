import argparse
import matplotlib.pyplot as plt
from Bio import Phylo


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tree",
        help="Path to the tree file",
        type=str,
    )
    parser.add_argument(
        "--clade",
        help="list of strains in the clade",
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "--pair",
        help="selected pair of isolates",
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "--fig",
        help="output figure",
        type=str,
    )
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


if __name__ == "__main__":

    args = parse_args()

    # load tree
    tree = Phylo.read(args.tree, "newick")

    # prune tree
    tree = prune_tree(tree, args.clade)

    # get isolates
    i1, i2 = args.pair

    # plot the tree, highlighting the pair
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    Phylo.draw(
        tree,
        axes=ax,
        label_func=lambda x: x.name if x.name in args.clade else "",
        do_show=False,
        label_colors={i1: "C0", i2: "C1"},
    )
    ax.set_yticks([])
    ax.set_ylabel("")
    for k in ["top", "right", "left"]:
        ax.spines[k].set_visible(False)
    plt.tight_layout()
    plt.savefig(args.fig, dpi=150, facecolor="white")
    plt.close(fig)
