import argparse
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as spc


from itertools import product

import pypangraph as pp


def parse_args():
    args = argparse.ArgumentParser()
    args.add_argument("--pangraph", type=str, required=True)
    args.add_argument("--anchor_block", type=str, required=True)
    args.add_argument("--leaves_colors", type=str, required=True)
    args.add_argument("--fig_paths", type=str, required=True)
    args.add_argument("--fig_matrix", type=str, required=True)
    args.add_argument("--block_colors", type=str, required=True)
    args.add_argument("--shared_len_df", type=str, required=True)
    return args.parse_args()


def shared_path_from_anchor(p1, p2, Ls, anchor):
    """given two paths, a dictionary of block lengths and an anchor block,
    return the length of shared co-path starting from the anchor block"""
    # list of block ids and strandedness
    B1, B2 = [list(p.block_ids) for p in [p1, p2]]
    S1, S2 = [list(p.block_strands) for p in [p1, p2]]

    # index of anchor block in each path
    a1, a2 = [l.index(anchor) for l in [B1, B2]]

    # flip paths if anchor block is on the reverse strand
    if not S1[a1]:
        B1 = B1[::-1]
        S1 = [not s for s in reversed(S1)]
        a1 = B1.index(anchor)
    if not S2[a2]:
        B2 = B2[::-1]
        S2 = [not s for s in reversed(S2)]
        a2 = B2.index(anchor)

    shared_L = 0

    # rev comparison, move backward from anchor block
    for k in range(0, min(a1, a2) + 1):
        i1, i2 = a1 - k, a2 - k
        if (B1[i1] == B2[i2]) and (S1[i1] == S2[i2]):
            shared_L += Ls[B1[i1]]
        else:
            break
    # fwd comparison, move forward from anchor block
    for k in range(1, min(len(B1) - a1, len(B2) - a2)):
        i1, i2 = a1 + k, a2 + k
        if (B1[i1] == B2[i2]) and (S1[i1] == S2[i2]):
            shared_L += Ls[B1[i1]]
        else:
            break

    return shared_L


def hierarchical_clustering_order(shared_L):
    """given the dataframe of shared lengths performs hierarchical clustering
    and returns the order of the paths"""
    # compute distance matrix
    S = shared_L.pivot(index="p1", columns="p2", values="shared_L")

    # dictionary of self-similarity
    diag = {c: S.loc[c, c] for c in S.columns}

    # divide element i,j by the minimum of the diagonal elements of i and j
    # and transform into a distance matrix
    D = S / S.apply(
        lambda x: np.minimum(diag[x.name], [diag[i] for i in x.index]),
        axis=1,
        result_type="broadcast",
    )
    D = 1 - D.to_numpy()

    # perform hierarchical clustering
    linkage = spc.linkage(spc.distance.squareform(D), method="complete")
    # order of paths
    order = spc.leaves_list(linkage)
    path_order = S.index[order]

    return path_order


def plot_paths(pan, paths_order, leaves_colors, core_anchor, fig_savename):
    """given a pangraph, a path order, a dictionary of leaves colors and a core
    block to be used as an anchor, plots all of the paths in a linear representation,
    assigning colors to blocks occurring multiple times."""

    # dictionary of block lengths
    df = pan.to_blockstats_df().sort_values("count", ascending=False)
    Ls = df["len"].to_dict()

    # generate colors for blocks
    def __color_generator():
        cm = mpl.colormaps.get_cmap("tab20")
        colors = [cm(i) for i in range(20)]
        colors.pop(15)  # remove light gray
        i = 0
        while True:
            yield colors[i]
            i += 1

    cg = __color_generator()
    block_colors = {
        b: next(cg) if df["count"][b] > 1 else "lightgray" for b in df.index
    }

    fig, ax = plt.subplots(1, 1, figsize=(12, 7))
    for n, p_name in enumerate(paths_order):
        p = pan.paths[p_name]
        l = 0
        cid = list(p.block_ids).index(core_anchor)
        s = p.block_strands[cid]
        bids = list(p.block_ids) if s else list(reversed(p.block_ids))
        for b in bids:
            lb = Ls[b]
            plt.plot([l, l + lb], [n, n], color=block_colors[b], lw=5)
            l += lb
        # add isolate name, with color extracted from leaves_colors
        strain = p_name.split("-")[0]
        plt.text(
            -300,
            n,
            strain,
            ha="right",
            va="center",
            fontsize=10,
            color=leaves_colors[strain],
        )
    plt.yticks([])
    for sp in ["top", "right", "left"]:
        ax.spines[sp].set_visible(False)
    plt.tight_layout()
    plt.savefig(fig_savename, dpi=300)
    plt.close(fig)

    return block_colors


def plot_shared_length_matrix(df, path_order, leaves_colors, fig_savename):
    """Plots the matrix of pairwise shared path length"""

    # create matrix with entries in selected order
    M = df.pivot(index="p1", columns="p2", values="shared_L")
    M = M.loc[path_order, path_order].to_numpy() / 1000

    N = M.shape[0]

    # plot the matrix
    fig, ax = plt.subplots(1, 1, figsize=(7, 4))
    g = ax.matshow(M)
    ax.set_xticks(range(N))
    ax.set_xticklabels(["" for _ in range(N)])
    ax.set_yticks(range(N))
    ax.set_yticklabels(["" for _ in range(N)])

    for path_n, path_name in enumerate(path_order):
        isolate = path_name.split("-")[0]
        ax.text(
            -1,
            path_n,
            isolate,
            ha="right",
            va="center",
            fontsize=10,
            color=leaves_colors[isolate],
        )

    for sp in ["top", "right", "left", "bottom"]:
        ax.spines[sp].set_visible(False)
    plt.colorbar(g, ax=ax, fraction=0.046, pad=0.04, label="shared length (kbp)")
    plt.tight_layout()
    plt.savefig(fig_savename, dpi=300)
    plt.close(fig)


if __name__ == "__main__":

    args = parse_args()

    # load pangraph
    pan = pp.Pangraph.load_json(args.pangraph)

    # block stats dataframe
    bdf = pan.to_blockstats_df()

    # load anchor block
    anchor = args.anchor_block
    assert bdf["core"][anchor], "the anchor block must be core"

    # diciotnary of block lengths
    Ls = bdf["len"].to_dict()

    # evaluate pairwise shared length from anchor block for all path pairs
    shared_L_df = []
    for p1, p2 in product(pan.paths, repeat=2):
        shared_L = shared_path_from_anchor(p1, p2, Ls, anchor)
        shared_L_df.append([p1.name, p2.name, shared_L])
    shared_L_df = pd.DataFrame(shared_L_df, columns=["p1", "p2", "shared_L"])

    # save shared length dataframe
    shared_L_df.to_csv(args.shared_len_df, index=False)

    # perform hierarchical clustering and find optial path order
    path_order = hierarchical_clustering_order(shared_L_df)

    # load leaves colors
    leaves_colors = pd.read_csv(args.leaves_colors, index_col=0)["color"].to_dict()

    # plot paths in linear representation
    block_colors = plot_paths(pan, path_order, leaves_colors, anchor, args.fig_paths)

    # plot shared length matrix
    plot_shared_length_matrix(shared_L_df, path_order, leaves_colors, args.fig_matrix)

    # save block colors
    hex_colors = {k: mpl.colors.to_hex(v) for k, v in block_colors.items()}
    hex_colors = pd.Series(hex_colors, name="Colour")
    hex_colors.to_csv(args.block_colors, index_label="Name")
