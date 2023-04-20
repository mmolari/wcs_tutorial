import argparse
import numpy as np
import matplotlib.pyplot as plt
import pypangraph as pp


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pangraph", type=str, help="pangraph file")
    parser.add_argument("--fig", type=str, help="output figure")
    return parser.parse_args()


def plot_block_distr(Ls, Ns, fig_savename):
    """plot block length and frequency distributions"""

    Nmax = Ns.max()

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))

    kwargs = {"cumulative": True, "histtype": "step"}

    # block length distribution
    ax = axs[0]
    bins = np.logspace(0, 6, 100)

    ax.hist(Ls, bins=bins, color="C0", **kwargs)
    ax.set_ylabel("n. blocks", color="C0")
    ax.set_xlabel("block length (bp)")
    ax.set_title("cumulative block length distr.")

    axt = ax.twinx()
    axt.hist(Ls, weights=Ls, bins=bins, color="C3", **kwargs)
    axt.set_ylabel("block size", color="C3")

    ax.set_xscale("log")
    ax.set_xlim(1e1, 1e6)

    # block n. isolates distribution
    ax = axs[1]
    bins = np.arange(Nmax + 2) - 0.5

    ax.hist(Ns, bins=bins, color="C0", **kwargs)
    ax.set_ylabel("n. blocks", color="C0")
    ax.set_xlabel("n. isolates")
    ax.set_title("cumulative block frequency distr.")

    axt = ax.twinx()
    axt.hist(Ns, weights=Ls, bins=bins, color="C3", **kwargs)
    axt.set_ylabel("block size", color="C3")
    ax.set_xlim(0, Nmax + 1)

    plt.tight_layout()
    plt.savefig("figs/block_distr.png", facecolor="white", dpi=300)
    plt.close(fig)


if __name__ == "__main__":

    args = parse_args()

    # load pangraph
    pan = pp.Pangraph.load_json(args.pangraph)

    # dataframe with summary statistics for each block
    df = pan.to_blockstats_df()

    # vectors of block lengths and n. of strains in which the block is present
    Ls = df["len"].to_numpy()
    Ns = df["n. strains"].to_numpy()

    # plot block length and frequency distributions
    plot_block_distr(Ls, Ns, args.fig)
