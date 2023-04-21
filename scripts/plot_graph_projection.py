import argparse
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

import pypangraph as pp
from pypangraph.pangraph_projector import PanProjector
from pypangraph.visualization_projection import draw_projection


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pangraph",
        help="Path to the pangraph file",
        type=str,
    )
    parser.add_argument(
        "--fig",
        help="Path to the figure file",
        type=str,
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    pan = pp.Pangraph.load_json(args.pangraph)
    i1, i2 = pan.strains()

    # create projector
    ppj = PanProjector(pan)

    # project over the pair
    pr = ppj.project(i1, i2, exclude_dupl=False)

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))

    cdict = defaultdict(lambda: plt.get_cmap("rainbow")(np.random.rand()))
    draw_projection(
        pr,
        ax=ax,
        color_dict=cdict,
    )
    ax.legend()
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

    plt.tight_layout()
    plt.savefig(args.fig, dpi=300, facecolor="white")
    plt.close(fig)
