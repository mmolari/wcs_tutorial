import argparse
import numpy as np
import pandas as pd
import pypangraph as pp
from itertools import combinations


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pangraph", type=str, help="pangraph file")
    parser.add_argument(
        "--dist_df", type=str, help="output pairwise distance dataframe"
    )
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()

    # load pangraph
    pan = pp.Pangraph.load_json(args.pangraph)

    # create a presence-absence dataframe for blocks
    pa_df = pan.to_blockcount_df() > 0

    # create a vector of block lengths. The order is the same as in the presence-absence dataframe
    Ls = pan.to_blockstats_df()["len"].to_dict()
    Ls = [Ls[b] for b in pa_df.columns]

    # compute pairwise private sequence distance
    dist_df = []
    for p1, p2 in combinations(pan.paths, 2):

        # names of isolates
        iso1 = p1.name
        iso2 = p2.name

        # presence-absence binary vector for each isolate
        vi = pa_df.loc[iso1].to_numpy()
        vj = pa_df.loc[iso2].to_numpy()

        # total length of private sequence
        d = np.sum((vi ^ vj) * Ls)

        # append to dataframe
        dist_df.append([iso1, iso2, d])
        dist_df.append([iso2, iso1, d])

    # add distance zero for self-comparisons
    for p in pan.paths:
        dist_df.append([p.name, p.name, 0])

    # transform to dataframe
    dist_df = pd.DataFrame(dist_df, columns=["p1", "p2", "private_seq"])

    # save to file
    dist_df.to_csv(args.dist_df, index=False)
