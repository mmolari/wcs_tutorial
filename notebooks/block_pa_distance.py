# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import Phylo
import yaml

from itertools import combinations

import pypangraph as pp

tree_file = "../data/coretree.nwk"
config_file = "../config/config.yaml"
pangraph_file = "../results/pangraph/subset.json"
# %%

# Parse the tree file (in Newick format) and get the root node
tree = Phylo.read(tree_file, "newick")

# List of leaf labels to keep
with open(config_file, "r") as f:
    config = yaml.safe_load(f)
selected_isolates = config["strains-subset"]

# Prune the tree to keep only the selected leaves
to_prune = []
for leaf in tree.get_terminals():
    if leaf.name not in selected_isolates:
        to_prune.append(leaf)

for leaf in to_prune:
    tree.prune(leaf)

# Write the pruned tree to a new file
Phylo.draw(tree)

# %%
pan = pp.Pangraph.load_json(pangraph_file)

# %%
df = pan.to_blockstats_df()
pa_df = pan.to_blockcount_df() > 0

# %%


fig, axs = plt.subplots(1, 2, figsize=(10, 5))

kwargs = {"cumulative": True, "histtype": "step"}

ax = axs[0]
bins = np.logspace(0, 6, 100)
ax.hist(df["len"], bins=bins, color="C0", **kwargs)
axt = ax.twinx()
axt.hist(df["len"], weights=df["len"], bins=bins, color="C3", **kwargs)
ax.set_ylabel("n. blocks", color="C0")
axt.set_ylabel("block size", color="C3")
ax.set_xscale("log")
ax.set_title("cumulative block length distr.")
ax.set_xlim(1e1, 1e6)
ax.set_xlabel("block length (bp)")

ax = axs[1]
bins = np.arange(len(selected_isolates) + 2) - 0.5
ax.hist(df["n. strains"], bins=bins, color="C0", **kwargs)
axt = ax.twinx()
axt.hist(df["n. strains"], weights=df["len"], bins=bins, color="C3", **kwargs)
ax.set_ylabel("n. blocks", color="C0")
axt.set_ylabel("block size", color="C3")
ax.set_title("cumulative block frequency distr.")
ax.set_xlim(0, 11)
ax.set_xlabel("n. isolates")


plt.tight_layout()
plt.savefig("figs/block_distr.png", facecolor="white", dpi=300)
plt.show()
# %%

Ls = df["len"].to_dict()
Ls = [Ls[b] for b in pa_df.columns]

dist_df = []
for i, j in combinations(selected_isolates, 2):
    vi = pa_df.loc[i].to_numpy()
    vj = pa_df.loc[j].to_numpy()
    d = np.sum((vi ^ vj) * Ls)
    d_tree = tree.distance(i, j)
    dist_df.append([i, j, d, d_tree])
    dist_df.append([j, i, d, d_tree])
dist_df = pd.DataFrame(dist_df, columns=["i", "j", "d_pa", "d_tree"])
dist_M = dist_df.pivot(index="i", columns="j", values="d_pa").fillna(0)

order = [x.name for x in tree.get_terminals()]
dist_M = dist_M.loc[order, order]
# %%

fig, axs = plt.subplots(1, 2, figsize=(10, 5))

ax = axs[0]
Phylo.draw(
    tree,
    axes=ax,
    do_show=False,
    label_func=lambda x: x.name if x.name in order else "",
    show_confidence=False,
)
for k in ["top", "right", "left"]:
    ax.spines[k].set_visible(False)
ax.set_ylabel("")
ax.set_yticks([])

ax = axs[1]
g = ax.matshow(dist_M / 1000)
ax.set_xticks(np.arange(len(order)))
ax.set_xticklabels(order, rotation=90)
ax.set_yticks(np.arange(len(order)))
ax.set_yticklabels([])
plt.colorbar(g, ax=ax, label="private seq. (kbp)", fraction=0.046, pad=0.04)
plt.tight_layout()
fig.subplots_adjust(wspace=-0.2, top=0.65)
plt.savefig("figs/private_seq_vs_tree.png", facecolor="white", dpi=300)
plt.show()

# %%
plt.figure(figsize=(4, 3))
sdf = dist_df[dist_df["i"] > dist_df["j"]]
plt.scatter(sdf["d_tree"], sdf["d_pa"] / 1000, alpha=0.4)
plt.xlabel("core genome tree distance")
plt.ylabel("private sequence distance (kbp)")
for s in ["top", "right"]:
    plt.gca().spines[s].set_visible(False)
plt.tight_layout()
plt.savefig("figs/pa_dist_vs_tree_scatter.png", facecolor="white", dpi=300)
plt.show()
# %%
