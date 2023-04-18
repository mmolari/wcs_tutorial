# %%
import matplotlib.pyplot as plt
import numpy as np
import pypangraph as pp
import yaml

from Bio import Phylo
from collections import defaultdict
from pypangraph.pangraph_projector import PanProjector
from pypangraph.visualization_projection import draw_projection


tree_file = "../data/coretree.nwk"
config_file = "../config/config.yaml"
pangraph_file = "../results/pangraph/marginalized.json"
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

Phylo.draw(tree)
# %%

pan = pp.Pangraph.load_json(pangraph_file)
i1, i2 = pan.strains()
# %%

# create projector
ppj = PanProjector(pan)

# project over the pair
pr = ppj.project(i1, i2, exclude_dupl=False)


fig, axs = plt.subplots(1, 2, figsize=(10, 6), gridspec_kw={"width_ratios": [1, 4]})

ax = axs[0]
Phylo.draw(
    tree,
    axes=ax,
    label_func=lambda x: x.name if x.name in selected_isolates else "",
    do_show=False,
    label_colors={i1: "C0", i2: "C1"},
)
ax.set_xticks([])
ax.set_yticks([])
for k in ["top", "right", "left", "bottom"]:
    ax.spines[k].set_visible(False)

ax = axs[1]
cdict = defaultdict(lambda: plt.get_cmap("rainbow")(np.random.rand()))
draw_projection(
    pr,
    ax=ax,
    color_dict=cdict,
)
ax.legend()
ax.set_xticks([])
ax.set_yticks([])
for k in ["top", "right", "left", "bottom"]:
    ax.spines[k].set_visible(False)

plt.tight_layout()
plt.savefig("figs/marginal_graph.png", dpi=300, facecolor="white")
plt.show()

# %%
