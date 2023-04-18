# %%
from Bio import Phylo
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


# load newick tree
tree = Phylo.read("../data/coretree.nwk", "newick")
# %%

# load paf file
paf = pd.read_csv("../results/bla15/map.paf", sep="\t", header=None)

# %%
# get all leave names
leaves = [n.name for n in tree.get_terminals()]
# generate dictionary of colors for leaves
selected = paf[0].unique()


# yield sequential colors from colormap
def color_generator():
    cmap = plt.get_cmap("Dark2")
    for x in np.linspace(0, 1, len(selected)):
        yield cmap(x)


cg = color_generator()
colors = {l: next(cg) if l in selected else "black" for l in leaves}

# %%

# plot tree coloring only selected leaves
plt.figure(figsize=(6, 12))
ax = plt.gca()
Phylo.draw(
    tree,
    label_func=lambda x: x.name if x.name in leaves else "",
    do_show=False,
    axes=ax,
    label_colors={n.name: colors[n.name] for n in tree.get_terminals()},
)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig("figs/bla15_tree.png", facecolor="white", dpi=300)
plt.savefig("figs/bla15_tree.pdf")
plt.show()

cdf = pd.Series({k: mpl.colors.to_hex(v) for k, v in colors.items()}, name="color")
cdf.to_csv("leaves_colors.csv", index_label="strain")

# %%
