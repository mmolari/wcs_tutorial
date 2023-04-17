# %%
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypangraph as ppg
from Bio import Phylo
from collections import defaultdict
from itertools import combinations, product

# %%
window = 5000
pangraph_file = f"../results/pangraph/window__bla15__{window}.json"
pan = ppg.Pangraph.load_json(pangraph_file)

# %%

path_d = pan.to_paths_dict()
pids = list(path_d.keys())
df = pan.to_blockstats_df()
core_anchor = df[df.core].index[0]
L = df["len"].to_dict()


def path_anchor_distance(pi, pj, a, len_dict):
    li, lj = [list(p.block_ids) for p in [pi, pj]]
    si, sj = [list(p.block_strands) for p in [pi, pj]]

    xi, xj = [l.index(a) for l in [li, lj]]
    if not si[xi]:
        li = li[::-1]
        si = [not s for s in reversed(si)]
        xi = li.index(a)
    if not sj[xj]:
        lj = lj[::-1]
        sj = [not s for s in reversed(sj)]
        xj = lj.index(a)

    S = 0
    # fwd comparison
    for k in range(0, min(xi, xj)):
        I, J = xi - k, xj - k
        if (li[I] == lj[J]) and (si[I] == sj[J]):
            S += len_dict[li[I]]
        else:
            break
    # rev comparison
    for k in range(1, min(len(li) - xi, len(lj) - xj)):
        I, J = xi + k, xj + k
        if (li[I] == lj[J]) and (si[I] == sj[J]):
            S += len_dict[li[I]]
        else:
            break

    return S


# load phylogenetic tree
tree = Phylo.read("../data/coretree.nwk", "newick")

# dictionary of leaves names to leaves
leaves = {n.name: n for n in tree.get_terminals()}


res = []
for pi, pj in product(pan.paths, pan.paths):
    S = path_anchor_distance(pi, pj, core_anchor, L)
    si = pi.name.split("-")[0]
    sj = pj.name.split("-")[0]
    res.append(
        {
            "pi": pi.name,
            "pj": pj.name,
            "shared_path": S,
            "tree_dist": tree.distance(leaves[si], leaves[sj]),
        }
    )
dist_df = pd.DataFrame(res)
S = dist_df.pivot(index="pi", columns="pj", values="shared_path")


# %%
sdf = dist_df[dist_df.pi > dist_df.pj]
plt.scatter(sdf.shared_path, sdf.tree_dist, alpha=0.2)
plt.xlabel("Shared path length")
plt.ylabel("Tree distance")
plt.tight_layout()
plt.savefig("figs/bla15_window_scatter.png", facecolor="white", dpi=300)
plt.savefig("figs/bla15_window_scatter.pdf")
plt.show()

# load hex colors from csv
leaves_colors = pd.read_csv("leaves_colors.csv", index_col=0)["color"].to_dict()

# %%
import scipy.cluster.hierarchy as spc

D = S.to_numpy()
T = spc.linkage(S)
d = spc.dendrogram(T)
order = S.index[d["leaves"]]
plt.matshow(S.loc[order, order])


# %%


def color_generator():
    cm = plt.cm.get_cmap("tab20")
    c_idx = 0
    while True:
        yield cm(c_idx)
        c_idx += 1


cg = color_generator()
colors = defaultdict(lambda: next(cg))

plt.figure(figsize=(12, 7))
ax = plt.gca()
for n, p_name in enumerate(order):
    p = pan.paths[p_name]
    l = 0
    cid = list(p.block_ids).index(core_anchor)
    s = p.block_strands[cid]
    bids = list(p.block_ids) if s else list(reversed(p.block_ids))
    for b in bids:
        lb = L[b]
        c = "lightgray" if df["count"][b] == 1 else colors[b]
        plt.plot([l, l + lb], [n, n], color=c, lw=5)
        l += lb
    # add text
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

    # plt.scatter(-300, n, marker="o", c=leaves_colors[strain])
# plt.yticks(range(len(order)), [o.split("-")[0] for o in order])
plt.yticks([])
for sp in ["top", "right", "left"]:
    ax.spines[sp].set_visible(False)
plt.tight_layout()
plt.savefig("figs/bla15_window.png", facecolor="white", dpi=300)
plt.savefig("figs/bla15_window.pdf")
plt.show()
# %%
cdf = pd.Series({k: mpl.colors.to_hex(v) for k, v in colors.items()}, name="Colour")
for b in pan.blocks:
    if not b.id in cdf:
        cdf[b.id] = mpl.colors.to_hex("lightgray")
cdf.to_csv("bandage_colors.csv", index_label="Name")
# %%
