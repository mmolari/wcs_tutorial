# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypangraph as ppg
from collections import defaultdict
from itertools import combinations, product

# %%

pan = ppg.Pangraph.load_json("../pangraph_bla.json")

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


res = []
for pi, pj in product(pan.paths, pan.paths):
    S = path_anchor_distance(pi, pj, core_anchor, L)
    res.append(
        {
            "pi": pi.name,
            "pj": pj.name,
            "d": S,
        }
    )
S = pd.DataFrame(res).pivot(index="pi", columns="pj", values="d")

# %%
import scipy.cluster.hierarchy as spc

D = S.to_numpy()
T = spc.linkage(S)
d = spc.dendrogram(T)
order = S.index[d["leaves"]]
plt.matshow(S.loc[order, order])


# %%
cz = 0
def c_idx():
    global cz
    cz += 1
    return cz
cm = plt.cm.get_cmap("tab20")
colors = defaultdict(lambda: cm(c_idx()))


plt.figure(figsize=(10, 7))

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
plt.yticks(range(len(order)), [o.split("-")[0] for o in order])

plt.show()
# %%
