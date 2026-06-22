"""
audit_constrained_permutation.py — Phase 6

Robustness of the cliff-rarity and hub-concentration findings under three null
models of increasing strictness:
  (1) unrestricted   — pIC50 shuffled across the whole Tanimoto>=0.8 subgraph
  (2) scaffold-const — pIC50 shuffled only WITHIN each Murcko scaffold series
  (3) assay-const    — pIC50 shuffled only WITHIN each assay-mechanism class

For each: severe-cliff count (rarity) and top-4 cliff-degree share (hubness),
with empirical p, z-score, and fold-change vs observed.

Output: outputs/audit/null_model_comparison.md  (+ .csv)
"""
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd

ROOT = Path('/home/nidhal/PAD4-db_V2')
rng = np.random.default_rng(42)
N = 10000

assets = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
pairs  = pd.read_parquet(ROOT / 'data/processed/activity_pairs_with_sali.parquet')

pic   = dict(zip(assets.inchi_key, assets.pIC50))
scaf  = dict(zip(assets.inchi_key, assets.murcko_smiles))
mech  = dict(zip(assets.inchi_key, assets.mechanism_class))

hi = pairs[pairs.tanimoto >= 0.8]
edges = list(zip(hi.inchi_key_a, hi.inchi_key_b))
nodes = sorted(set(hi.inchi_key_a) | set(hi.inchi_key_b))

def metrics(pmap):
    deg = Counter(); n = 0
    for a, b in edges:
        if abs(pmap[a] - pmap[b]) >= 2.0:
            deg[a] += 1; deg[b] += 1; n += 1
    share = sum(d for _, d in deg.most_common(4)) / n if n else 0.0
    return n, share

obs_n, obs_share = metrics(pic)

# group index for constrained shuffles (only groups with >=2 members are permutable)
def grouped_perm(label_map):
    groups = {}
    for nd in nodes:
        groups.setdefault(label_map.get(nd), []).append(nd)
    base = {nd: pic[nd] for nd in nodes}
    def one():
        pm = dict(base)
        for g, members in groups.items():
            if len(members) > 1:
                vals = rng.permutation([base[m] for m in members])
                for m, v in zip(members, vals):
                    pm[m] = v
        return pm
    return one

def unrestricted_perm():
    vals = np.array([pic[nd] for nd in nodes])
    def one():
        return dict(zip(nodes, rng.permutation(vals)))
    return one

def run(perm_fn, name):
    cnt = np.empty(N); shr = np.empty(N)
    for i in range(N):
        pm = perm_fn(); cnt[i], shr[i] = metrics(pm)
    p_n = (np.sum(cnt <= obs_n) + 1) / (N + 1)
    p_s = (np.sum(shr >= obs_share) + 1) / (N + 1)
    z_s = (obs_share - shr.mean()) / (shr.std() + 1e-9)
    print(f"  [{name}] cliffs null {cnt.mean():.0f}±{cnt.std():.0f} "
          f"(obs {obs_n}, {cnt.mean()/obs_n:.1f}x, p<{p_n:.4f}) | "
          f"hub-share null {shr.mean()*100:.1f}±{shr.std()*100:.1f}% "
          f"(obs {obs_share*100:.1f}%, {obs_share/shr.mean():.1f}x, z={z_s:.1f}, p<{p_s:.4f})")
    return {
        'null_model': name,
        'obs_cliffs': obs_n, 'null_cliffs_mean': round(cnt.mean(), 1), 'null_cliffs_sd': round(cnt.std(), 1),
        'cliff_fold_depletion': round(cnt.mean() / obs_n, 1), 'cliff_p': p_n,
        'obs_hub_share_pct': round(obs_share * 100, 1),
        'null_hub_share_mean_pct': round(shr.mean() * 100, 1), 'null_hub_share_sd_pct': round(shr.std() * 100, 1),
        'hub_fold_enrichment': round(obs_share / shr.mean(), 1), 'hub_z': round(z_s, 1), 'hub_p': p_s,
    }

print("=" * 64)
print("PHASE 6 — NULL MODEL COMPARISON (10,000 permutations each)")
print(f"  observed: {obs_n} severe cliffs, top-4 hub share {obs_share*100:.1f}%")
print("=" * 64)
res = [
    run(unrestricted_perm(), 'unrestricted'),
    run(grouped_perm(scaf), 'scaffold_constrained'),
    run(grouped_perm(mech), 'assay_constrained'),
]
df = pd.DataFrame(res)
df.to_csv(ROOT / 'outputs/audit/null_model_comparison.csv', index=False)

lines = ["# Phase 6 — Null-model comparison for cliff rarity & hub concentration", "",
         f"Observed on the Tanimoto≥0.8 subgraph ({len(nodes)} compounds, {len(edges)} eligible pairs): "
         f"**{obs_n} severe cliffs**, **top-4 hub share {obs_share*100:.1f}%**. 10,000 permutations per model.", "",
         "| Null model | Cliffs (null) | Cliff depletion | p | Hub share (null) | Hub enrichment | z | p |",
         "|---|---|---|---|---|---|---|---|"]
for r in res:
    lines.append(f"| {r['null_model']} | {r['null_cliffs_mean']}±{r['null_cliffs_sd']} | "
                 f"{r['cliff_fold_depletion']}× | <{r['cliff_p']:.4f} | "
                 f"{r['null_hub_share_mean_pct']}±{r['null_hub_share_sd_pct']}% | "
                 f"{r['hub_fold_enrichment']}× | {r['hub_z']} | <{r['hub_p']:.4f} |")
lines += ["",
          "**Interpretation.** The scaffold-constrained null is the strictest test: it shuffles "
          "potency only within Murcko scaffold series, so any apparent hub concentration arising purely "
          "from scaffold structure is absorbed into the null. Hub enrichment and cliff rarity that survive "
          "this null are properties of the potency assignment, not of the similarity topology."]
(ROOT / 'outputs/audit/null_model_comparison.md').write_text("\n".join(lines))
print("\nSaved: outputs/audit/null_model_comparison.{csv,md}")
print("DONE")
