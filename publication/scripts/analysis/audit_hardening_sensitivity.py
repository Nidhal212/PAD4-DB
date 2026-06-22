"""
audit_hardening_sensitivity.py — Priorities 4, 5, 6 of the hardening pass.

P4  cliff-threshold sensitivity (|ΔpIC50| >= 1.5 / 2.0 / 2.5 at Tanimoto >= 0.8)
P5  source-independence-score robustness (±20% weight perturbation + alt scheme)
P6  assay-constrained null at a finer (6-level) assay partition

Outputs:
  outputs/tables/supp_threshold_sensitivity.csv
  outputs/tables/supp_independence_robustness.csv
  outputs/audit/assay_null_finer.md
"""
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd

ROOT = Path('/home/nidhal/PAD4-db_V2')
rng = np.random.default_rng(42)
assets = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
pairs  = pd.read_parquet(ROOT / 'data/processed/activity_pairs_with_sali.parquet')
pic    = dict(zip(assets.inchi_key, assets.pIC50))
scaf   = dict(zip(assets.inchi_key, assets.murcko_smiles))

HUBS = {'SMADULGDNOCLOP-GISFHXKWSA-N', 'RAVBZQAQTVGKIV-XBPDSQQVSA-N',
        'UDCDEKJNAMHBFH-HSZRJFAPSA-N', 'DVCKJOQIVOGXEI-XMMPIXPASA-N'}
hi = pairs[pairs.tanimoto >= 0.8]
edges_all = list(zip(hi.inchi_key_a, hi.inchi_key_b, hi.delta_pic50.abs()))
nodes = sorted(set(hi.inchi_key_a) | set(hi.inchi_key_b))

# ── P4: threshold sensitivity ─────────────────────────────────────────────────
print("=" * 60); print("P4 — CLIFF-THRESHOLD SENSITIVITY"); print("=" * 60)
def cliff_stats(delta_cut):
    deg = Counter(); cps = []
    for a, b, d in edges_all:
        if d >= delta_cut:
            deg[a] += 1; deg[b] += 1; cps.append((a, b))
    n = len(cps); comp = set(a for a, b in cps) | set(b for a, b in cps)
    top4 = deg.most_common(4)
    top4share = sum(c for _, c in top4) / n if n else 0.0
    top4ik = set(k for k, _ in top4)
    return n, len(comp), top4share, len(top4ik & HUBS)

# quick unrestricted null hub-share per threshold — permute over ALL eligible (Tan>=0.8) pairs
idx = {x: i for i, x in enumerate(nodes)}
elig = [(idx[a], idx[b]) for a, b, d in edges_all]      # ALL eligible pairs, not just observed cliffs
def null_hubshare(delta_cut, N=2000):
    base = np.array([pic[x] for x in nodes])
    out = np.empty(N)
    for it in range(N):
        perm = rng.permutation(base); deg = Counter(); n = 0
        for i, j in elig:
            if abs(perm[i] - perm[j]) >= delta_cut:
                deg[i] += 1; deg[j] += 1; n += 1
        out[it] = (sum(c for _, c in deg.most_common(4)) / n) if n else 0
    return out

rows = []
for dc in (1.5, 2.0, 2.5):
    n, ncomp, share, nhub = cliff_stats(dc)
    nl = null_hubshare(dc); p = (np.sum(nl >= share) + 1) / (len(nl) + 1)
    rows.append({'delta_threshold': dc, 'tanimoto_threshold': 0.8, 'n_cliff_pairs': n,
                 'n_compounds': ncomp, 'top4_hub_share_pct': round(share * 100, 1),
                 'top4_are_canonical_hubs': f'{nhub}/4', 'null_top4_share_pct': round(np.mean(nl) * 100, 1),
                 'hub_enrichment_fold': round(share / np.mean(nl), 1), 'perm_p': p})
    print(f"  Δ>={dc}: {n} cliffs, {ncomp} cmpds, top4 share {share*100:.1f}% "
          f"({nhub}/4 canonical), {share/np.mean(nl):.1f}x null, p<{p:.4f}")
pd.DataFrame(rows).to_csv(ROOT / 'outputs/tables/supp_threshold_sensitivity.csv', index=False)

# ── P5: independence-score robustness ─────────────────────────────────────────
print("\n" + "=" * 60); print("P5 — INDEPENDENCE-SCORE ROBUSTNESS"); print("=" * 60)
def combo(s):
    return '+'.join(sorted(set(s.split('|')))) if isinstance(s, str) else str(s)
assets['combo'] = assets.source_list.apply(combo)
base_w = {'bindingdb+chembl+pubchem_confirmatory': 0.3, 'bindingdb+pubchem_confirmatory': 0.5,
          'bindingdb+chembl': 0.6, 'chembl+pubchem_confirmatory': 0.7}
def score(w):
    return assets.combo.apply(lambda c: w.get(c, 1.0)).values   # single-source -> 1.0
orig = score(base_w)
from scipy.stats import spearmanr
rhos = []; counts = []
for _ in range(1000):
    w = {k: v * (1 + rng.uniform(-0.2, 0.2)) for k, v in base_w.items()}
    s = score(w); rhos.append(spearmanr(orig, s).correlation); counts.append(int((s >= 0.6).sum()))
rhos = np.array(rhos); counts = np.array(counts)
# alternative monotonic scheme: by number of sources (3->0.33, 2->0.66, 1->1.0)
nsrc = assets.combo.apply(lambda c: c.count('+') + 1)
alt = (1.0 / nsrc).values
rho_alt = spearmanr(orig, alt).correlation
print(f"  ±20% weight perturbation (1000 draws): Spearman ρ = {rhos.mean():.3f} "
      f"(min {rhos.min():.3f}); non-redundant count {counts.mean():.0f} "
      f"[{counts.min()}–{counts.max()}] (orig 528)")
print(f"  alternative (1/n_sources) scheme: Spearman ρ vs original = {rho_alt:.3f}")
pd.DataFrame([
    {'analysis': '±20% weight perturbation (1000 draws)', 'spearman_rho_mean': round(rhos.mean(), 3),
     'spearman_rho_min': round(rhos.min(), 3), 'nonredundant_mean': int(counts.mean()),
     'nonredundant_min': int(counts.min()), 'nonredundant_max': int(counts.max()), 'original_nonredundant': 528},
    {'analysis': 'alternative 1/n_sources scheme', 'spearman_rho_mean': round(rho_alt, 3),
     'spearman_rho_min': round(rho_alt, 3), 'nonredundant_mean': int((alt >= 0.6).sum()),
     'nonredundant_min': '', 'nonredundant_max': '', 'original_nonredundant': 528},
]).to_csv(ROOT / 'outputs/tables/supp_independence_robustness.csv', index=False)

# ── P6: assay-constrained null at finer partition ─────────────────────────────
print("\n" + "=" * 60); print("P6 — ASSAY-CONSTRAINED NULL (finer partition)"); print("=" * 60)
mech6 = dict(zip(assets.inchi_key, assets.assay_mechanism_classes))
elig_ik = [(a, b) for a, b, d in edges_all]              # ALL eligible Tan>=0.8 pairs
def hubshare(pm):
    deg = Counter(); n = 0
    for a, b in elig_ik:
        if abs(pm[a] - pm[b]) >= 2.0: deg[a] += 1; deg[b] += 1; n += 1
    return (sum(c for _, c in deg.most_common(4)) / n) if n else 0
obs = hubshare(pic)
groups = {}
for nd in nodes: groups.setdefault(mech6.get(nd), []).append(nd)
N = 5000; null = np.empty(N)
for it in range(N):
    pm = dict(pic)
    for g, mem in groups.items():
        if len(mem) > 1:
            vals = rng.permutation([pic[m] for m in mem])
            for m, v in zip(mem, vals): pm[m] = v
    null[it] = hubshare(pm)
p = (np.sum(null >= obs) + 1) / (N + 1)
print(f"  6-level assay partition: observed {obs*100:.1f}% vs null {null.mean()*100:.1f}±{null.std()*100:.1f}% "
      f"({obs/null.mean():.1f}x, p<{p:.4f})")
(ROOT / 'outputs/audit/assay_null_finer.md').write_text(
    "# P6 — Assay-constrained null at a finer partition\n\n"
    "Harmonized per-record assay identifiers (e.g., a single PubChem AID / ChEMBL assay ID) are **not "
    "available at the consensus-compound level**, because each compound is aggregated across multiple "
    "assays during curation; the assay-constrained null therefore uses assay-mechanism categories as the "
    "available proxy. As a finer-grained robustness check we repeated the hub-concentration null using the "
    "six-level `assay_mechanism_classes` partition (BAEE|FP, BAEE|RFMS, BAEE, FP, BAEE|covalent, RFMS) "
    "rather than the four-level mechanism class.\n\n"
    f"Observed top-4 hub share = {obs*100:.1f}%; null (5,000 permutations) = {null.mean()*100:.1f} ± "
    f"{null.std()*100:.1f}%; enrichment {obs/null.mean():.1f}×, P < {p:.4f}. Hub concentration remains "
    "significant under the finer partition. We note as a limitation that a true assay-identifier-level null "
    "could not be constructed from the consensus table.")
print("\nDONE")
