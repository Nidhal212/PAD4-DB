"""
audit_within_assay_and_mw.py — Reviewer-hardening Tasks 1 & 2.

T1  Within-assay cliff robustness: of the 94 severe cliffs, how many persist
    (|ΔpIC50| >= 2.0) using per-AID potencies from the SAME assay?
T2  ΔMW sensitivity: are severe cliffs driven by small or large mass changes?

Outputs:
  outputs/audit/within_assay_feasibility.md
  outputs/tables/supp_within_assay.csv          (Table S-assay)
  outputs/tables/supp_mw_sensitivity.csv         (Table S-mw)
"""
from pathlib import Path
from collections import Counter, defaultdict
import numpy as np
import pandas as pd

ROOT = Path('/home/nidhal/PAD4-db_V2')
assets = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
ac = pd.read_parquet(ROOT / 'data/processed/activity_cliffs.parquet')
sev = ac[ac.cliff_tier == 'severe'].copy()
HUBS = {'SMADULGDNOCLOP-GISFHXKWSA-N', 'RAVBZQAQTVGKIV-XBPDSQQVSA-N',
        'UDCDEKJNAMHBFH-HSZRJFAPSA-N', 'DVCKJOQIVOGXEI-XMMPIXPASA-N'}

# ── T1: within-assay ──────────────────────────────────────────────────────────
ps = pd.read_parquet(ROOT / 'data/interim/normalized/potency_space.parquet')
# per-compound: aid -> pic50_aid (mean if duplicated)
aid_pic = defaultdict(dict)
for r in ps.itertuples():
    if pd.notna(r.pic50_aid):
        aid_pic[r.inchi_key][r.aid] = r.pic50_aid
fam = dict(zip(assets.inchi_key, assets.assay_mechanism_classes))

rows = []
for r in sev.itertuples():
    a, b = r.inchi_key_a, r.inchi_key_b
    common = set(aid_pic.get(a, {})) & set(aid_pic.get(b, {}))
    best = max((abs(aid_pic[a][k] - aid_pic[b][k]) for k in common), default=np.nan)
    fa = set(str(fam.get(a, '')).split('|')); fb = set(str(fam.get(b, '')).split('|'))
    rows.append({'a': a[:14], 'b': b[:14], 'consensus_delta': round(abs(r.delta_pic50), 3),
                 'n_common_aids': len(common), 'within_aid_max_delta': round(best, 3) if common else np.nan,
                 'within_aid_confirmed': bool(common) and best >= 2.0,
                 'shares_assay_family': bool(fa & fb and fa != {''}),
                 'involves_hub': a in HUBS or b in HUBS})
W = pd.DataFrame(rows)
n = len(W)
n_common = int((W.n_common_aids > 0).sum())
n_conf = int(W.within_aid_confirmed.sum())
n_fam = int(W.shares_assay_family.sum())

feasible = n_common > 0
fz = (ROOT / 'outputs/audit/within_assay_feasibility.md')
fz.write_text(
    "# T1 — Within-assay cliff robustness: feasibility\n\n"
    "**Outcome A — FEASIBLE.** Per-assay potencies (`pic50_aid`, keyed by PubChem AID / ChEMBL assay) are "
    "retained in `data/interim/normalized/potency_space.parquet`, so each severe-cliff pair can be tested for "
    "persistence within a single shared assay. Note: the consensus pIC50 used elsewhere aggregates across "
    "assays, but the pre-aggregation per-AID layer survives and is used here.\n\n"
    f"- Severe cliff pairs: {n}\n"
    f"- Pairs with >=1 assay (AID) measuring BOTH compounds: {n_common} ({n_common/n*100:.0f}%)\n"
    f"- Pairs that remain severe (|ΔpIC50| >= 2.0) within a shared assay: {n_conf} ({n_conf/n*100:.0f}%)\n"
    f"- Pairs whose two compounds share at least one assay *family*: {n_fam} ({n_fam/n*100:.0f}%)\n\n"
    "Limitation: many cliff pairs draw their two members from different assays/sources, so a shared-AID test "
    "is only possible for the subset measured in a common assay; for that subset it is a strong, "
    "aggregation-free check.")

# Table S-assay
tbl = pd.DataFrame([
    {'metric': 'severe cliff pairs', 'value': n, 'pct_of_94': 100.0},
    {'metric': 'pairs with a common assay (AID)', 'value': n_common, 'pct_of_94': round(n_common/n*100, 1)},
    {'metric': 'pairs confirmed severe within a common assay', 'value': n_conf, 'pct_of_94': round(n_conf/n*100, 1)},
    {'metric': 'of pairs with a common assay, % still severe',
     'value': round(n_conf/n_common*100, 1) if n_common else 0, 'pct_of_94': ''},
    {'metric': 'pairs sharing >=1 assay family', 'value': n_fam, 'pct_of_94': round(n_fam/n*100, 1)},
])
tbl.to_csv(ROOT / 'outputs/tables/supp_within_assay.csv', index=False)
W.to_csv(ROOT / 'outputs/audit/within_assay_per_pair.csv', index=False)

# hub robustness among confirmed
conf = W[W.within_aid_confirmed]
hub_share_conf = (conf.involves_hub.sum() / len(conf) * 100) if len(conf) else 0
print("=" * 60); print("T1 — WITHIN-ASSAY"); print("=" * 60)
print(f"  {n_conf}/{n} severe cliffs confirmed within a common assay "
      f"({n_common} had a common assay; {n_conf}/{n_common}={n_conf/max(n_common,1)*100:.0f}% of those persist)")
print(f"  hub-involved among confirmed: {conf.involves_hub.sum()}/{len(conf)} ({hub_share_conf:.0f}%)")

# ── T2: ΔMW sensitivity ───────────────────────────────────────────────────────
mw = dict(zip(assets.inchi_key, assets.mol_weight))
recs = []
for r in sev.itertuples():
    dmw = abs(mw[r.inchi_key_a] - mw[r.inchi_key_b])
    cat = 'A: <20 Da' if dmw < 20 else ('B: 20-50 Da' if dmw < 50 else 'C: >=50 Da')
    recs.append({'dmw': dmw, 'dpic50': abs(r.delta_pic50), 'cat': cat})
D = pd.DataFrame(recs)
mwt = []
for cat in ['A: <20 Da', 'B: 20-50 Da', 'C: >=50 Da']:
    g = D[D.cat == cat]
    mwt.append({'dMW_category': cat, 'n_cliffs': len(g), 'pct': round(len(g)/len(D)*100, 1),
                'mean_dpIC50': round(g.dpic50.mean(), 3) if len(g) else np.nan,
                'median_dpIC50': round(g.dpic50.median(), 3) if len(g) else np.nan,
                'max_dpIC50': round(g.dpic50.max(), 3) if len(g) else np.nan})
MW = pd.DataFrame(mwt)
MW.to_csv(ROOT / 'outputs/tables/supp_mw_sensitivity.csv', index=False)
print("\n" + "=" * 60); print("T2 — ΔMW SENSITIVITY"); print("=" * 60)
print(MW.to_string(index=False))
small = D[D.dmw < 20]
print(f"\n  severe cliffs with ΔMW < 20 Da: {len(small)}/{len(D)} ({len(small)/len(D)*100:.0f}%), "
      f"mean |ΔpIC50| {small.dpic50.mean():.2f} (max {small.dpic50.max():.2f})")
print(f"  Spearman ΔMW vs |ΔpIC50|: {D[['dmw','dpic50']].corr(method='spearman').iloc[0,1]:.3f}")
print("DONE")
