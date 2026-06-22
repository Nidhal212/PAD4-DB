"""
audit_final_hardening.py — Final reviewer-risk pass (M5 within-assay, M18 ΔMW).

No new modeling. Uses existing standardized data + assay provenance.

Outputs:
  outputs/audit/within_assay_robustness.csv     (Table S-assay)
  outputs/audit/mw_sensitivity_analysis.csv      (Table S-mw)
  outputs/audit/representative_cliffs_for_lit.csv (Task 3 selection, data-only)
  prints feasibility + interpretation
"""
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT = ROOT / 'outputs/audit'; OUT.mkdir(parents=True, exist_ok=True)
HUBS = {'SMADULGDNOCLOP-GISFHXKWSA-N', 'RAVBZQAQTVGKIV-XBPDSQQVSA-N',
        'UDCDEKJNAMHBFH-HSZRJFAPSA-N', 'DVCKJOQIVOGXEI-XMMPIXPASA-N'}

ac    = pd.read_parquet(ROOT / 'data/processed/activity_cliffs.parquet')
sev   = ac[ac.cliff_tier == 'severe'].copy().reset_index(drop=True)
norm  = pd.read_parquet(ROOT / 'data/interim/normalized/normalized_activities.parquet')
pad   = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
assert len(sev) == 94

def top4_share(pairs_df):
    deg = Counter()
    for _, r in pairs_df.iterrows():
        deg[r.inchi_key_a] += 1; deg[r.inchi_key_b] += 1
    if len(pairs_df) == 0:
        return 0.0, []
    top = deg.most_common(4)
    return sum(d for _, d in top) / len(pairs_df), top

print("=" * 64)
print("TASK 1 — WITHIN-ASSAY CLIFF ROBUSTNESS (M5)")
print("=" * 64)

# compound -> set of assay identifiers (AID as string; numeric=PubChem AID,
# 'CHEMBL6111'/'Q9UM07' = coarse ChEMBL/BindingDB source buckets)
norm_q = norm[norm['aid'].notna()].copy()
norm_q['aid'] = norm_q['aid'].astype(str)
ik2aids = norm_q.groupby('inchi_key')['aid'].apply(lambda s: set(s.dropna())).to_dict()
# distinguish genuine PubChem assays (numeric) from coarse source buckets
def numeric_aids(s): return {a for a in s if a.isdigit()}
ik2pcaids = {ik: numeric_aids(a) for ik, a in ik2aids.items()}
sev_iks = set(sev.inchi_key_a) | set(sev.inchi_key_b)
with_pc = sum(1 for ik in sev_iks if ik2pcaids.get(ik))
print(f"  Feasibility: {len(norm_q)} measurement rows carry an assay identifier.")
print(f"  Of {len(sev_iks)} severe-cliff compounds, {with_pc} have >=1 genuine PubChem AID "
      f"(numeric); ChEMBL/BindingDB use single coarse source buckets (CHEMBL6111 / Q9UM07).")

def shared(a, b, m):
    return len(m.get(a, set()) & m.get(b, set())) > 0

# strict = share a genuine PubChem AID; any = share any assay identifier incl. coarse buckets
sev['same_pc_aid'] = sev.apply(lambda r: shared(r.inchi_key_a, r.inchi_key_b, ik2pcaids), axis=1)
sev['same_assay']  = sev.apply(lambda r: shared(r.inchi_key_a, r.inchi_key_b, ik2aids), axis=1)
n_pc   = int(sev.same_pc_aid.sum())
n_same = int(sev.same_assay.sum())
print(f"  Severe cliffs sharing a genuine PubChem AID (strict same-assay): {n_pc}/94 ({n_pc/94*100:.0f}%)")
print(f"  Severe cliffs sharing any assay identifier incl. source buckets: {n_same}/94 ({n_same/94*100:.0f}%)")

obs_share, _ = top4_share(sev)
sa_share, sa_top = top4_share(sev[sev.same_assay])
print(f"  Hub top-4 share — all 94: {obs_share*100:.1f}% ; same-assay subset (n={n_same}): {sa_share*100:.1f}%")
hubs_present = [ik[:14] for ik, _ in sa_top if ik in HUBS]
print(f"  Hubs retained in same-assay subset top-4: {hubs_present}")

pd.DataFrame([
    {'metric': 'Severe cliff pairs', 'all_94': 94, 'same_assay_any': n_same, 'same_pubchem_aid': n_pc},
    {'metric': 'Percent of 94', 'all_94': 100.0, 'same_assay_any': round(n_same/94*100, 1), 'same_pubchem_aid': round(n_pc/94*100, 1)},
    {'metric': 'Top-4 hub share (%)', 'all_94': round(obs_share*100, 1), 'same_assay_any': round(sa_share*100, 1), 'same_pubchem_aid': '—'},
]).to_csv(OUT / 'within_assay_robustness.csv', index=False)
print(f"  Saved: {OUT/'within_assay_robustness.csv'}")

print("\n" + "=" * 64)
print("TASK 2 — ΔMW SENSITIVITY (M18)")
print("=" * 64)
mw = dict(zip(pad.inchi_key, pad.mol_weight))
sev['dMW'] = sev.apply(lambda r: abs(mw.get(r.inchi_key_a, np.nan) - mw.get(r.inchi_key_b, np.nan)), axis=1)
sev['dpic'] = sev.delta_pic50.abs()

def cat(d):
    if d < 20: return 'A: <20 Da'
    if d < 50: return 'B: 20-50 Da'
    return 'C: >=50 Da'
sev['mw_cat'] = sev.dMW.apply(cat)

rows = []
for c in ['A: <20 Da', 'B: 20-50 Da', 'C: >=50 Da']:
    sub = sev[sev.mw_cat == c]
    rows.append({'dMW category': c, 'count': len(sub), 'pct': round(len(sub)/94*100, 1),
                 'mean |dpIC50|': round(sub.dpic.mean(), 2) if len(sub) else np.nan,
                 'median |dpIC50|': round(sub.dpic.median(), 2) if len(sub) else np.nan,
                 'max |dpIC50|': round(sub.dpic.max(), 2) if len(sub) else np.nan})
mwt = pd.DataFrame(rows)
print(mwt.to_string(index=False))
mwt.to_csv(OUT / 'mw_sensitivity_analysis.csv', index=False)
small = sev[sev.dMW < 20]
print(f"\n  KEY: {len(small)}/94 severe cliffs ({len(small)/94*100:.0f}%) involve ΔMW < 20 Da; "
      f"their mean |ΔpIC50| = {small.dpic.mean():.2f} (max {small.dpic.max():.2f}).")
print(f"  Median ΔMW across all 94 cliffs: {sev.dMW.median():.1f} Da.")
print(f"  Saved: {OUT/'mw_sensitivity_analysis.csv'}")

print("\n" + "=" * 64)
print("TASK 3 — REPRESENTATIVE CLIFF SELECTION (data-only; literature = author curation)")
print("=" * 64)
sev['sali'] = sev.dpic / (1 - sev.tanimoto)
sev['involves_hub'] = sev.apply(lambda r: r.inchi_key_a in HUBS or r.inchi_key_b in HUBS, axis=1)
# pick: top-3 SALI hub-associated + top-3 SALI non-hub + top-2 highest dpIC50
pick = pd.concat([
    sev[sev.involves_hub].nlargest(3, 'sali'),
    sev[~sev.involves_hub].nlargest(3, 'sali'),
    sev.nlargest(2, 'dpic'),
]).drop_duplicates(subset=['inchi_key_a', 'inchi_key_b']).head(8)
out3 = pick[['inchi_key_a', 'inchi_key_b', 'tanimoto', 'dpic', 'sali', 'involves_hub']].copy()
out3.columns = ['cpd_A', 'cpd_B', 'tanimoto', 'delta_pIC50', 'SALI', 'hub_associated']
out3['literature_source'] = '[author to curate]'
out3['agreement_status'] = '[author to curate]'
out3.to_csv(OUT / 'representative_cliffs_for_lit.csv', index=False)
print(f"  Selected {len(out3)} representative cliffs (hub + non-hub, high SALI/ΔpIC50).")
print(f"  Literature columns left for manual author curation (no fabrication).")
print(f"  Saved: {OUT/'representative_cliffs_for_lit.csv'}")
print("\nDONE")
