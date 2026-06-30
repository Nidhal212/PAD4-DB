"""
E3 Cliff-Landscape Permutation Test
Null: permute 3,093 pIC50 values (similarity fixed), 10,000 iterations.
Outputs: outputs/audit/E3_permutation_results.json
         publication/figures/supplementary/fig_s06_permutation.{png,pdf}
"""

import json
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

ROOT   = Path('/home/nidhal/PAD4-db_V2')
SEED   = 42
N_PERM = 10_000

print("=" * 65)
print("E3 CLIFF-LANDSCAPE PERMUTATION TEST")
print("=" * 65)

# ── Load data ────────────────────────────────────────────────────────────────
print("\n[1] Loading data ...")
cpd_df = pd.read_parquet(ROOT / 'publication/data/pad4_compounds.parquet')
pairs_df = pd.read_parquet(ROOT / 'publication/data/activity_pairs_with_sali.parquet')

print(f"    Compounds: {len(cpd_df)}")
print(f"    All pairs: {len(pairs_df)}")

# ── Build compound index ─────────────────────────────────────────────────────
ik_to_idx = {ik: i for i, ik in enumerate(cpd_df['inchi_key'].values)}
n_cpd = len(cpd_df)
pic50_arr = cpd_df['pic50_consensus'].values.astype(np.float64)

# ── Filter high-sim pairs (Tanimoto ≥ 0.8) ───────────────────────────────────
print("\n[2] Filtering Tanimoto ≥ 0.8 ...")
hs = pairs_df[pairs_df['tanimoto'] >= 0.8].copy()
print(f"    High-sim pairs: {len(hs)}")
if len(hs) != 12071:
    print(f"    WARNING: expected 12,071, got {len(hs)}")

# Map pairs to integer indices
idx_a = np.array([ik_to_idx[k] for k in hs['inchi_key_a']], dtype=np.int32)
idx_b = np.array([ik_to_idx[k] for k in hs['inchi_key_b']], dtype=np.int32)
n_pairs = len(idx_a)

# ── SANITY CHECK ─────────────────────────────────────────────────────────────
print("\n[3] Sanity check — identity assignment ...")
obs_deltas = np.abs(pic50_arr[idx_a] - pic50_arr[idx_b])
n_obs_cliffs = int((obs_deltas >= 2.0).sum())
print(f"    Observed cliff pairs (|Δ| ≥ 2.0): {n_obs_cliffs}  (expected 94)")
if n_obs_cliffs != 94:
    raise RuntimeError(f"SANITY FAIL: got {n_obs_cliffs} cliff pairs, expected 94")
print("    SANITY PASS ✓")

# ── Verify observed hub concentration (Stat 2) ────────────────────────────────
print("\n[4] Computing observed hub concentration ...")
HUB_IKS = {
    'A1': 'SMADULGDNOCLOP-GISFHXKWSA-N',
    'A2': 'RAVBZQAQTVGKIV-XBPDSQQVSA-N',
    'B1': 'UDCDEKJNAMHBFH-HSZRJFAPSA-N',
    'B2': 'DVCKJOQIVOGXEI-XMMPIXPASA-N',
}
hub_idx_set = {ik_to_idx[v] for v in HUB_IKS.values()}

cliff_mask = obs_deltas >= 2.0
cliff_a_obs = idx_a[cliff_mask]
cliff_b_obs = idx_b[cliff_mask]
assert len(cliff_a_obs) == 94

# Edges incident to any hub
hub_incident_obs = int(np.sum(
    np.isin(cliff_a_obs, list(hub_idx_set)) |
    np.isin(cliff_b_obs, list(hub_idx_set))
))
obs_hub_pct = hub_incident_obs / 94
print(f"    Hub-incident edges (obs): {hub_incident_obs} / 94 = {obs_hub_pct:.4f}")
if hub_incident_obs != 50:
    print(f"    WARNING: expected 50, got {hub_incident_obs}")
else:
    print("    Hub concentration MATCHES expected 50/94 ✓")

# ── NULL DISTRIBUTION — 10,000 permutations ──────────────────────────────────
print(f"\n[5] Running {N_PERM:,} permutations (seed={SEED}) ...")
rng = np.random.default_rng(SEED)

null_n_cliffs = np.empty(N_PERM, dtype=np.int32)
null_hub_conc = []          # only for perms where n_cliffs ≥ 94
null_hub_conc_all = []      # (permutation index, value) for counting

t0 = time.time()
for perm_i in range(N_PERM):
    perm = rng.permutation(n_cpd)
    perm_pic50 = pic50_arr[perm]

    deltas_p = np.abs(perm_pic50[idx_a] - perm_pic50[idx_b])
    cliff_flags = deltas_p >= 2.0
    nc = int(cliff_flags.sum())
    null_n_cliffs[perm_i] = nc

    # Stat 2: only when enough cliffs to subsample 94
    if nc >= 94:
        cliff_indices = np.where(cliff_flags)[0]
        sub_idx = rng.choice(cliff_indices, size=94, replace=False)
        ca = idx_a[sub_idx]
        cb = idx_b[sub_idx]

        # 4 highest-degree nodes in this null sub-graph
        degree = np.bincount(np.concatenate([ca, cb]), minlength=n_cpd)
        top4_nodes = np.argpartition(degree, -4)[-4:]
        top4_set = set(top4_nodes.tolist())

        incident = int(np.sum(
            np.isin(ca, list(top4_set)) | np.isin(cb, list(top4_set))
        ))
        conc = incident / 94
        null_hub_conc.append(conc)
        null_hub_conc_all.append((perm_i, conc))

    if (perm_i + 1) % 1000 == 0:
        elapsed = time.time() - t0
        print(f"    {perm_i+1:>5,} / {N_PERM:,}  ({elapsed:.0f}s elapsed, "
              f"{elapsed/(perm_i+1)*N_PERM:.0f}s estimated total)")

elapsed_total = time.time() - t0
print(f"    Done in {elapsed_total:.1f}s")

# ── STATISTICS ───────────────────────────────────────────────────────────────
print("\n[6] Computing statistics ...")

null_n_cliffs_arr = null_n_cliffs.astype(float)
null_mean_nc  = float(null_n_cliffs_arr.mean())
null_sd_nc    = float(null_n_cliffs_arr.std())
depletion_ratio = 94 / null_mean_nc
# one-sided depletion p: fraction of perms with n_cliffs ≤ 94
p_depletion = float((null_n_cliffs <= 94).sum()) / N_PERM

null_hub_arr = np.array(null_hub_conc, dtype=float)
n_hub_perms = len(null_hub_arr)
if n_hub_perms > 0:
    null_mean_hub = float(null_hub_arr.mean())
    null_sd_hub   = float(null_hub_arr.std())
    # one-sided enrichment p: fraction with hub_conc ≥ observed
    p_hub_enrichment = float((null_hub_arr >= obs_hub_pct).sum()) / n_hub_perms
else:
    null_mean_hub = null_sd_hub = p_hub_enrichment = float('nan')

print(f"\n    STAT 1 — n_cliffs:")
print(f"      Observed: {n_obs_cliffs}")
print(f"      Null mean ± SD: {null_mean_nc:.2f} ± {null_sd_nc:.2f}")
print(f"      Depletion ratio: {depletion_ratio:.4f}")
print(f"      One-sided p (depletion): {p_depletion:.4f}")

print(f"\n    STAT 2 — hub concentration:")
print(f"      Observed: {hub_incident_obs}/94 = {obs_hub_pct:.4f} ({obs_hub_pct*100:.1f}%)")
print(f"      N perms with n_cliffs ≥ 94: {n_hub_perms}")
if n_hub_perms > 0:
    print(f"      Null mean ± SD: {null_mean_hub:.4f} ± {null_sd_hub:.4f}")
    print(f"      One-sided p (enrichment): {p_hub_enrichment:.4f}")

# ── Save JSON ─────────────────────────────────────────────────────────────────
print("\n[7] Saving results ...")
results = {
    "description": "E3 cliff-landscape permutation test — PAD4-DB v2",
    "seed": SEED,
    "n_permutations": N_PERM,
    "high_sim_pairs": int(n_pairs),
    "n_compounds": int(n_cpd),
    "sanity_check_obs_cliffs": int(n_obs_cliffs),
    "sanity_pass": (n_obs_cliffs == 94),
    "stat1_n_cliffs": {
        "observed": int(n_obs_cliffs),
        "null_mean": round(null_mean_nc, 4),
        "null_sd":   round(null_sd_nc, 4),
        "depletion_ratio": round(depletion_ratio, 4),
        "p_one_sided_depletion": round(p_depletion, 5),
        "n_perms_leq_94": int((null_n_cliffs <= 94).sum()),
        "null_min": int(null_n_cliffs.min()),
        "null_max": int(null_n_cliffs.max()),
        "null_p5":  float(np.percentile(null_n_cliffs, 5)),
        "null_p95": float(np.percentile(null_n_cliffs, 95)),
    },
    "stat2_hub_concentration": {
        "observed_hub_incident": int(hub_incident_obs),
        "observed_total":        94,
        "observed_fraction":     round(float(obs_hub_pct), 5),
        "n_perms_with_geq94_cliffs": int(n_hub_perms),
        "null_mean_hub_conc":    round(float(null_mean_hub), 5) if n_hub_perms else None,
        "null_sd_hub_conc":      round(float(null_sd_hub), 5)   if n_hub_perms else None,
        "p_one_sided_enrichment": round(float(p_hub_enrichment), 5) if n_hub_perms else None,
        "n_perms_geq_observed": int((null_hub_arr >= obs_hub_pct).sum()) if n_hub_perms else None,
    },
    "hub_compounds": HUB_IKS,
    "runtime_seconds": round(elapsed_total, 1),
}

OUT_AUDIT = ROOT / 'outputs/audit'
OUT_AUDIT.mkdir(parents=True, exist_ok=True)
json_path = OUT_AUDIT / 'E3_permutation_results.json'
with open(json_path, 'w') as f:
    json.dump(results, f, indent=2)
print(f"    Saved: {json_path}")

# ── Figure ────────────────────────────────────────────────────────────────────
print("\n[8] Generating figure ...")

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Liberation Sans', 'DejaVu Sans', 'Arial'],
    'axes.linewidth': 0.75,
    'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'xtick.major.size': 3,    'ytick.major.size': 3,
    'axes.labelsize': 7,
    'xtick.labelsize': 6,     'ytick.labelsize': 6,
    'legend.fontsize': 6,
    'figure.dpi': 300,
})

BLUE = '#0077BB'
RED  = '#CC3311'
GRAY = '#BBBBBB'
DARK = '#555555'

def fmt_p(p):
    """Format p-value: show 'p < 0.0001' when p rounds to 0.0000."""
    return 'p < 0.0001' if p < 0.0001 else f'p = {p:.4f}'

fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(7.205, 2.8))  # 183 mm double-col

# ── Panel a: null n_cliffs distribution ──────────────────────────────────────
ax_a.hist(null_n_cliffs, bins=50, color=GRAY, edgecolor='white', lw=0.3,
          density=False, label=f'Null ({N_PERM:,} perms)')
ax_a.axvline(94, color=RED, lw=1.2, ls='-',
             label=f'Observed = 94\n({fmt_p(p_depletion)})')
ax_a.set_xlabel('Number of cliff pairs (|ΔpIC50| ≥ 2.0)', fontsize=7)
ax_a.set_ylabel('Permutations', fontsize=7)
ax_a.text(0.97, 0.95,
          f'Null mean = {null_mean_nc:.1f} ± {null_sd_nc:.1f}\n'
          f'Observed = 94\n'
          f'Depletion ratio = {depletion_ratio:.3f}\n'
          f'{fmt_p(p_depletion)}',
          transform=ax_a.transAxes, fontsize=5.5, va='top', ha='right',
          color=DARK)
ax_a.legend(fontsize=5.5, frameon=False, loc='upper left')
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)
ax_a.tick_params(direction='in', length=3, width=0.5)
ax_a.text(-0.12, 1.04, 'a', transform=ax_a.transAxes, fontsize=8, fontweight='bold')

# ── Panel b: null hub-concentration distribution ──────────────────────────────
if n_hub_perms > 30:
    ax_b.hist(null_hub_arr * 100, bins=40, color=BLUE, alpha=0.7,
              edgecolor='white', lw=0.3,
              label=f'Null ({n_hub_perms} perms with ≥94 cliffs)')
    ax_b.axvline(obs_hub_pct * 100, color=RED, lw=1.2, ls='-',
                 label=f'Observed = {obs_hub_pct*100:.1f}%\n'
                       f'({fmt_p(p_hub_enrichment)})')
    ax_b.set_xlabel('Hub concentration (% edges incident to top-4 nodes)', fontsize=7)
    ax_b.set_ylabel('Permutations', fontsize=7)
    ax_b.text(0.97, 0.95,
              f'Null mean = {null_mean_hub*100:.1f}% ± {null_sd_hub*100:.1f}%\n'
              f'Observed = {obs_hub_pct*100:.1f}%\n'
              f'{fmt_p(p_hub_enrichment)}',
              transform=ax_b.transAxes, fontsize=5.5, va='top', ha='right',
              color=DARK)
    ax_b.legend(fontsize=5.5, frameon=False, loc='upper left')
else:
    ax_b.text(0.5, 0.5, f'Insufficient perms\nwith ≥94 cliffs\n(n={n_hub_perms})',
              transform=ax_b.transAxes, ha='center', va='center', fontsize=7)
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.tick_params(direction='in', length=3, width=0.5)
ax_b.text(-0.12, 1.04, 'b', transform=ax_b.transAxes, fontsize=8, fontweight='bold')

fig.tight_layout(pad=0.8)

OUT_SUPP = ROOT / 'publication/figures/supplementary'
OUT_SUPP.mkdir(parents=True, exist_ok=True)
for ext in ['png', 'pdf']:
    out = OUT_SUPP / f'fig_s06_permutation.{ext}'
    fig.savefig(out, dpi=600, bbox_inches='tight', facecolor='white')
    print(f"    Saved: {out}")

plt.close(fig)

print("\n" + "=" * 65)
print("E3 COMPLETE")
print("=" * 65)
print(f"\nStat 1 (n_cliffs):        obs={n_obs_cliffs}  null={null_mean_nc:.1f}±{null_sd_nc:.1f}  "
      f"ratio={depletion_ratio:.3f}  p={p_depletion:.4f}")
if n_hub_perms:
    print(f"Stat 2 (hub conc):        obs={obs_hub_pct*100:.1f}%  "
          f"null={null_mean_hub*100:.1f}%±{null_sd_hub*100:.1f}%  "
          f"p={p_hub_enrichment:.4f}  (based on {n_hub_perms} eligible perms)")
