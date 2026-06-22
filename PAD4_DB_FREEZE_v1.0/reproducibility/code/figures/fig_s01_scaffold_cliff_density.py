"""
fig_s01_scaffold_cliff_density.py — Supplementary Figure S8.

Which PAD4 chemotypes produce rugged SAR? Ranks Murcko scaffold series (n>=4)
by within-scaffold severe-cliff count and density, and contrasts equally-sampled
rugged vs smooth series — showing ruggedness is scaffold-intrinsic, not a
sampling artifact.

Outputs:
  publication/figures/supplementary/fig_s01_scaffold_cliff_density.{png,pdf}
  outputs/tables/supp_scaffold_cliff_density.csv   (Scaffold | Compounds | Cliffs | Density)
"""
import sys
sys.path.insert(0, 'publication/scripts/figures')
from figure_style import set_style, panel_label, save_fig, SEM, C, DOUBLE

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

set_style()
ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT  = ROOT / 'publication/figures/supplementary'
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("FIGURE S8 — SCAFFOLD CLIFF-DENSITY RANKING")
print("=" * 60)

assets = pd.read_parquet(ROOT / 'data/interim/shared_assets.parquet')
scaff  = pd.read_csv(ROOT / 'outputs/tables/05_scaffold_summary.csv')
scaff  = scaff.sort_values('n_compounds', ascending=False).reset_index(drop=True)
scaff['rank'] = scaff.index + 1
ac  = pd.read_parquet(ROOT / 'data/processed/activity_cliffs.parquet')
sev = ac[ac.cliff_tier == 'severe'].copy()

ik2s = dict(zip(assets.inchi_key, assets.murcko_smiles))
sev['sa'] = sev.inchi_key_a.map(ik2s); sev['sb'] = sev.inchi_key_b.map(ik2s)
within = sev[sev.sa == sev.sb]
cc = within.groupby('sa').size().rename('n_cliffs')

S = scaff[scaff.n_compounds >= 4].copy()
S = S.merge(cc.reset_index().rename(columns={'sa': 'scaffold_smiles'}),
            on='scaffold_smiles', how='left')
S['n_cliffs'] = S.n_cliffs.fillna(0).astype(int)
S['poss'] = S.n_compounds * (S.n_compounds - 1) / 2
S['density'] = S.n_cliffs / S.poss
S = S.sort_values(['n_cliffs', 'density'], ascending=False).reset_index(drop=True)

# ── Deposit the ranking table requested by reviewer ───────────────────────────
tbl = S[['rank', 'n_compounds', 'n_cliffs', 'density', 'std_pic50', 'mean_pic50',
         'contains_patent_exclusive', 'scaffold_smiles']].copy()
tbl.columns = ['Scaffold rank', 'Compounds', 'Severe cliffs', 'Cliff density',
               'Intra-scaffold sigma', 'Mean pIC50', 'Contains patent-excl', 'Scaffold SMILES']
tbl.to_csv(ROOT / 'outputs/tables/supp_scaffold_cliff_density.csv', index=False)
n_with = int((S.n_cliffs > 0).sum())
n_smooth_large = int(((S.n_compounds >= 10) & (S.n_cliffs == 0)).sum())
print(f"  scaffolds (n>=4): {len(S)} | with >=1 cliff: {n_with} | smooth large (n>=10): {n_smooth_large}")

# ── Figure: 2 panels ──────────────────────────────────────────────────────────
fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(DOUBLE, 3.2),
                                 gridspec_kw={'width_ratios': [1.25, 1.0]},
                                 constrained_layout=True)

# Panel a: the cliff-bearing scaffolds ranked by cliff count
cb = S[S.n_cliffs > 0].copy().sort_values('n_cliffs')
ypos = np.arange(len(cb))
colors = [SEM['patent'] if p else C['blue'] for p in cb.contains_patent_exclusive]
ax_a.barh(ypos, cb.n_cliffs, color=colors, height=0.66, edgecolor='white', lw=0.4)
for y, r in zip(ypos, cb.itertuples()):
    ax_a.text(r.n_cliffs + 0.3, y, f'n={r.n_compounds}, ρ={r.density:.3f}',
              va='center', ha='left', fontsize=5.3, color=C['gray_dark'])
ax_a.set_yticks(ypos)
ax_a.set_yticklabels([f'S{r}' for r in cb['rank']], fontsize=5.5)
ax_a.set_xlabel('Within-scaffold severe cliffs')
ax_a.set_ylabel('Scaffold rank')
ax_a.set_xlim(0, cb.n_cliffs.max() * 1.45)
ax_a.set_title(f'{n_with}/{len(S)} scaffolds (n≥4) harbor any cliff', fontsize=6.5, pad=3)
panel_label(ax_a, 'a', x=-0.16, y=1.04)

# Panel b: rugged vs smooth at matched sampling density (n vs cliffs, colour=density)
ax_b.scatter(S.n_compounds, S.n_cliffs, s=18, c=C['grey'],
             edgecolors='none', zorder=2, label='all series (n≥4)')
rugged = S[S.n_cliffs > 0]
ax_b.scatter(rugged.n_compounds, rugged.n_cliffs, s=26, c=SEM['cliff'],
             edgecolors='white', lw=0.3, zorder=3, label='rugged (≥1 cliff)')
# annotate the matched-density contrast: rank1 vs rank2
r1 = S[S['rank'] == 1].iloc[0]; r2 = S[S['rank'] == 2].iloc[0]
ax_b.annotate(f'S1: n={int(r1.n_compounds)}, {int(r1.n_cliffs)} cliffs (rugged)',
              xy=(r1.n_compounds, r1.n_cliffs), xytext=(r1.n_compounds - 150, r1.n_cliffs - 3),
              fontsize=5.3, color=SEM['cliff'],
              arrowprops=dict(arrowstyle='->', lw=0.5, color=SEM['cliff']))
ax_b.annotate(f'S2: n={int(r2.n_compounds)}, 0 cliffs (smooth, σ={r2.std_pic50:.2f})',
              xy=(r2.n_compounds, 0), xytext=(r2.n_compounds - 60, 6),
              fontsize=5.3, color=C['gray_dark'],
              arrowprops=dict(arrowstyle='->', lw=0.5, color=C['gray_dark']))
ax_b.set_xlabel('Scaffold series size')
ax_b.set_ylabel('Within-scaffold severe cliffs')
ax_b.set_title('Ruggedness is scaffold-intrinsic, not a sampling artifact', fontsize=6.5, pad=3)
ax_b.legend(fontsize=5.5, loc='upper right', framealpha=0.0)
panel_label(ax_b, 'b', x=-0.16, y=1.04)

save_fig(fig, str(OUT / 'fig_s01_scaffold_cliff_density'))
plt.close(fig)
print("Figure S8 complete.")
