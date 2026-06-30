"""
supp_s05_scaffold_structures.py — S5: Top 20 Scaffold 2D Structures
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Liberation Sans', 'DejaVu Sans', 'Arial'],
    'figure.dpi': 600,
    'savefig.dpi': 600,
    'savefig.bbox': 'tight',
})

ROOT = Path('/home/nidhal/PAD4-db_V2')
OUT = ROOT / 'outputs/figures/nature_v2'
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("SUPPLEMENTARY S5 — TOP SCAFFOLD 2D STRUCTURES")
print("=" * 60)

scaffold_sum = pd.read_csv(ROOT / 'outputs/tables/05_scaffold_summary.csv')
scaffold_sum = scaffold_sum.sort_values('n_compounds', ascending=False).reset_index(drop=True)
scaffold_sum['scaffold_rank'] = scaffold_sum.index + 1
top20 = scaffold_sum.head(20)

try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem.Draw import rdMolDraw2D
    from PIL import Image
    import io

    mols = []
    legends = []
    for _, row in top20.iterrows():
        mol = Chem.MolFromSmiles(row['scaffold_smiles'])
        mols.append(mol)
        mean_pic50 = row.get('mean_pic50', float('nan'))
        legends.append(f"Rank {int(row['scaffold_rank'])} | n={int(row['n_compounds'])} | pIC50={mean_pic50:.2f}")

    img = Draw.MolsToGridImage(
        mols, molsPerRow=4, subImgSize=(350, 280),
        legends=legends,
        returnPNG=False
    )

    # Save PIL image
    outpath_png = OUT / 'supp_s05_scaffold_structures.png'
    img.save(str(outpath_png), dpi=(600, 600))
    print(f"  Saved: {outpath_png}")

    # Also save as PDF via matplotlib
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.imshow(np.array(img))
    ax.axis('off')
    ax.set_title('Top 20 PAD4 Inhibitor Scaffold 2D Structures',
                  fontsize=12, fontweight='bold', pad=10)
    outpath_pdf = OUT / 'supp_s05_scaffold_structures.pdf'
    fig.savefig(str(outpath_pdf), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {outpath_pdf}")

except Exception as e:
    print(f"  RDKit Draw failed: {e}")
    print("  Creating placeholder figure ...")

    fig, axes = plt.subplots(5, 4, figsize=(8.5, 11))
    for i, (ax, (_, row)) in enumerate(zip(axes.flatten(), top20.iterrows())):
        mean_pic50 = row.get('mean_pic50', float('nan'))
        ax.text(0.5, 0.6, f"Rank {int(row['scaffold_rank'])}",
                transform=ax.transAxes, ha='center', va='center',
                fontsize=10, fontweight='bold')
        ax.text(0.5, 0.35, f"n={int(row['n_compounds'])}\npIC50={mean_pic50:.2f}",
                transform=ax.transAxes, ha='center', va='center', fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)

    fig.suptitle('Top 20 PAD4 Inhibitor Scaffold 2D Structures',
                  fontsize=12, fontweight='bold')
    plt.tight_layout()

    for ext in ['png', 'pdf']:
        outpath = OUT / f'supp_s05_scaffold_structures.{ext}'
        fig.savefig(outpath, dpi=600 if ext == 'png' else 300, bbox_inches='tight')
        print(f"  Saved: {outpath}")
    plt.close(fig)

print("S5 complete.")
