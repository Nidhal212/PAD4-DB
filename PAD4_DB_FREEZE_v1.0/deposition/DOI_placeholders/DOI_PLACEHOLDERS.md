# DOI Placeholders — PAD4-DB v1.0

These DOIs are required before final deposition and submission.

## Zenodo DOI
- Status: NOT YET ASSIGNED
- Action: Pre-reserve DOI at zenodo.org before deposition
- Used in: manuscript reference [27] (currently URL-only: https://www.rdkit.org)
- Manuscript line: ~163 ("RDKit version 2025.09.5 [27]")

## Repository DOI (GitHub via Zenodo)
- Status: NOT YET ASSIGNED
- Action: Tag v1.0 release on GitHub → Zenodo GitHub integration assigns DOI
- Used in: Methods data availability statement

## Instructions
1. Log into zenodo.org
2. New Upload → Reserve DOI (before uploading)
3. Note the DOI
4. Update manuscript ref [27]: replace URL with: Landrum, G. RDKit: Open-Source Cheminformatics Software. https://doi.org/10.5281/zenodo.XXXXXXX (2025).
5. Rebuild manuscript: conda activate pad4bench && python reproducibility/code/build_manuscript_from_md.py
6. Upload deposition/zenodo/ contents to Zenodo
7. Publish — DOI becomes active
