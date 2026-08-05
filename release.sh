#!/usr/bin/env zsh
set -euo pipefail

VERSION="${1:-v1.0.0}"
MESSAGE="${2:-PAD4-DB publication release}"

echo "==> Repository root"
git rev-parse --show-toplevel

echo "==> Updating .gitignore"

cat <<'EOF' >> .gitignore

# ---- Publication temporary files ----
*_dryrun.parquet
*_pre_columns_*.parquet
*_pre_remediation_*.parquet
*debug*

__pycache__/
*.pyc
EOF

echo "==> Removing ignored files from index (if any)"
git rm --cached -r --ignore-unmatch \
data/processed/*dryrun*.parquet \
data/processed/*pre_columns*.parquet \
data/processed/*pre_remediation*.parquet || true

echo "==> Staging files"
git add .

echo "==> Git status"
git status --short

echo
read "?Continue with commit? [y/N] " ans

[[ "$ans" =~ ^[Yy]$ ]] || exit 0

echo "==> Commit"
git commit -m "$MESSAGE" || echo "Nothing to commit."

echo "==> Push"
git push origin main

echo "==> Creating tag $VERSION"

git tag -a "$VERSION" -m "$VERSION" 2>/dev/null || \
echo "Tag already exists."

echo "==> Pushing tag"

git push origin "$VERSION"

if command -v gh >/dev/null; then
    echo "==> Creating GitHub Release"

    gh release create "$VERSION" \
        --generate-notes \
        --latest \
    || echo "Release already exists."
else
    echo
    echo "GitHub CLI (gh) not installed."
    echo "Install with:"
    echo "sudo apt install gh"
fi

echo
echo "Done!"
echo
echo "If Zenodo is connected to GitHub,"
echo "a DOI will be minted automatically."
