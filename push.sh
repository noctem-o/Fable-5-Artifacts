#!/usr/bin/env bash
# One-command upload of the One-Exponent Program artifacts to your repo.
# Run this from a machine where your GitHub auth lives (gh, SSH, or a cached HTTPS token).
set -euo pipefail
REPO="https://github.com/noctem-o/Fable-5-Artifacts.git"
TMP="$(mktemp -d)"
git clone "$REPO" "$TMP"
# copy everything sitting next to this script (except the script itself) into the repo
cp "$(dirname "$0")"/*.md "$(dirname "$0")"/*.py "$(dirname "$0")"/*.png "$TMP"/ 2>/dev/null || true
cd "$TMP"
git add -A
git commit -m "Add One-Exponent Program: manuscript, theorems, proofs, pipeline, falsification record"
git push origin main
echo "Done — pushed to $REPO"
