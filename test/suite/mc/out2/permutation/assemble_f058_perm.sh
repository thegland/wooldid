#!/bin/zsh
# Assemble review_process/research_2026-09-06/F058_permutation.md from the prose skeleton,
# the analyzer's tables and the conclusions block.
set -u
REPO=/Users/hegland/Documents/Research/wooldid_dev
MC=$REPO/test/suite/mc
OUT=$REPO/review_process/research_2026-09-06/F058_permutation.md
mkdir -p "$(dirname "$OUT")"
{
  cat "$MC/out2/permutation/F058_perm_prose.md"
  echo
  python3 "$MC/analyze_mcri2.py" "$MC/out2/permutation" | tail -n +2 | python3 "$MC/out2/permutation/renumber.py"
  echo
  cat "$MC/out2/permutation/F058_perm_conclusions.md"
} > "$OUT"
echo "wrote $OUT ($(wc -l < "$OUT") lines)"
