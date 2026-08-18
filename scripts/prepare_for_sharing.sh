#!/usr/bin/env bash
# Untrack files that are specific to one lab's datasets, so the repo can be shared without
# publishing unpublished project paths. The files STAY ON DISK — this only stops git tracking them
# and adds them to .gitignore.
#
#   bash scripts/prepare_for_sharing.sh --dry-run   # show what would change (default)
#   bash scripts/prepare_for_sharing.sh --apply
#
# What it removes from version control and why:
#   config/*.yaml except the template  absolute paths into manuscript directories
#   inputs/                            symlinks pointing at those same paths (machine-specific)
#   ref_materials/                     real annotation spreadsheets from unpublished work
#   INITIAL_PLAN*.MD                   planning notes for this project's own development
#   docs/DESIGN_NOTES.md               development history and per-dataset measurements
#
# Kept: app/, scripts/, preprocess/, tests/, run_app.sh, requirements.txt, environment.yml,
# README.md, and config/dataset.template.yaml — everything a new user needs.
set -euo pipefail
cd "$(dirname "$0")/.."
APPLY=false
[[ "${1:-}" == "--apply" ]] && APPLY=true

PATTERNS=(
  'config/*.yaml'
  'inputs'
  'ref_materials'
  'INITIAL_PLAN.MD'
  'INITIAL_PLAN_ANSWERS.MD'
  'docs/DESIGN_NOTES.md'
)

echo "Currently tracked files that would be untracked:"
TOTAL=0
for pat in "${PATTERNS[@]}"; do
    n=$(git ls-files -- "$pat" | grep -v 'dataset.template.yaml' | wc -l | tr -d ' ')
    [[ "$n" == "0" ]] && continue
    TOTAL=$((TOTAL + n))
    printf '  %5s  %s\n' "$n" "$pat"
done
echo "  -----"
printf '  %5s  total\n' "$TOTAL"
echo
echo "Kept (what a new user needs):"
git ls-files | grep -vE '^(inputs/|ref_materials/|INITIAL_PLAN|docs/DESIGN_NOTES)' \
  | grep -vE '^config/.*\.yaml$' | sed 's/^/  /'
echo "  config/dataset.template.yaml"

if ! $APPLY; then
    echo
    echo "Dry run — nothing changed. Re-run with --apply to do it."
    exit 0
fi

for pat in "${PATTERNS[@]}"; do
    for f in $(git ls-files -- "$pat" | grep -v 'dataset.template.yaml'); do
        git rm --cached -q "$f"
    done
done

# Keep them ignored so they don't get re-added by a later `git add -A`.
python3 - <<'PYEOF'
import pathlib
p = pathlib.Path(".gitignore"); s = p.read_text()
block = """
# --- Lab-specific, kept out of the shared repo (see scripts/prepare_for_sharing.sh) ---
# Real dataset configs hold absolute paths into unpublished project directories.
config/*.yaml
!config/dataset.template.yaml
# Symlinks to source data; machine-specific and point at the same paths.
inputs/
# Annotation spreadsheets from unpublished work.
ref_materials/
# This project's own planning notes and development history.
INITIAL_PLAN.MD
INITIAL_PLAN_ANSWERS.MD
docs/DESIGN_NOTES.md
"""
if "prepare_for_sharing" not in s:
    p.write_text(s.rstrip() + "\n" + block)
    print("  .gitignore updated")
PYEOF

echo
echo "Done. Files are still on disk; git no longer tracks them."
echo "Review with:  git status --short"
echo "Then commit:  git commit -m 'Untrack lab-specific datasets and history for sharing'"
