#!/usr/bin/env bash
#
# Clone or update CombinatorialMultigrid.jl into the laplacian-bench repo so the
# Julia environment develops it locally (setup.jl auto-detects ./CombinatorialMultigrid.jl).
#
# Pure git — no Julia — so it is safe on a cluster LOGIN node that has internet
# but restricts heavy compute (e.g. NJIT Wulver). Run the Julia setup on a
# COMPUTE node afterwards; with the checkout present it resolves fully offline:
#
#   ./fetch_cmg.sh                 # clone/update to main
#   ./fetch_cmg.sh --rev v1.2.3    # a specific branch/tag/commit
#
set -euo pipefail

REV="main"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --rev) REV="$2"; shift 2 ;;
        -h|--help) sed -n '2,14p' "$0"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 1 ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"          # laplacian-bench root
DEST="$ROOT/CombinatorialMultigrid.jl"
URL="https://github.com/ikoutis/CombinatorialMultigrid.jl.git"

if [[ -d "$DEST/.git" ]]; then
    echo "updating $DEST -> $REV"
    git -C "$DEST" fetch origin "$REV"
    git -C "$DEST" checkout "$REV"
    git -C "$DEST" pull --ff-only origin "$REV" || true   # no-op for a detached tag/commit
else
    echo "cloning $URL -> $DEST"
    git clone "$URL" "$DEST"
    git -C "$DEST" checkout "$REV"
fi

echo "CombinatorialMultigrid.jl ready at $DEST ($(git -C "$DEST" rev-parse --short HEAD))"
if [[ ! -f "$DEST/src/elimination.jl" ]]; then
    echo "WARNING: src/elimination.jl not found — this checkout predates the degree-1/2" >&2
    echo "         elimination branch; the cmg-*-elim columns will not build." >&2
fi
