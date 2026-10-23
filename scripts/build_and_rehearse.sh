#!/usr/bin/env bash
# Build the nsforest image from the committed code, check it, and rehearse the whole workflow on a small file,
# with the real images. Nothing is pushed unless you say so. Everything is written to rehearsal_run_<date>/ and the
# verdict is in rehearsal_run_<date>/REPORT.txt.
#
#   scripts/build_and_rehearse.sh /path/to/adata_sampled_AIBS_human_neocortex_plainstrings.h5ad
#   scripts/build_and_rehearse.sh /path/to/...plainstrings.h5ad push      # also push the image when everything passed
#
# Run it from the repository, on the branch you committed. Keep the computer awake: caffeinate -i scripts/build_and_rehearse.sh ...
set -uo pipefail

H5AD="${1:?give the path to the plain-strings h5ad}"
PUSH="${2:-}"
PYTHON="${PYTHON:-python3}"
IMAGE="ghcr.io/nih-nlm/sc-nsforest-qc-nf/nsforest"

cd "$(git rev-parse --show-toplevel)"
SHA="$(git rev-parse HEAD)"
RUN="rehearsal_run_$(date +%Y%m%d_%H%M)"
mkdir -p "$RUN"
REPORT="$RUN/REPORT.txt"
say() { echo "$@" | tee -a "$REPORT"; }
stop() { say ""; say "STOPPED: $*"; say "RESULT: FAIL"; exit 1; }

say "rehearsal started $(date)"
say "commit $SHA   branch $(git branch --show-current)"
if [ -n "$(git status --porcelain -- main.nf nextflow.config modules container params scripts)" ]; then
    say "WARNING: uncommitted changes in the pipeline files; the image is built from the working folder, not from the commit"
fi

# 1. build
say ""; say "[1/5] docker build (amd64; emulated on Apple Silicon, this is the slow part)"
docker build --platform linux/amd64 -t "$IMAGE:$SHA" container/nsforest/ > "$RUN/build.log" 2>&1 || stop "docker build failed, see $RUN/build.log"
say "      built $IMAGE:$SHA"

# 2. check the image
say ""; say "[2/5] checks on the image"
V="$(docker run --rm --platform linux/amd64 "$IMAGE:$SHA" pip show nsforest 2>/dev/null | awk '/^Version:/{print $2}')"
say "      nsforest version: $V"; [ "$V" = "4.0.0" ] || stop "NS-Forest is not 4.0.0"
WHERE="$(docker run --rm --platform linux/amd64 "$IMAGE:$SHA" python -c "import importlib.util as u; print(u.find_spec('nsforest').origin)" 2>/dev/null)"
say "      nsforest loads from: $WHERE"; case "$WHERE" in */site-packages/*) ;; *) stop "nsforest does not load from site-packages";; esac
docker run --rm --platform linux/amd64 "$IMAGE:$SHA" python -c "import importlib.util as u; print(u.find_spec('nsforest_cli.nsf_compat').origin)" > "$RUN/compat.txt" 2>&1 || stop "nsf_compat.py is not in the image"
HELP="$(docker run --rm --platform linux/amd64 "$IMAGE:$SHA" nsforest-cli --help 2>/dev/null)"
for c in combinations compare-published dendrogram prep run-nsforest merge-nsforest-results plots; do
    echo "$HELP" | grep -q "$c" || stop "nsforest-cli has no command $c"
done
say "      commands present, nsf_compat.py present"

# 3. a small file with the real structure
say ""; say "[3/5] small file from $H5AD"
"$PYTHON" scripts/make_rehearsal_subsample.py "$H5AD" "$RUN/rehearsal_subsample.h5ad" --cluster-header CrossArea_cluster --per-cluster 20 2>&1 | grep -v -i warn | tee -a "$REPORT" || stop "the subsample could not be made"
[ -s "$RUN/rehearsal_subsample.h5ad" ] || stop "no subsample file"

# 4. the whole workflow
say ""; say "[4/5] nextflow, both cluster headers, the real images"
nextflow run main.nf -c params/rehearsal.config \
    --datasets_csv params/neocortex_datasets.csv --organ neocortex \
    --custom_h5ad "$PWD/$RUN/rehearsal_subsample.h5ad" \
    --nsforest_image "$IMAGE:$SHA" \
    --outdir "$PWD/$RUN/out" -work-dir "$PWD/$RUN/work" > "$RUN/nextflow.log" 2>&1
say "      nextflow exit status: $?"

# 5. the verdict
say ""; say "[5/5] report"
"$PYTHON" scripts/rehearsal_report.py "$RUN/out" | tee -a "$REPORT"
VERDICT=${PIPESTATUS[0]}

if [ "$VERDICT" = "0" ] && [ "$PUSH" = "push" ]; then
    say ""; say "pushing $IMAGE:$SHA"
    docker push "$IMAGE:$SHA" 2>&1 | tee -a "$REPORT"
    say "use the digest printed above (digest: sha256:...) as --nsforest_image on Lifebit"
fi
say ""; say "finished $(date); everything is in $RUN/"
exit "$VERDICT"
