# Cluster Annotation Browser

A web tool for annotating single-cell / single-nucleus cluster identities. It puts all the
evidence for one cluster on screen — marker genes, marker motifs, spatial projections, UMAP
highlight, feature plots — asks Claude for a cited cell-type hypothesis, and records your decision
in a form that exports to a spreadsheet.

Runs on a server, viewed in your own browser over an SSH tunnel. One dataset at a time,
configured by a single YAML file.

---

## 1. Install

You need **Python 3.11** (via conda) and **R** (for `Rscript`, used once per dataset to export the
marker tables).

```bash
git clone <this-repo> cluster_annotation
cd cluster_annotation

conda create -y -n cluster_annotation python=3.11
~/.conda/envs/cluster_annotation/bin/pip install -r requirements.txt
```

Or, equivalently, `conda env create -f environment.yml`.

The launcher calls the environment's interpreter directly, so you never need `conda activate`.
If your env lives elsewhere, set `CLUSTER_ANNOTATION_ENV=/path/to/env`.

### Add an API key (only needed for the AI annotation features)

Generate a key at <https://console.claude.com>. A dedicated key with its own spend limit is
worth it — you can then see exactly what this tool costs. Store it outside the repo, readable
only by you:

```bash
mkdir -p ~/.config/cluster_annotation && chmod 700 ~/.config/cluster_annotation
printf 'ANTHROPIC_API_KEY=sk-ant-...\n' > ~/.config/cluster_annotation/secrets.env
chmod 600 ~/.config/cluster_annotation/secrets.env
```

Without a key everything except the AI sections still works.

### Check the install

```bash
./run_app.sh doctor
```

This verifies the interpreter, the Python packages, `Rscript`, your config's input paths, and
that the secrets file exists with owner-only permissions. Fix anything it marks `✗` before going
further.

---

## 2. Point it at your dataset

Copy the template and edit the paths:

```bash
cp config/dataset.template.yaml config/mydataset.yaml
```

You need five things:

| Field | What it is |
|---|---|
| `name` | Short identifier. Becomes the cache directory and the export filename — **must be unique per dataset**, or two datasets will overwrite each other's results. |
| `markers_rds` | `FindAllMarkers()` output for genes, saved as `.rds` |
| `motif_markers_rds` | The same for motifs (uses `avg_diff` instead of `avg_log2FC`) |
| `tangram_pdf` | Optional. Spatial projection PDF, one page per cluster |
| `*_featureplot_glob` | Per-cluster feature-plot PDFs, with `{cluster}` where the cluster id goes |
| `biological_context` | **The most important field.** Describe the sample in prose — species, developmental stage, tissue, assay, caveats. Everything the AI produces is built on this, so be specific. |

The marker `.rds` files are expected to be a length-1 list wrapping a `FindAllMarkers` data frame
(the usual shape when markers are saved per-object). Gene markers may carry the feature in
`rownames`; motif markers in a `gene` column. Both are handled.

Feature-plot PDFs are assumed to be page 1 = UMAP highlight, then one page per marker in the order
your plotting script emitted them. The `featureplot:` block in the config describes that ordering
(how many top markers, sorted by what, with what separator) so page numbers can be matched to
feature names. **If your plotting script used a different order, edit that block to match** — this
is the one part of the config that has to mirror code you wrote elsewhere.

Then verify before launching:

```bash
./run_app.sh doctor config/mydataset.yaml
```

---

## 3. Launch

```bash
./run_app.sh start
```

It prints the exact SSH command to run on your laptop:

```
started (pid 12345, port 8501), logging to logs/app.log

  On your laptop, open the tunnel then browse:
    ssh -N -L 8501:localhost:8501 you@yourserver
    http://localhost:8501
```

The app runs detached, so it survives you logging out.

```bash
./run_app.sh status     # running? which port? tunnel command again
./run_app.sh logs       # follow the log
./run_app.sh stop
./run_app.sh restart
```

If port 8501 is busy the next free port is used and reported — you can run alongside a colleague
on the same server.

The first load of a new dataset runs `Rscript` to export the marker tables to `.cache/<name>/`.
That takes a minute; subsequent loads are instant. Use **Re-preprocess** if the source `.rds`
files change.

---

## 4. Annotate

Pick a cluster in the sidebar, then work through the sections:

| Section | Use |
|---|---|
| **Annotation** | Where you record the decision, and where you export from |
| **All clusters** | Progress across the whole dataset, plus search over the AI annotations |
| **Marker genes / motifs** | Sortable marker tables. Tick **★** on the markers that convinced you — these become the `key_marker_genes` / `select_marker_motifs` columns in the export |
| **Spatial / UMAP / Feature plots** | The images. Feature plots step with **◀ ▶** or the **← →** arrow keys |
| **Other annotations** | Any extra per-cluster PDFs you configured (label transfers, QC overlays) |
| **AI insights** | Claude's cited hypothesis for this cluster |
| **Cohort review** | A whole-dataset critique: likely over-split clusters, inconsistencies, expected cell types that are missing |
| **Comments** | Free-text notes, kept per cluster |

The cluster picker shows state as you go — `12 ✓ Forebrain progenitors 🚩 ★4 💬2` — so you can see
what's done, what's flagged, and what's untouched. **Next unreviewed** jumps to the next gap.

Your notes, stars and annotations live in `annotations/<name>.sqlite`. That file is your work —
**back it up**; it is not regenerable.

---

## 5. AI annotation

Nothing costs money until you click a button, and every paid action asks you to confirm.

**Run it in the app:** open **AI insights** and use *Build reference primer*, then
*Regenerate ALL originals*. Or from the command line, which is better for a whole dataset:

```bash
~/.conda/envs/cluster_annotation/bin/python scripts/generate_all.py config/mydataset.yaml 8
```

(The last argument is how many clusters to process at once.) This does three things:

1. **Builds a reference primer** — one literature-search call that produces a reference sheet for
   your dataset: the cell types plausibly present, their discriminating markers, and a citation
   library where every entry has a PubMed ID or DOI.
2. **Annotates every cluster** against that sheet. These calls cite *by key* into the library, and
   the keys are resolved to real references in code — so a citation that isn't in the library
   cannot appear in your output.
3. **Runs a cohort review** over the finished set and flags cross-cluster problems.

Then, optionally, escalate the clusters where the identity is genuinely uncertain to their own
full literature search:

```bash
python scripts/reannotate_flagged.py config/mydataset.yaml 5 --clusters 3,17,22
```

This keeps the original annotation and writes a separate revised one, so you can compare. Worth
doing for low-confidence or composite-looking clusters. **Not** worth doing for clusters flagged
only as redundant with a neighbour — that's a question about your clustering, and more reading
won't answer it.

### Cost

Every artifact records its own token usage and an estimated cost under `_meta.usage`. You'll see
it per cluster in the app, as a dataset total at the top of **AI insights**, and in the script
output. Expect the primer to be the largest single call and each cluster to be a small fraction of
it. The estimate uses published list prices — check the Console for actual billing.

Knobs in the config, if you need them:

- `effort` — `low`/`medium`/`high`; lower is cheaper and shallower
- `research_mode` — `primer` (default) or `per_cluster`, which gives every cluster its own
  literature search: much slower and dearer, occasionally better on an unusual dataset
- `num_references` — how many citations per cluster
- `primary_model` / `fallback_model` — the fallback is used automatically if the primary declines

Regenerating overwrites an annotation, but the previous version is always copied to
`.cache/<name>/ai_insights/_backups/` first.

---

## 6. Export

In the **Annotation** section, **Prefill from AI annotation** copies the AI's cell-type call and
its PMIDs into the form — then edit it into your own words. The AI output is a first draft, not
the answer.

**Build xlsx export** writes a spreadsheet to `annotations/`. Columns come from three places:

- the annotation form (labels, origin, figure order, references, cell counts)
- your **★** starred markers (`key_marker_genes`, `select_marker_motifs`)
- your sidebar notes (`comments`)

Rows come out in the `annot_order` you assigned, with unnumbered clusters last. By default only
clusters you've labelled are included; toggle that to see what's outstanding.

---

## Troubleshooting

| Symptom | Cause |
|---|---|
| `./run_app.sh doctor` reports a missing import | Env not built, or `CLUSTER_ANNOTATION_ENV` points at the wrong place |
| "Rscript not on PATH" | R isn't installed or loaded. Needed once per dataset to export markers |
| Load fails on a new dataset | Run `doctor` with your config — usually a wrong path or a `{cluster}` token missing from a featureplot glob |
| Feature plot labels don't match the images | The `featureplot:` block doesn't match the ordering your plotting script used |
| AI sections say no API key | Check the secrets file path and that it contains `ANTHROPIC_API_KEY=` |
| `APITimeoutError` during generation | A research call is exceeding the HTTP timeout. Confirm nothing has changed the research call away from streaming |
| Two datasets overwriting each other | They share a `name:` — it drives the cache directory |
| Everything reads as needing regeneration after an upgrade | The prompt version changed, so old annotations aren't comparable to new ones. They still display; only the "needs regeneration" flag changed |

Run the offline test suites any time — they need no API key and no network:

```bash
for t in tests/test_*.py; do ~/.conda/envs/cluster_annotation/bin/python $t; done
```

---

## Layout

```
config/dataset.template.yaml   copy this per dataset
app/                           app.py (UI), data.py, pdf.py, insights.py (AI), store.py, export.py
preprocess/export_markers.R    RDS -> marker tables + feature-plot page index
scripts/generate_all.py        build primer + annotate every cluster + cohort review
scripts/reannotate_flagged.py  escalate specific clusters to their own literature search
scripts/doctor.py              preflight checks (./run_app.sh doctor)
tests/                         offline tests, no API key needed
.cache/<name>/                 derived data + AI annotations (regenerable)
annotations/<name>.sqlite       YOUR notes, stars and annotations — back this up
```
