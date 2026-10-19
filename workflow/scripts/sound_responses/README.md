# Bayesian-block sound responses

This optional stage estimates sound-aligned temporal modulation and response
boundaries. It depends only on `meta.h5` and `units.h5`; pose, speaker geometry,
task performance and spatial tuning are not prerequisites. It does not change
existing unit-response metrics, PSTH analyses or downstream unit selection.

Enable `blocks.sound_responses_bayesian_blocks: true` in a user configuration,
or request one of its outputs explicitly. It works for active and passive
sessions. All existing configurations leave this new block disabled by default.

## Environment

Create the dedicated environment using the repository specification (from the
repository root):

```bash
conda env create \
  --prefix /mnt/nevermind.data-share/ag-grothe/AG_Pecka/envs/sound_responses_bayesian_blocks \
  --file workflow/envs/sound_responses_bayesian_blocks.yaml
```

Run Snakemake with `--use-conda`. The configurable
`sound_responses_bayesian_blocks.conda_env` selects this prefix. No installation
into `pipe-dev` is required. The default CPU backend is `hepstats==0.9.2`, whose
Bayesian-block implementation is compared against the older Scikit-HEP backend
used by the reference. Astropy's similarly named implementation is not silently
substituted. Package versions and implementation hashes are saved with results.
The environment sets `PYTHONNOUSERSITE=1` to prevent user-site packages from
overriding its dependencies. Activate it before running the tests below.
The adjacent `.lock.yaml` records the resolved Linux-64 package builds for
recreating this installation exactly; the shorter YAML records direct dependencies.

## Optional single-GPU backend

Set `sound_responses_bayesian_blocks.backend: gpu` to use the CuPy/CUDA
implementation. `cpu` remains the default, and neither backend changes the
pipeline's existing response metrics or unit selection. GPU fitting uses
unbinned timestamps and float64 arithmetic without rounding, subsampling,
fast-math or fused multiply-add. It implements the same unweighted recurrence,
prior and first-index tie handling as the hepstats reference. Both engines feed
the same onset/offset selection function.

Create the separate environment (or use its adjacent Linux-64 lock file):

```bash
conda env create \
  --prefix /mnt/nevermind.data-share/ag-grothe/AG_Pecka/envs/sound_responses_bayesian_blocks_gpu \
  --file workflow/envs/sound_responses_bayesian_blocks_gpu.yaml
```

This leaves `pipe-dev` and the CPU environment unchanged. The GPU environment
pins CuPy 13.6 and CUDA 12.9; a compatible NVIDIA driver is still required on the
host. The ordinary CPU environment remains in use for timing selection and the
diagnostic PDF rule; neither step requires GPU access.

```yaml
sound_responses_bayesian_blocks:
  backend: gpu
  gpu_conda_env: /mnt/nevermind.data-share/ag-grothe/AG_Pecka/envs/sound_responses_bayesian_blocks_gpu
  gpu_memory_mb: 4096
```

Launch Snakemake with **exactly one** GPU exposed through
`CUDA_VISIBLE_DEVICES`. Prefer a UUID from
`nvidia-smi --query-gpu=index,uuid,name --format=csv`; if a scheduler already
provides one visible GPU, retain its assignment. For example, prepend
`CUDA_VISIBLE_DEVICES=GPU-<selected-UUID>` to the normal Snakemake command with
`--use-conda`. There is no automatic device selection or silent CPU fallback.

The fitting rule requests `sound_responses_gpu=1`, with a workflow-wide capacity of one
even when `--resources` is omitted. A host-local, per-user/per-device file lock
also serializes this stage across separate local Snakemake invocations. It does
not reserve the GPU against unrelated programs or other users. Plotting jobs
request no GPU.

Unit/condition cases are batched across the session. `gpu_memory_mb` bounds the
device allocation pool; the CUDA context requires additional memory. Batches
are split conservatively to fit this budget. One case exceeding the budget
fails with an explicit error instead of changing its data or switching engines.
Aligned spikes are retained in host memory until their cache entries are written.
The GPU supports 1–16 configured penalty values; longer lists require CPU.
Numerically invalid cells mark only the affected case unavailable, as in the
CPU path. Other CUDA failures stop the job.

Outputs retain the existing schema, adding `backend` and `backend_provenance`
attributes. These record the engine, device, CuPy/runtime/driver versions,
compiler options, memory budget, batch count and elapsed batch-analysis time.
Implementation hashes include the GPU Python source and CUDA kernel. `fits.h5`
stores every configured penalty level. The evaluated fits in `responses.h5`
still stop at the same level as the reference timing analysis.

Tests without CUDA cover input validation and the CPU path. Activate the GPU
environment and opt in explicitly to run the GPU equivalence tests:

```bash
RUN_PPLSIT_GPU_TESTS=1 CUDA_VISIBLE_DEVICES=GPU-<selected-UUID> \
  python -B -m unittest discover -s tests -p 'test_sound_responses*.py' -v
```

These compare frozen reference cases, custom penalty lists, repeated timestamps,
forced multi-batch execution, numerical failures, QC and session HDF5/CSV outputs.
Exact agreement on these tests and the A100 pilot is validation evidence, not a
promise of bitwise agreement for every possible dataset or future CUDA version.

Before the cache split, the integrated Snakemake pilot
(`017276_SIT_2026-03-04_14-14-44`, 179 units,
537 unit/condition results) completed its HDF5/CSV job in 3 min 23 s on one A100,
including 165.3 s of batch analysis and a 1.82 GiB peak device allocation pool.
With the 240 ms observation window and duration + (0, 20) ms offset search, all
537 results, 2,430 stored fits and 6,482 datasets matched the CPU reference
exactly. This timing excludes plot generation.

## Example configuration

```yaml
blocks:
  sound_responses_bayesian_blocks: true

sound_responses_bayesian_blocks:
  window_s: null
  include_pooled: true
  conditions: null
  unit_names: null
  time_ranges_s: null
```

The full defaults are in `config/default.yaml`. Null `window_s` uses the session's
`sound.latency`, which is its nominal sound-onset interval, not neural latency.
An explicit window is required if that metadata is absent.
The Miguel 017276 configuration sets this analysis's window to 240 ms to retain
more complete windows despite onset jitter in its nominal 250 ms schedule.
This setting also controls this stage's diagnostic PSTHs; it does not change
other PSTH analyses, stimulus durations or synchronization.

Null `conditions` analyzes each observed positive sound ID separately, plus a
pooled condition. Names, frequencies, durations and channels are retained from
the stimulus catalog, with a fallback to the session parameters for older meta
files. For example, `conditions: {background: [1], target: [2]}` selects those
stimuli; the pooled condition then uses their union. IDs are *sound IDs*, never
speaker/channel numbers. Synchronization channel 8 in a stimulus's channel list
does not become a separate condition. Silence/noise IDs 0 and -1 are excluded.

Optional `time_ranges_s: [[10, 300], [600, 1200]]` restricts analysis to complete
windows inside disjoint, sorted session-time intervals. Null uses the timeline's
recording bounds. `unit_names` can restrict a pilot run without modifying units.h5.

## What the estimates mean

- `is_modulated`: -1 means unavailable, 0 means no accepted within-window change,
  and 1 means an accepted change in firing rate after alignment to sound onset.
  This is not a generic auditory-responsive label or a comparison to silence.
  Sustained sound responses, for example, need not have an internal rate change.
- `onset_latency_s` and `offset_latency_s` are both measured from the sound onset
  stored in `meta.h5`. The latter is a neural-response boundary, not a measurement
  of the physical stimulus offset. Increases and decreases in firing can qualify;
  the method does not classify excitation versus suppression.
- `modulation_p0` is the last successful segmentation setting. `onset_edge_p0`
  records the level supplying the selected latency, which can differ.
  `offset_p0` records the selected offset level. These are regularization settings,
  not calibrated p-values for a responsiveness test, and should not be used as
  replacements for the pipeline's existing significance measures.
- Missing timing is NaN with an explicit status. Insufficient events/spikes or an
  exceeded computation guard are never encoded as a negative modulation result.

## Reference and explicit adaptations

The reference is `code_07062021/neurofunctions/general.py`, functions
`determine_latency_bayesianblocks` and `determine_offset`. Its SHA-256 is saved in
the output. This identifies the provided dated code; it does not establish the
authors' exact published software environment.

Preserved: three copies of unbinned aligned spikes at `-window`, zero and
`+window`; the p0 sequence 0.05, 0.01, 0.001, 0.0001, 0.00001; one-pass deletion
of the right edge of blocks shorter than 1 ms; boundaries strictly inside
1 ms to `window - 1 ms`; previous-successful-level latency selection; and the
reference's special handling of latencies below 7 ms. Offset fitting reuses the
same segmentations rather than recomputing them.

Adaptations are deliberate and recorded in settings/statuses:

- The original fixes the cycle at 250 ms. Here the observation window is
  configurable, and spike intervals are half-open `[onset, onset + window)`.
- Only complete, equally exposed windows are fitted. Windows extending beyond
  the recording/time ranges or into the next logged event are excluded. This
  can remove many events when actual intervals jitter around the nominal period.
  Inspect `n_events_used` and the exclusion flags; an explicitly shorter window
  (for example 240 ms for a 250 ms schedule) may retain more events. No variable
  exposure correction, padding, jitter tolerance or silent subsampling is applied.
  Missing/wrong event times in the input are not repaired by this stage.
- The reference's offset function hardcodes a search 50–70 ms after neural
  onset for 50 ms sounds, ignoring its duration argument. Here the default
  bounds are neural onset + actual sound duration + [0, 20] ms, with both
  endpoints excluded as in the reference. This reproduces the reference search
  for 50 ms sounds and adapts it for other durations. The earlier pilot used
  [-7, +13] ms; that remains available through
  `offset_search_relative_to_duration_s: [-0.007, 0.013]`. Saved outputs record
  their effective settings. Pooled conditions with mixed or unknown durations
  have no offset estimate.
- Reference fallback boundaries outside the analysis window are flagged and
  withheld. Offset searches extending outside the window are also unavailable.
- Minimum counts (20 events, 20 spikes) are configurable QC additions, not power
  guarantees. The optional `max_spikes` guard marks a result unavailable instead
  of subsampling. Both backends are quadratic in unique spike times after
  tripling; the GPU accelerates the recurrence and batches independent cases.
  Start with a few units when evaluating large recordings on CPU.
- Backend numerical failures (for example nearly coincident times collapsing
  a cell width to zero) are marked unavailable with a diagnostic. Spike times
  are not silently rounded to work around such a failure.
- PSTHs are descriptive firing-rate plots only. They do not enter Bayesian-block
  fitting and do not reproduce the paper's bootstrap confidence intervals.

These outputs do not establish synchronization quality or control for behavior,
trial state, frequency, sound level or unequal background/target sampling. Those
remain relevant when interpreting between-condition differences and developing
the later spatial-tuning stage.

## Outputs and reuse

All files live under a session's `analysis/sound_responses/bayesian_blocks/`:

- `fits.h5`: persistent segmentation cache, including **all** configured p0
  levels before short-block filtering or onset-based stopping. Each unit has
  raw edges, per-level fitting statuses, counts and compressed aligned spikes.
  The file records fitted inputs, event eligibility, condition/duration metadata,
  fitting settings, backend provenance and a cache schema/completion marker.
- `responses.h5`: method/configuration/provenance attributes, source sound events
  and event exclusion flags; condition metadata, selected event indices/times;
  and a `units/<unit>/` group per condition containing scalar result attributes,
  descriptive PSTHs and raw/filtered edges for each evaluated p0.
- `summary.csv`: one row per condition/unit, with timing, counts and statuses.
- `figures/responses.pdf`: diagnostic PSTHs with sound duration and detected
  response boundaries. Empty inputs produce an explicit placeholder report.
- `figures/psth_bgr_tgt_bar.pdf` (active sessions): the existing micro-PSTH
  BGR/TGT bar layout with matching-color dashed onset and dotted offset lines.
  It reuses `utils.psth.get_spike_counts`, `psth.micro` bins/limits and the
  shank/depth unit order. Bars use all logged BGR/TGT events, as in the original
  plot; timing markers use the saved eligible-window estimates. The header
  reports both event counts, and each panel reports timing or missing statuses.
  A first-row PNG preview is saved alongside it. Original PSTH PDFs are unchanged.
  This plotting rule uses the base pipeline environment, like `psth_micro`,
  and does not rerun the Bayesian-block fit.

Snakemake runs `fit_sound_responses_bayesian_blocks` →
`sound_responses_bayesian_blocks` → plotting rules. The first rule creates
`fits.h5` on the configured CPU/GPU backend. The second reads only that cache,
selects onset/offset boundaries and builds diagnostic PSTHs on CPU. `fits.h5`
is retained rather than marked as a temporary output.

| Change | Work repeated |
| --- | --- |
| `offset_search_relative_to_duration_s`, `short_latency_s`, `edge_margin_s`, `min_block_s` | Timing selection and dependent plots |
| `psth_bin_s` | Selection/PSTH output and dependent plots; no fitting |
| Observation window, p0 values, conditions, event/time ranges, units, fitting QC limits | Fitting, selection and plots |
| Source `meta.h5` or `units.h5`, fitting implementation or backend/execution configuration | Fitting, selection and plots |
| Plot appearance/code only | Affected plots |

Stage-specific Snakemake parameters prevent a timing-only configuration edit
from invalidating the fit cache. Keep Snakemake's default rerun triggers, which
include `params`; **`--rerun-triggers mtime` alone does not detect configuration
changes**. Backend selection applies only to fitting. When the cache is current,
timing changes can be run without `CUDA_VISIBLE_DEVICES` or an available GPU.

Changing fitting QC limits deliberately invalidates the cache because those
limits determine which unit/condition cases were fitted. A backend change also
rebuilds the cache so its provenance remains accurate. Both CPU and GPU cache
all configured levels; this can cost more than the old early-stopping CPU run,
but avoids refitting when changed criteria require previously unused levels.

`responses.h5` records the cache path/SHA-256, original fitting provenance,
selection implementation hashes and effective timing settings. It keeps the
existing result/PSTH layout for downstream consumers. The selector rejects
incomplete caches or incompatible fitting settings and never modifies the cache.
PSTH rebinning uses its stored aligned spikes, so it needs no source spike-file
reads or fitting. A numerical failure at an unused later penalty does not change
an earlier valid result; it becomes unavailable only if selection reaches it.

Older `responses.h5` files are **not** accepted as caches: their early stopping
may have omitted necessary penalty levels. Existing sessions need one fresh
fitting run to create `fits.h5`; subsequent timing changes reuse it.

The pilot cache contains all 2,685 raw fits in a 37.7 MiB file. Selection from
this cache took approximately 3.5 seconds on CPU, excluding plot generation.
Switching between the old (-7, +13) ms and default (0, 20) ms offset searches
reproduced both prior reference outputs exactly for all 537 unit/condition cases,
with the cache hash and modification time unchanged. This was tested through
Snakemake with every GPU hidden: offset changes scheduled only selection/plots;
a changed observation window scheduled fitting as well.

The numerical fit/data preparation code is in `segmentation.py`, `gpu.py` and
`gpu_kernel.cu`; `fits_session.py` writes the cache. `bayesian_blocks.py` selects
boundaries and `session.py` writes response outputs. `settings.py` defines the
stage parameter sets. Thin scripts connect these helpers to Snakemake. No
archive code is imported at runtime. The Python convenience API `run_session`
runs both stages and writes a sibling `<response-stem>.fits.h5`; use
`select_session` directly for cache-only reuse.

Tests include frozen results from the dated reference functions and the legacy
Scikit-HEP backend on flat, excited and suppressed synthetic spike trains. All
five p0 segmentations and the timing selections agree for these cases. Tests
also cover event exclusions, edge cases, unavailable results, HDF5/CSV/PDF I/O,
cache completeness, criteria requiring later penalty levels, rebinning without
fitting, incompatible cache settings and cache immutability.

Run regression tests with the dedicated environment's Python:

```bash
python -m unittest discover -s tests -p 'test_sound_responses*.py' -v
```
