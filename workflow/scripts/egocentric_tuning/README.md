# Neural egocentric tuning

Optional `blocks.egocentric_tuning` stage. Uses `geometry.h5`, Bayesian-block
`responses.h5`, and the **same units export used to compute those responses**.
The response input is `<session>/analysis/bayesian_blocks/responses.h5`.
It does not introduce an auditory responsiveness test or alter existing analyses.
Existing `pipe-dev` NumPy, h5py and matplotlib suffice; no environment changes.

```yaml
blocks:
  egocentric_tuning: True
egocentric_tuning:
  timing_condition: pooled
  conditions: [sound_1, sound_2]
  references: [head, body]  # or moseq_heading if included in geometry.h5
  unit_names: null
  comparisons: [[sound_1, sound_2]]
  n_bootstrap: 1000
  seed: 0
```

The stage is off by default. Missing geometry/response inputs become upstream
Snakemake dependencies even if their blocks are not selected as final targets.
Thus enabling tuning can schedule upstream tracking or Bayesian-block fitting
when their outputs are absent or stale. The tuning and plotting rules themselves
use one CPU each and never run inference or GPU computation.

`conditions: null` uses all nonpooled conditions from the response artifact;
`references: null` uses all saved geometry references. `comparisons: null` compares
the sole pair when exactly two conditions are selected, otherwise no pair is
chosen implicitly. Supply explicit pairs for more conditions or `[]` for none.
Each pair `[A, B]` reports **B minus A** and must contain disjoint events.
There is no assumption that A/B use different speakers, or that the task is
localization. Sound identity, frequency, amplitude and task role remain available
in the saved stimulus catalog. Passive tasks use the same generic calculation
provided tracking geometry and valid single-source events are available.

## Timing and event selection

By default every condition uses each unit's `pooled` response timing, as in the
reference's pooled-sound latency selection. Set `timing_condition: null` to use
each condition's own timing, or give another saved condition name. Missing
requested timing conditions cause an error; missing unit latencies produce
explicit `timing_unavailable` results.

Let L be detected neural onset, O detected offset, D that condition's sound
duration, and W the response artifact's observation window:

* **Onset:** `[L, O)`, or `[L, L+D)` when offset is unavailable.
* **Offset:** the next D seconds after the onset phase ends.
* **Late:** `[0, L)` plus the time after the offset phase up to W.

“Late” follows the original remainder definition: it includes pre-onset time.
All intervals are half-open. Invalid, zero-duration or out-of-observation phases
are unavailable; no clipping, inferred onset at zero, or guessed sound duration.
Mixed-duration conditions are unavailable and should be split upstream.
Raw rates divide phase counts by phase duration.

Only events eligible for the existing response analysis are used (complete,
nonoverlapping observation windows and any configured response time ranges).
They are joined to geometry by **original sound-event index**, with source event
tables and mappings checked for equality. Invalid mathematical geometry is
excluded per reference. No tracking confidence, native-video gap, speed,
center-of-arena, task-state or catch-trial filter is added here.
The stage checks units size/mtime against response provenance and checks each
unit's full-window spike counts against the saved counts. If a units file has
been re-exported or copied without preserving its timestamp, use the original
export or regenerate responses with that file. This prevents silently mixing
responses and spikes from different synchronization exports.

## Reference method and explicit adaptations

Reference: `code_07062021/neurofunctions/general.py`, particularly
`calc_CI_centered`, `calc_vs_centered`, `calc_vectorsum`,
`check_nospikes_asymmetrypresentation`, `compare_vs_checkup_table`, and
`get_spikes_segments`; `calc_psth.py::polar_plot` supplies the phase definitions.
This is the final centered-bootstrap method, not the older 36-bin shuffle code.
Source hashes are saved in `low_spike_thresholds.json` and output provenance.
The lookup values are an exact export of the reference `vs_lookuptable.pkl`.
The archived source files are only read, never changed.

* Eight 45-degree bins centered at 0, 45, ..., 315 degrees; 0 is in front,
  positive is counterclockwise in the calibrated arena. Geometry coverage plots
  have their own configurable bins, which do not affect neural tuning.
* Sum spikes in each bin, then divide by presentations in that bin, including
  zero-spike presentations. Normalize the resulting curve to sum one for vector
  strength (VS) and doubled-angle orientation strength (VS2).
* Resample whole presentations with replacement, 1,000 replicates by default.
  Multinomial sampling of identical `(angular bin, spike count)` pairs has
  exactly the same distribution as drawing event indices. This avoids large
  intermediate spike tables. Seeds are deterministic per unit/condition/reference/
  phase, unaffected by processing order or selecting a subset.
* Save pointwise 95% rate intervals and VS/VS2 15.9th, 50th and 84.1st percentiles.
  Bootstrap unsampled bins are zero-filled for reference metric compatibility;
  observed unsampled bins and their reported rate intervals are NaN.
* Eligibility requires at least 10 presentations in every bin, maximum/minimum
  bin occupancy at most 10, and more than 20 spikes in the phase.
* With at least 75 spikes, the point strength must be >=0.2 and its bootstrap
  15.9th percentile >=0.15. With 21–74 spikes the lower percentile must strictly
  exceed the lookup threshold from the last row with `No_spikes < n_spikes`.
* Choose the larger of VS/VS2 (directional wins ties); classify only if that
  chosen metric passes. Eligibility failures produce `unavailable`, not
  `none`. This is an **Amaro criterion**, not a p-value or auditory
  responsiveness decision.

Adaptations: D is the actual sound duration rather than the plotting code's
hardcoded 57 ms; W comes from the response artifact (240 ms in the pilot rather
than 250 ms). Report Hz and spikes/presentation rather than normalization by the
original attentive-silence baseline, which is not defined generally in pplSIT.
Preferred directions/axes of numerically zero-strength vectors are NaN.
Low-count curves remain visible descriptively but cannot receive a tuning label.

Condition comparisons are an extension: independent bootstrap samples supply
pointwise rate-difference intervals and an equal-angular-weight mean rate
difference interval. Preferred direction/axis shifts and VS/VS2 differences are
descriptive. The mean requires all eight bins to be observed; comparison status
is `descriptive_only` if either curve fails eligibility. Pointwise intervals
are not corrected for multiple bins, units, phases or conditions and do not
constitute an overall test of different tuning. Presentation bootstrap also
does not model temporal dependence between successive sounds. Differences may
reflect frequency, amplitude, movement or task state as well as sound bearing.

## Outputs

Under `<session>/analysis/spatial_tuning/egocentric/`:

* `tuning.h5`: reusable curves, intervals, spike counts and provenance.
* `summary.csv`: one row per unit/condition/reference/phase, including missing
  timing, coverage/spike eligibility, thresholds and classification.
* `comparisons.csv`: one row per pair/unit/reference/phase, scalar differences.
* `figures/tuning.pdf`: one page per unit/reference, with three phase curves and
  condition differences; `tuning.png` previews the first page.

HDF5 root attributes include schema/method, config, source and implementation
hashes, condition/reference/unit order, stimulus catalog, window and conventions.
`complete` is false until calculations and CSV writes succeed.

```
angular_centers_rad
conditions/<condition>/
  event_indices, event_times_s
  geometry/<reference>/bearing_rad, valid
  units/<unit>/<phase>/
    intervals_s
    event_spike_counts              # aligned to condition event_indices
    <reference>/
      presentations, spikes, spikes_per_presentation
      rate_hz, rate_ci95_hz         # CI axis: lower, upper
      vector_strength_interval     # rows: 15.9, 50, 84.1; columns: VS, VS2
comparisons/<pair_index>/<unit>/<reference>/<phase>/
  rate_difference_hz, rate_difference_ci95_hz
  equal_angle_mean_difference_ci95_hz
```

Scalar results/statuses live in group attributes. Unavailable phases retain their
status and empty intervals but have no counts or curve datasets. `preferred_axis_rad`
is modulo pi; `preferred_orientation_rad` is the equivalent axis direction closest
to the maximum-rate bin, matching the reference plotting convention.
Event arrays and saved seeds allow bootstrap reproduction without storing every
replicate. Inputs and existing standard pipeline outputs are never edited.

Validation:

```
/opt/miniconda3/envs/pipe-dev/bin/python -m unittest discover -s tests -p 'test_egocentric_tuning.py' -v
```

Tests include fixed numerical fixtures executed from the archived reference,
exhaustive small-event bootstrap enumeration, occupancy normalization, timing
boundaries, strict low-count thresholds, provenance rejection and an end-to-end
two-condition synthetic calculation and plot.

Pilot validation (2026-10-20): session `017276_SIT_2026-03-04_14-14-44`,
179 units, head/body/corrected MoSeq references, 1,000 bootstrap replicates.
Numerical calculation took 35.5 seconds on one CPU. All 2,916 available curves
were independently reconstructed (maximum VS error 2.3e-16); phase counts summed
to the saved full-window response counts. Seventeen units lacked pooled onset
timing and remained unavailable. Input hashes were unchanged. The report has
537 pages and the completed workflow was up to date under normal rerun triggers.
The geometry and tuning test suites passed all 22 tests.
