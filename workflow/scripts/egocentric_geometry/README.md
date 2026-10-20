# Egocentric sound geometry

This optional stage provides animal direction, speaker bearing, distance and
angular event coverage for later spatial-tuning analyses. It uses the existing
`meta.h5`, `DLC_100Hz.h5`, and optionally corrected `MoSeq.h5` outputs.
It runs in the existing pipeline environment and needs no additional packages.
The block is off by default.

Enable it in an animal/project configuration:

```yaml
blocks:
  egocentric_geometry: true

egocentric_geometry:
  references: [head, body]
```

To compare MoSeq's modeled orientation too:

```yaml
egocentric_geometry:
  references: [head, body, moseq_heading]
```

This option requires a `MoSeq.h5` carrying the `arena_circular_v1` heading
marker from the corrected synchronization script. Older files must be
regenerated with `synch_moseq_100Hz`; no inference is necessary.
See [the synchronization notes](../README_moseq_sync.md).
The DLC and MoSeq time grids must match.

## Geometry definitions

- `head`: ear midpoint toward nose.
- `body`: tail base toward neck.
- `moseq_heading`: MoSeq's modeled orientation, converted into arena axes by
  synchronization. This is not assumed to equal anatomical head direction.

The position origin is the **same for all references**, by default the
midpoint of the ears. Thus changing direction reference does not also change
the animal position used to calculate source bearing and distance.
Landmarks and origin can be configured:

```yaml
egocentric_geometry:
  origin_landmarks: [left_ear, right_ear]
  axes:
    head:
      front: [nose]
      back: [left_ear, right_ear]
    body:
      front: [neck]
      back: [tail_base]
  angular_bins: 12
```

Each landmark list denotes its mean position. Directions run from `back`
to `front`. Coordinates and distances are in arena meters. Bearing is the
angle from the selected animal direction to the speaker, in `[-pi, pi)`;
zero means in front and positive means counterclockwise in arena x/y.
These are explicit configurable geometry definitions, not a claim of exact
reproduction of the Amaro spatial-tuning analysis.

For each sound onset, use the nearest **existing synchronized sample**, with
ties assigned to the earlier sample. The selected tracking time and signed
offset from sound onset are saved. On a complete 100-Hz grid the nearest-sample
offset is at most 5 ms. This stage does not reinterpolate tracking.

No new confidence threshold, native-frame screening, gap detection,
smoothing or pose repair is applied. Upstream tracking is used as supplied.
Undefined geometry (nonfinite coordinates or a zero-length direction vector)
and times outside the available synchronized recording are marked unavailable.
These availability checks are not a tracking-quality assessment.

## Task and speaker handling

All observed positive stimulus IDs are processed using the existing stimulus
catalog. Frequency tasks, localization tasks and passive tasks share the same
code. Background and target may use the same physical speaker.

Physical speaker positions come from the dated tables in `utils/setup.py`.
Channel 8 is excluded as a synchronization channel by default:

```yaml
egocentric_geometry:
  ignored_channels: [8]
  speaker_positions_m: null
```

Optional per-configuration `speaker_positions_m` overrides update individual
channels, for example `{3: [0.0, 0.5]}`. Coordinates must use the same arena
axes and units as synchronized tracking. Effective positions are saved in the
output. The existing speaker estimates and their historical uncertainty
remain relevant; the stage does not recalibrate them.

Exactly one nonignored physical channel is required for a source bearing.
Unknown channels, missing channel metadata and simultaneous multiple sources
receive explicit statuses. The code never silently picks the first speaker.
Their events remain in the output with unavailable bearings.

## Outputs

Under each processed session:

```text
analysis/spatial_tuning/geometry/
  geometry.h5
  summary.csv
  figures/
    coverage.pdf
    coverage.png
```

`geometry.h5` contains:

- All source sound events, plus `event_indices` linking positive-stimulus
  rows back to `meta.h5` and the Bayesian-block response stage.
- Event times, stimulus IDs, speaker coordinates, selected tracking times
  and sample-time offsets.
- Per-reference timeline origins, directions and availability flags.
- Per-reference event origins, directions, bearings, distances and flags.
- Per-stimulus angle-bin counts and direction-reference comparisons.
- Effective settings, stimulus catalog, speaker map, input hashes,
  implementation hashes and a completion marker.

`summary.csv` reports usable events, bin coverage and overlapping exclusion
counts for each reference/stimulus. Flag definitions are stored in the HDF5
attributes. Invalid angles/distances are NaN, not zero.

`coverage.pdf` shows event counts and normalized angular distributions, plus
direction-reference comparisons. These are **sampling distributions**, not
neural tuning curves. `coverage.png` previews the first page.

There is no spike input, auditory-responsiveness selection, response-window
filter, Bayesian-block dependency, bootstrap or tuning significance test.
Future tuning analyses must join events by their source indices and account
for the different background/target sampling distributions.

The rules are `egocentric_geometry` and `plot_egocentric_geometry`.
Keep Snakemake's normal rerun triggers so configuration and implementation
changes are detected. Targeting the files directly also works without enabling
the block. Selecting only head/body does not add a MoSeq dependency.

## Tests

From the repository root, using the pipeline environment:

```sh
python -m unittest discover -s tests -p 'test_egocentric_geometry.py' -v
```

Tests cover angle conventions, existing-sample selection, unchanged tracking
values, undefined axes, speaker ambiguity, stale MoSeq rejection, and complete
synthetic active/passive output generation.
