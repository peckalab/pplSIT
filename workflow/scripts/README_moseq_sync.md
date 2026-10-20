# MoSeq heading synchronization

`synch_moseq_100Hz.py` resamples the existing MoSeq CSV onto the session's
100-Hz timeline. It does not run inference or change the MoSeq model.

The `heading` column in `MoSeq.h5` is in radians, in the **same arena
coordinate system as the centroid and estimated bodypart positions**.
It uses `atan2(y, x)`: zero points along +x and positive angles turn from
+x toward +y. Values are between -pi and pi.

Heading is unwrapped before linear interpolation, then converted from
image axes to arena axes and wrapped again. This follows the shortest
angular displacement between consecutive video frames. Rotation greater
than pi between frames cannot be resolved from wrapped headings alone.
Outside the video's time interval, the nearest endpoint is retained,
matching the other synchronized columns. Invalid/nonfinite heading samples
or non-increasing video timestamps fail explicitly.

The axis signs follow the existing position conversion:
`(arena_origin - pixel) * pixel_size * (-1 if flip else 1)`.
Thus a heading-vector component keeps its sign when that axis has `flip=True`,
and changes sign otherwise. No additional `angle_compensation` is applied:
that is not part of the position conversion in this script.

Corrected files have these attributes on the `moseq` dataset:

- `heading_processing: arena_circular_v1`
- `heading_coordinate_system: arena`
- `heading_units: radians`
- `heading_interpolation: unwrap_linear_wrap`
- `heading_convention`, `heading_flip_x`, `heading_flip_y`

## Existing results

Earlier outputs interpolated across the +/-pi discontinuity as ordinary
numbers and left heading in image coordinates while transforming positions.
They lack the `heading_processing` marker. Analyses using that heading
should be regenerated after rerunning `synch_moseq_100Hz`. Syllable labels,
positions, latent states and the time grid are unchanged by this fix.
This correction does not turn MoSeq's modeled orientation into anatomical
head direction; that is a separate choice for downstream analyses.

Target the desired session's `MoSeq.h5` explicitly, using the usual pipeline
configuration and `--force` if needed to refresh an existing file. With its
existing `moseq/moseq.csv` available, `--allowed-rules synch_moseq_100Hz`
limits the run to synchronization. Inspect a dry run first. Model fitting,
inference, DLC and package installation are not required.

Regression tests:

```sh
python -m unittest discover -s tests -p 'test_moseq_sync.py' -v
```
