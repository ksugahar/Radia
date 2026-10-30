# Neutral CAD model capture

The C++ command captures the current CAD model without generating a mesh:

```text
export model_capture "model.json" units mm
```

`units mm` is a mandatory declaration by the model author. Cubit coordinates
are unitless; this command does not discover units or convert their values.
Every CAD volume must have exactly one explicitly named block, with one volume
per block. Sidesets must contain named CAD surfaces. Mesh-only sidesets and
unsupported curve types fail instead of being omitted. Existing output files
are never replaced. A failed write can leave an explicitly reported partial
file; retain it for diagnosis and choose a reviewed fresh output for a new run.

The `cubit.model-capture.v1` JSON contains source volume/face/curve/vertex IDs,
curve endpoint incidence and 33 parameter-fraction samples, exact queried sphere
parameters for a one-face whole sphere, block identity and CAD sideset coverage.
Its `labels` object uses `cubit.source-face-labels.v1`. Source IDs are not target
application IDs. Shared surface uses remain visible across source volumes.
Rational splines are explicitly rejected; sampled rational geometry is never
declared to be a polynomial spline.

This is source capture, not arbitrary exact CAD translation. The receiving
adapter must validate its supported constructions, preserve source-to-output
maps, mark sampled reconstruction as approximate, and check saved/reopened
geometry and memberships independently. Planar/analytic primitive support does
not certify a boundary-only loft's interior. Unsupported constructions must
fail before target-model mutation.

Materials and physical boundary conditions are never guessed from names. An
external adapter may take an explicit physical mapping bound to the SHA256 of
these exact capture bytes. That mapping declares units, material values and
each sideset's physical role; it is not inserted into this file with a circular
self-hash. This package contains no commercial application runtime or adapter.

AI operation uses the existing batch/nographics journal route. It must verify
the loaded plugin identity and output hash, and preserve the originating CAD.
The geometry capture does not need meshing, a field solve or a GUI session.
