"""
Regression tests: field evaluation on TrfOrnt-transformed containers/elements
and RadiaField CoefficientFunction evaluation during NGSolve assembly.

Locks the 2026-07-10/11 fixes for two native-crash mechanisms (0xC0000374 /
0xC0000005 heap corruption, hard process death without a Python traceback):

1. radTGroup::B_genComp propagated the group's transforms by MUTATING every
   child's g3dListOfTransform (push_front + restore).  The batch field
   kernel (rad.Fld with an (N,3) array) evaluates points from parallel
   TaskManager threads, so a TrfOrnt-wrapped container corrupted the heap.
   Groups now inherit the non-mutating radTg3d::B_genComp (NestedFor_B).
   The same fix removed radTPolyhedron::B_genComp, which silently IGNORED
   the element's own transform list (TrfOrnt on a hexahedron evaluated the
   untransformed field).

2. RadiaFieldCF::Evaluate round-tripped through Python rad.Fld, whose batch
   kernel starts an internal ngcore ParallelFor.  NGSolve assembly calls
   Evaluate from INSIDE a running TaskManager job; ngcore job state is
   static, so the nested CreateJob (and the GIL save/restore on worker
   threads) corrupted the running job.  The CF now evaluates GIL-free via
   the serial C entries (RadFldBatchSerial etc.) -- parallelism stays in
   NGSolve's element loop.

Crash-prone paths run in a SUBPROCESS so a regression can never take down
the pytest process.
"""

import os
import subprocess
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../src"))

import radia as rad

REPO_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "../src"))

ROT_CEN = [0.0, -0.1, 0.0]
ROT_AXIS = [0.0, 0.0, 1.0]
ROT_ANGLE = -0.349  # rad, about +z


def _build_current_container():
	"""Container of TrfOrnt-positioned current elements (CoilBuilder-style)."""
	objs = []
	bar1 = rad.ObjRecCur([0, 0, 0], [0.01, 0.2, 0.01], [0, 3.0e6, 0])
	objs.append(rad.TrfOrnt(bar1, rad.TrfTrsl([0.08, 0.0, 0.03])))
	bar2 = rad.ObjRecCur([0, 0, 0], [0.01, 0.2, 0.01], [0, -3.0e6, 0])
	objs.append(rad.TrfOrnt(bar2, rad.TrfTrsl([-0.08, 0.0, 0.03])))
	arc1 = rad.ObjArcCur([0, 0, 0], [0.075, 0.085], [0.0, np.pi], 0.01,
	                     10, "auto", "z", 3.0e6)
	objs.append(rad.TrfOrnt(arc1, rad.TrfTrsl([0.0, 0.1, 0.03])))
	arc2 = rad.ObjArcCur([0, 0, 0], [0.075, 0.085], [np.pi, 2.0 * np.pi],
	                     0.01, 10, "auto", "z", 3.0e6)
	objs.append(rad.TrfOrnt(arc2, rad.TrfTrsl([0.0, -0.1, 0.03])))
	return rad.ObjCnt(objs)


def _rotation_matrix(axis, angle):
	ax = np.asarray(axis, float)
	ax = ax / np.linalg.norm(ax)
	c, s = np.cos(angle), np.sin(angle)
	K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
	return np.eye(3) * c + s * K + (1 - c) * np.outer(ax, ax)


def _assert_vec_close(actual, expected, rtol, msg):
	"""Vector-difference comparison (norm(a-b) <= rtol*norm(b)), per the
	field-comparison policy -- componentwise rtol fails on symmetric-zero
	components that only carry floating-point noise."""
	actual = np.asarray(actual, float)
	expected = np.asarray(expected, float)
	ref = np.linalg.norm(expected)
	assert ref > 0.0, "%s: reference vector is zero" % msg
	err = np.linalg.norm(actual - expected) / ref
	assert err <= rtol, "%s: rel err %.3e > %.1e (%r vs %r)" % (
		msg, err, rtol, actual, expected)


def _run_subprocess(script):
	env = dict(os.environ)
	env["MKL_NUM_THREADS"] = "1"
	env["OMP_NUM_THREADS"] = "1"
	return subprocess.run([sys.executable, "-c", script],
	                      capture_output=True, text=True, timeout=600, env=env)


class TestBatchFldTransformedContainer:
	"""Batch rad.Fld on a TrfOrnt-wrapped container (the crash path).

	Runs in a subprocess: before the 2026-07-10 fix this died with
	0xC0000374 / 0xC0000005 inside the batch ParallelFor.
	"""

	def test_batch_fld_no_crash_and_matches_reference(self):
		script = f"""
import os, sys
sys.path.insert(0, {REPO_SRC!r})
import numpy as np
import radia as rad

ROT_CEN = np.array({ROT_CEN!r})
ROT_ANGLE = {ROT_ANGLE!r}

def build():
	objs = []
	bar1 = rad.ObjRecCur([0, 0, 0], [0.01, 0.2, 0.01], [0, 3.0e6, 0])
	objs.append(rad.TrfOrnt(bar1, rad.TrfTrsl([0.08, 0.0, 0.03])))
	bar2 = rad.ObjRecCur([0, 0, 0], [0.01, 0.2, 0.01], [0, -3.0e6, 0])
	objs.append(rad.TrfOrnt(bar2, rad.TrfTrsl([-0.08, 0.0, 0.03])))
	arc1 = rad.ObjArcCur([0, 0, 0], [0.075, 0.085], [0.0, np.pi], 0.01,
	                     10, "auto", "z", 3.0e6)
	objs.append(rad.TrfOrnt(arc1, rad.TrfTrsl([0.0, 0.1, 0.03])))
	arc2 = rad.ObjArcCur([0, 0, 0], [0.075, 0.085], [np.pi, 2.0 * np.pi],
	                     0.01, 10, "auto", "z", 3.0e6)
	objs.append(rad.TrfOrnt(arc2, rad.TrfTrsl([0.0, -0.1, 0.03])))
	return rad.ObjCnt(objs)

c, s = np.cos(ROT_ANGLE), np.sin(ROT_ANGLE)
R = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])

rng = np.random.default_rng(0)
pts = rng.uniform(-0.2, 0.2, size=(20000, 3))

# Reference: PLAIN container evaluated at inverse-rotated points, fields rotated.
rad.UtiDelAll()
cnt_plain = rad.ObjCnt([build()])
pts_local = (pts - ROT_CEN) @ R + ROT_CEN  # row-vector form of R^-1 (p - c) + c
H_ref = np.asarray(rad.Fld(cnt_plain, "h", pts_local)) @ R.T

# Transformed container: batch evaluation (crashed before the fix).
rad.UtiDelAll()
cnt = rad.ObjCnt([build()])
rad.TrfOrnt(cnt, rad.TrfRot(ROT_CEN.tolist(), [0.0, 0.0, 1.0], ROT_ANGLE))
H_batch = np.asarray(rad.Fld(cnt, "h", pts))

# Batch (parallel) must agree with single-point (serial classic path).
for i in (0, 1234, 7777, 19999):
	H_single = np.asarray(rad.Fld(cnt, "h", pts[i].tolist()))
	assert np.allclose(H_batch[i], H_single, rtol=1e-9, atol=1e-9), (
		"batch vs single mismatch at %d: %r vs %r" % (i, H_batch[i], H_single))

# Transform semantics: rotated container == rotated reference solution.
scale = np.abs(H_ref).max()
assert scale > 0.0
err = np.abs(H_batch - H_ref).max() / scale
assert err < 1e-9, "rotated-container field mismatch: rel err %.3e" % err
print("BATCH-TRF-OK maxrel=%.3e scale=%.6e" % (err, scale))
"""
		res = _run_subprocess(script)
		assert res.returncode == 0, (
			"batch rad.Fld on transformed container crashed or failed: "
			"exit=%s (0x%08X)\nstdout:\n%s\nstderr:\n%s"
			% (res.returncode, res.returncode & 0xFFFFFFFF,
			   res.stdout[-2000:], res.stderr[-2000:]))
		assert "BATCH-TRF-OK" in res.stdout


class TestRadiaFieldCFTransformedContainer:
	"""RadiaField CoefficientFunction on a TrfOrnt-wrapped container.

	The original crash site: NGSolve LinearForm assembly batch-evaluates the
	CF, which delegates to rad.Fld batch.  Runs in a subprocess.
	"""

	def test_linearform_assembly_and_pointwise_match(self):
		pytest.importorskip("ngsolve")
		script = f"""
import os, sys
sys.path.insert(0, {REPO_SRC!r})
import numpy as np
import radia as rad
from ngsolve import H1, LinearForm, Mesh, TaskManager, dx, grad
from netgen.occ import Box, OCCGeometry, Pnt

ROT_CEN = {ROT_CEN!r}
ROT_ANGLE = {ROT_ANGLE!r}

def build():
	objs = []
	bar1 = rad.ObjRecCur([0, 0, 0], [0.01, 0.2, 0.01], [0, 3.0e6, 0])
	objs.append(rad.TrfOrnt(bar1, rad.TrfTrsl([0.08, 0.0, 0.03])))
	bar2 = rad.ObjRecCur([0, 0, 0], [0.01, 0.2, 0.01], [0, -3.0e6, 0])
	objs.append(rad.TrfOrnt(bar2, rad.TrfTrsl([-0.08, 0.0, 0.03])))
	arc1 = rad.ObjArcCur([0, 0, 0], [0.075, 0.085], [0.0, np.pi], 0.01,
	                     10, "auto", "z", 3.0e6)
	objs.append(rad.TrfOrnt(arc1, rad.TrfTrsl([0.0, 0.1, 0.03])))
	arc2 = rad.ObjArcCur([0, 0, 0], [0.075, 0.085], [np.pi, 2.0 * np.pi],
	                     0.01, 10, "auto", "z", 3.0e6)
	objs.append(rad.TrfOrnt(arc2, rad.TrfTrsl([0.0, -0.1, 0.03])))
	return rad.ObjCnt(objs)

rad.UtiDelAll()
cnt = rad.ObjCnt([build()])
rad.TrfOrnt(cnt, rad.TrfRot(ROT_CEN, [0.0, 0.0, 1.0], ROT_ANGLE))

box = Box(Pnt(-0.2, -0.2, -0.2), Pnt(0.2, 0.2, 0.2))
box.mat("air")
with TaskManager():
	mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.08))
	Hs = rad.RadiaField(cnt, "h")
	fes = H1(mesh, order=2)
	u, v = fes.TnT()
	f = LinearForm(fes)
	f += Hs * grad(v) * dx
	f.Assemble()  # crashed with heap corruption before the fix
nrm = f.vec.Norm()
assert np.isfinite(nrm) and nrm > 0.0, "assembled RHS norm invalid: %r" % nrm

# CF pointwise values must match rad.Fld on the same transformed container.
probes = [(0.05, 0.02, 0.01), (-0.1, 0.05, 0.04), (0.0, -0.15, 0.02)]
for p in probes:
	cf_val = np.asarray(Hs(mesh(*p)))
	fld_val = np.asarray(rad.Fld(cnt, "h", list(p)))
	assert np.allclose(cf_val, fld_val, rtol=1e-9, atol=1e-9), (
		"RadiaField CF vs rad.Fld mismatch at %r: %r vs %r"
		% (p, cf_val, fld_val))
	assert np.abs(fld_val).max() > 0.0, "zero field at probe %r" % (p,)
print("CF-TRF-OK |f|=%.6e" % nrm)
"""
		res = _run_subprocess(script)
		assert res.returncode == 0, (
			"RadiaField CF assembly on transformed container crashed or "
			"failed: exit=%s (0x%08X)\nstdout:\n%s\nstderr:\n%s"
			% (res.returncode, res.returncode & 0xFFFFFFFF,
			   res.stdout[-2000:], res.stderr[-2000:]))
		assert "CF-TRF-OK" in res.stdout


class TestTrfOrntPolyhedronApplied:
	"""TrfOrnt on a polyhedron element must move its field (it was silently
	ignored by the removed radTPolyhedron::B_genComp override).

	Single-point evaluations are serial and safe in-process.
	"""

	def test_hexahedron_translation_applied(self):
		rad.UtiDelAll()
		mag = rad.magnet_box([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 954930.0])
		probe = np.array([0.02, 0.003, 0.015])
		B_orig = np.asarray(rad.Fld(mag, "b", probe.tolist()))
		shift = np.array([0.05, 0.0, 0.0])
		rad.TrfOrnt(mag, rad.TrfTrsl(shift.tolist()))
		B_moved = np.asarray(rad.Fld(mag, "b", (probe + shift).tolist()))
		_assert_vec_close(B_moved, B_orig, 1e-10,
		                  "translated hexahedron field does not follow the element")
		rad.UtiDelAll()

	def test_hexahedron_rotation_applied(self):
		rad.UtiDelAll()
		mag = rad.magnet_box([0.02, 0, 0], [0.01, 0.01, 0.01], [954930.0, 0, 0])
		probe = np.array([0.05, 0.004, 0.006])
		B_orig = np.asarray(rad.Fld(mag, "b", probe.tolist()))
		angle = 0.5 * np.pi
		R = _rotation_matrix([0, 0, 1], angle)
		rad.TrfOrnt(mag, rad.TrfRot([0, 0, 0], [0, 0, 1], angle))
		B_rot = np.asarray(rad.Fld(mag, "b", (R @ probe).tolist()))
		_assert_vec_close(B_rot, R @ B_orig, 1e-10,
		                  "rotated hexahedron field mismatch")
		rad.UtiDelAll()

	def test_container_of_hexahedra_translation_applied(self):
		rad.UtiDelAll()
		mag1 = rad.magnet_box([0, 0, 0], [0.008, 0.008, 0.008], [0, 0, 954930.0])
		mag2 = rad.magnet_box([0.02, 0, 0], [0.008, 0.008, 0.008], [0, 0, 954930.0])
		cnt = rad.ObjCnt([mag1, mag2])
		probe = np.array([0.01, 0.005, 0.02])
		B_orig = np.asarray(rad.Fld(cnt, "b", probe.tolist()))
		shift = np.array([0.0, 0.06, 0.0])
		rad.TrfOrnt(cnt, rad.TrfTrsl(shift.tolist()))
		B_moved = np.asarray(rad.Fld(cnt, "b", (probe + shift).tolist()))
		_assert_vec_close(B_moved, B_orig, 1e-10,
		                  "translated container field does not follow the container")
		rad.UtiDelAll()


if __name__ == "__main__":
	pytest.main([__file__, "-v"])
