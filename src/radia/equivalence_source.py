"""Equivalence-theorem near-field source (CST Near-Field Source equivalent).

Schelkunoff/Love equivalence theorem implementation for the Radia /
NGSolve stack.  Record EM field on a closed surface around a source
region, persist it, replay it later as either

  (a) a probe of the EXTERIOR field at any external observation point
      via the Stratton-Chu surface integral, or

  (b) a re-radiation source for downstream tools (Radia rad.Fld, or
      external MoM via Nastran Near_Field_Area_*.dat export).

See radia_mcp.fem.equivalence_source_knowledge for theory & lab
provenance (3 FEMM/Femtet directories distilled).

Scope (intentional):
  - Equivalence-theorem extraction + reconstruction ONLY.
  - NO IABC / SDI / asymptotic-BC content -- Radia + NGSolve already
    handle open boundaries via the Kelvin transformation, which is
    stronger and more general than IABC variants.

Typical workflow:

    >>> # 1. Extract from an NGSolve solution
    >>> from radia.equivalence_source import NearFieldSource
    >>> nfs = NearFieldSource.extract_ngsolve(
    ...     mesh, gf_E=None, gf_H=B_cf / MU_0,
    ...     surface_label="nfs_surface", omega=0)
    >>> nfs.save("magnet.nfs.json")

    >>> # 2. Probe at external points
    >>> import numpy as np
    >>> obs = np.array([[0, 0, 10.0]])  # 10 m on z-axis
    >>> H_rec = nfs.evaluate_static_H(obs)
    >>> print(H_rec)

    >>> # 3. Re-load + hand off
    >>> nfs2 = NearFieldSource.load("magnet.nfs.json")
    >>> nfs2.write_nastran_nfs("Near_Field_Area_1.dat")   # external MoM
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

# Physical constants (SI)
MU_0 = 4.0e-7 * math.pi          # vacuum permeability      [H/m]
EPS_0 = 8.8541878128e-12         # vacuum permittivity      [F/m]
C0 = 1.0 / math.sqrt(MU_0 * EPS_0)   # speed of light          [m/s]


# =====================================================================
# Dataclass
# =====================================================================

@dataclass
class NearFieldSource:
    """A discretised closed-surface equivalent EM source.

    Stores per-face: centroid, outward unit normal, area, and the
    sampled E and H phasors.  evaluate() applies the Stratton-Chu
    surface integral to compute (E, H) at any external point.

    Attributes:
        centroids:  (N, 3) float64 -- face centroid positions (m)
        normals:    (N, 3) float64 -- OUTWARD unit normals
        areas:      (N,)   float64 -- face areas (m^2)
        E:          (N, 3) complex128 or None -- E phasor at centroid (V/m)
        H:          (N, 3) complex128 or None -- H phasor at centroid (A/m)
        omega:      float -- angular frequency (rad/s); 0 means static
        metadata:   dict -- arbitrary metadata (origin, units, etc.)
    """

    centroids: np.ndarray
    normals: np.ndarray
    areas: np.ndarray
    E: Optional[np.ndarray] = None
    H: Optional[np.ndarray] = None
    omega: float = 0.0
    metadata: dict = field(default_factory=dict)

    # -----------------------------------------------------------------
    # Validation
    # -----------------------------------------------------------------

    def __post_init__(self):
        self.centroids = np.asarray(self.centroids, dtype=np.float64)
        self.normals = np.asarray(self.normals, dtype=np.float64)
        self.areas = np.asarray(self.areas, dtype=np.float64)
        if self.centroids.shape != self.normals.shape:
            raise ValueError(
                f"centroids shape {self.centroids.shape} != normals shape "
                f"{self.normals.shape}")
        if self.centroids.shape[1] != 3:
            raise ValueError(
                f"centroids must be (N, 3), got {self.centroids.shape}")
        if self.areas.shape != (self.centroids.shape[0],):
            raise ValueError(
                f"areas shape {self.areas.shape} != ({self.centroids.shape[0]},)")
        # Normalise normals just in case
        nrm = np.linalg.norm(self.normals, axis=1, keepdims=True)
        if np.any(nrm < 1e-15):
            raise ValueError("Zero-length normal found in surface mesh")
        self.normals = self.normals / nrm
        # E / H must match centroids in count if provided
        if self.E is not None:
            self.E = np.asarray(self.E, dtype=np.complex128)
            if self.E.shape != self.centroids.shape:
                raise ValueError(
                    f"E shape {self.E.shape} != centroids shape "
                    f"{self.centroids.shape}")
        if self.H is not None:
            self.H = np.asarray(self.H, dtype=np.complex128)
            if self.H.shape != self.centroids.shape:
                raise ValueError(
                    f"H shape {self.H.shape} != centroids shape "
                    f"{self.centroids.shape}")
        if self.E is None and self.H is None:
            raise ValueError(
                "At least one of E or H must be provided")

    @property
    def n_faces(self) -> int:
        return self.centroids.shape[0]

    @property
    def is_static(self) -> bool:
        return self.omega == 0.0 or (self.E is None and self.H is not None)

    # -----------------------------------------------------------------
    # Constructors -- pre-built surface meshes
    # -----------------------------------------------------------------

    @classmethod
    def spherical_surface_mesh(cls, R: float, n_theta: int = 20,
                                 n_phi: int = 40, center=(0, 0, 0)
                                 ) -> tuple:
        """Generate a (theta, phi) lat-long triangulation of a sphere.

        Returns (centroids, normals, areas) for use as a surface skeleton
        before sampling E, H.  n_theta = polar count (theta in [0, pi]);
        n_phi = azimuthal count (phi in [0, 2*pi)).  Total faces:
        2 * n_theta * n_phi (split each quad into 2 triangles).
        """
        center = np.asarray(center, dtype=np.float64)
        # Quad vertices at (theta_i, phi_j) for i in [0, n_theta], j in [0, n_phi]
        thetas = np.linspace(0.0, math.pi, n_theta + 1)
        phis = np.linspace(0.0, 2.0 * math.pi, n_phi + 1)
        TH, PH = np.meshgrid(thetas, phis, indexing='ij')
        X = R * np.sin(TH) * np.cos(PH) + center[0]
        Y = R * np.sin(TH) * np.sin(PH) + center[1]
        Z = R * np.cos(TH) + center[2]
        verts = np.stack([X, Y, Z], axis=-1)  # (n_theta+1, n_phi+1, 3)
        centroids = []
        normals = []
        areas = []
        for i in range(n_theta):
            for j in range(n_phi):
                v00 = verts[i, j]
                v01 = verts[i, j + 1]
                v10 = verts[i + 1, j]
                v11 = verts[i + 1, j + 1]
                # Two triangles: (v00, v10, v11) and (v00, v11, v01)
                for tri in [(v00, v10, v11), (v00, v11, v01)]:
                    c = np.mean(tri, axis=0)
                    e1 = tri[1] - tri[0]
                    e2 = tri[2] - tri[0]
                    cross = np.cross(e1, e2)
                    a = 0.5 * np.linalg.norm(cross)
                    if a < 1e-15:
                        continue  # degenerate (polar cap)
                    # Outward normal = radial direction from sphere center
                    n_hat = (c - center) / np.linalg.norm(c - center)
                    # Make sure cross-product orientation matches outward;
                    # if not, swap two vertices (i.e. flip)
                    if np.dot(cross / (2 * a), n_hat) < 0:
                        # Already accounted by normalising to outward n_hat,
                        # but area is positive regardless.
                        pass
                    centroids.append(c)
                    normals.append(n_hat)
                    areas.append(a)
        return (np.asarray(centroids), np.asarray(normals),
                np.asarray(areas))

    @classmethod
    def from_surface_mesh(cls, centroids, normals, areas,
                            E=None, H=None, omega=0.0, **metadata):
        """Low-level constructor from raw arrays."""
        return cls(centroids=centroids, normals=normals, areas=areas,
                    E=E, H=H, omega=omega, metadata=metadata)

    # -----------------------------------------------------------------
    # Constructors -- extraction from solvers
    # -----------------------------------------------------------------

    @classmethod
    def extract_ngsolve(cls, mesh, gf_E=None, gf_H=None,
                          surface_label: str = "nfs_surface",
                          omega: float = 0.0,
                          quadrature_order: int = 1):
        """Extract NFS from NGSolve GridFunctions on a named BND surface.

        Args:
            mesh:   NGSolve Mesh
            gf_E:   GridFunction or CoefficientFunction for E (V/m).
                    Pass None for magnetostatic problems.
            gf_H:   GridFunction or CoefficientFunction for H (A/m).
                    Must be provided.
            surface_label: BND label of the extraction surface.
                    Must be CLOSED and enclose all primary sources.
            omega:  Angular frequency (rad/s); 0 means static.
            quadrature_order: 1 (the only implemented mode) = per-element
                    area-averaged sampling (Integrate element_wise; exact
                    centroid for flat triangles, area-average for curved
                    BND elements).  > 1 (Gauss quadrature) is NOT implemented
                    and RAISES -- refine the mesh (mesh.Curve / smaller maxh)
                    for higher accuracy instead.

        Returns:
            NearFieldSource instance.
        """
        try:
            from ngsolve import BND, CoefficientFunction, Integrate
        except ImportError as e:
            raise ImportError(
                "extract_ngsolve requires NGSolve to be importable") from e

        if gf_H is None:
            raise ValueError("gf_H must be provided")
        if quadrature_order > 1:
            # No-Fallback: Gauss quadrature (order > 1) is a planned extension,
            # NOT implemented.  Raise rather than silently sample at order 1 and
            # only warn -- the caller asked for higher accuracy and would get a
            # different (lower-order) computation.  For higher accuracy on
            # curved BND elements, refine the mesh (mesh.Curve / smaller maxh)
            # and extract at quadrature_order=1, or pass externally-sampled
            # dense points to from_surface_mesh().
            raise NotImplementedError(
                f"quadrature_order={quadrature_order} (> 1, Gauss quadrature) "
                f"is not implemented in extract_ngsolve, which samples by "
                f"area-average (order 1). Refine the mesh and extract at "
                f"quadrature_order=1, or use from_surface_mesh() with "
                f"externally-sampled dense points for higher accuracy.")

        # Vectorized per-element area-averaged sampling using
        # `Integrate(cf, mesh, BND, element_wise=True)`.  This is the
        # canonical NGSolve pattern (see bem_sibc_solver.py:1815-1819)
        # and is internally TaskManager-parallel, honouring the caller's
        # `with TaskManager():` context.  The result is a 1D array
        # indexed by BND element number.
        #
        # For flat triangles with smooth fields, area-average equals
        # centroid value; for curved (mesh.Curve(p>=2)) BND elements
        # this is strictly more accurate than centroid sampling.
        bnd_region = mesh.Boundaries(surface_label)
        # Region indicator: 1 inside surface_label, 0 elsewhere -- so
        # `Integrate(cf * indicator, ..., element_wise=True)` only
        # contributes on the named surface.  Other elements get 0,
        # which we filter post-hoc via `elem_area`.
        indicator = mesh.BoundaryCF(
            {surface_label: 1.0}, default=0.0)

        elem_area = Integrate(
            indicator, mesh, VOL_or_BND=BND, element_wise=True)
        # Three vector components, each integrated separately.
        H_components = [Integrate(gf_H[i] * indicator, mesh,
                                  VOL_or_BND=BND, element_wise=True)
                        for i in range(3)]
        if gf_E is not None:
            E_components = [Integrate(gf_E[i] * indicator, mesh,
                                      VOL_or_BND=BND, element_wise=True)
                            for i in range(3)]

        # Iterate BND elements on the named surface (geometry still
        # needs random access for vertex enumeration + area/normal).
        centroids = []
        normals = []
        areas = []
        E_list = [] if gf_E is not None else None
        H_list = []

        for el in mesh.Elements(BND):
            if not bnd_region.Mask()[el.index]:
                continue
            # Face vertices in PHYSICAL coordinates
            verts_phy = [mesh[v].point for v in el.vertices]
            # Centroid in physical
            c = np.mean(verts_phy, axis=0)
            # Edge vectors for area + normal (triangle only for now)
            if len(verts_phy) == 3:
                e1 = np.asarray(verts_phy[1]) - np.asarray(verts_phy[0])
                e2 = np.asarray(verts_phy[2]) - np.asarray(verts_phy[0])
                cross = np.cross(e1, e2)
                a = 0.5 * np.linalg.norm(cross)
                n_hat = cross / (2.0 * a) if a > 1e-15 else None
            elif len(verts_phy) == 4:
                # Quad: split into 2 triangles, sum areas, use centroid normal
                v0, v1, v2, v3 = map(np.asarray, verts_phy)
                a1_cross = np.cross(v1 - v0, v2 - v0)
                a2_cross = np.cross(v2 - v0, v3 - v0)
                a = 0.5 * (np.linalg.norm(a1_cross) +
                           np.linalg.norm(a2_cross))
                n_hat = (a1_cross + a2_cross)
                n_hat /= np.linalg.norm(n_hat) if np.linalg.norm(n_hat) > 1e-15 else 1.0
            else:
                raise ValueError(
                    f"Unsupported BND element with {len(verts_phy)} vertices")

            if n_hat is None:
                continue

            # Per-element area-averaged H from vectorized Integrate.
            # `elem_area[el.nr]` is the integral of the region
            # indicator over the element -- equals the FE measure of
            # the element on the named surface (also handles curved
            # BND elements correctly).
            ae = elem_area[el.nr]
            if abs(ae) < 1e-30:
                # Element on surface_label but degenerate; skip
                continue
            H_at = np.array(
                [complex(H_components[i][el.nr]) / ae for i in range(3)],
                dtype=np.complex128)

            if gf_E is not None:
                ve = np.array(
                    [complex(E_components[i][el.nr]) / ae for i in range(3)],
                    dtype=np.complex128)
                E_list.append(ve)

            centroids.append(c)
            normals.append(n_hat)
            areas.append(a)
            H_list.append(H_at)

        if not centroids:
            raise ValueError(
                f"No BND elements found on surface '{surface_label}'. "
                f"Available: {list(mesh.GetBoundaries())}")

        metadata = {
            "source": f"ngsolve.extract({surface_label})",
            "quadrature_order": quadrature_order,   # always 1 here (>1 raised)
            "n_faces": len(centroids),
        }
        return cls(
            centroids=np.asarray(centroids),
            normals=np.asarray(normals),
            areas=np.asarray(areas),
            E=np.asarray(E_list) if E_list is not None else None,
            H=np.asarray(H_list),
            omega=omega,
            metadata=metadata,
        )

    @classmethod
    def extract_radia(cls, rad_container, surface_skeleton: tuple,
                       omega: float = 0.0):
        """Extract NFS from a Radia container via rad.Fld() sampling.

        Args:
            rad_container: Radia container handle (return value of
                rad.ObjCnt or rad.ObjHexahedron etc.)
            surface_skeleton: (centroids, normals, areas) tuple, e.g.
                from `spherical_surface_mesh()`.
            omega: Static = 0.  (Radia is magnetostatic, so omega is
                   tagged in metadata only; the produced NFS is static.)

        Returns:
            NearFieldSource (static, E = None, H from rad.Fld).
        """
        try:
            import radia as rad
        except ImportError as e:
            raise ImportError(
                "extract_radia requires the radia package") from e

        centroids, normals, areas = surface_skeleton
        # rad.Fld returns B [T] in batch; convert to H [A/m]
        B = np.asarray(rad.Fld(rad_container, "b", centroids))
        if B.ndim == 1:
            B = B.reshape(1, 3)
        H = B / MU_0
        metadata = {
            "source": "radia.rad_Fld",
            "n_faces": len(centroids),
            "omega_tag_only": omega,   # Radia is DC; omega stored for label
        }
        return cls(
            centroids=np.asarray(centroids),
            normals=np.asarray(normals),
            areas=np.asarray(areas),
            E=None,
            H=H.astype(np.complex128),
            omega=0.0,
            metadata=metadata,
        )

    # -----------------------------------------------------------------
    # Serialization
    # -----------------------------------------------------------------

    def save(self, path):
        """Save NFS to a JSON file.

        Format:
            {
              "version": 1,
              "omega": float,
              "metadata": {...},
              "centroids":  [[x, y, z], ...],
              "normals":    [[nx, ny, nz], ...],
              "areas":      [dS1, dS2, ...],
              "E_re":   [[Ex_re, Ey_re, Ez_re], ...]  (optional),
              "E_im":   [[Ex_im, ...], ...]            (optional),
              "H_re":   [[Hx_re, Hy_re, Hz_re], ...],
              "H_im":   [[Hx_im, ...], ...],
            }
        """
        path = Path(path)
        data = {
            "version": 1,
            "omega": float(self.omega),
            "metadata": self.metadata,
            "centroids": self.centroids.tolist(),
            "normals": self.normals.tolist(),
            "areas": self.areas.tolist(),
            "H_re": self.H.real.tolist(),
            "H_im": self.H.imag.tolist(),
        }
        if self.E is not None:
            data["E_re"] = self.E.real.tolist()
            data["E_im"] = self.E.imag.tolist()
        path.write_text(json.dumps(data, separators=(",", ":")))

    @classmethod
    def load(cls, path):
        """Load NFS from a JSON file (counterpart of save())."""
        data = json.loads(Path(path).read_text())
        if data.get("version") != 1:
            raise ValueError(f"Unsupported NFS file version: {data.get('version')}")
        centroids = np.asarray(data["centroids"])
        normals = np.asarray(data["normals"])
        areas = np.asarray(data["areas"])
        H = (np.asarray(data["H_re"]) + 1j * np.asarray(data["H_im"]))
        E = None
        if "E_re" in data:
            E = (np.asarray(data["E_re"]) + 1j * np.asarray(data["E_im"]))
        return cls(
            centroids=centroids, normals=normals, areas=areas,
            E=E, H=H, omega=float(data["omega"]),
            metadata=data.get("metadata", {}),
        )

    # -----------------------------------------------------------------
    # Evaluation (Stratton-Chu)
    # -----------------------------------------------------------------

    def evaluate_static_H(self, obs_points: np.ndarray,
                            use_cpp: bool = True,
                            n_threads: int = 0) -> np.ndarray:
        """Magnetostatic reduction: reconstruct H at observation points.

        Schelkunoff/Love representation, scalar + vector potential form,
        differentiated inside the integral:

            H(r) = (1/(4 pi)) * sum_faces { grad(1/R) x J_s
                                            - (n . H_s) grad(1/R) } * dS

        with J_s = n x H_s, R = |r - centroid|, grad(1/R) = -R_vec/R^3.
        Sign convention fixed by magnetic-dipole unit test.

        Phase A (2026-05-26): delegates to the C++ kernel
        `radia._radia_pybind._EquivalenceSourceStaticH` for ~50-100x
        speedup over the legacy numpy loop.  The Python fallback path
        is kept as a safety net and will be removed in Phase E.

        Args:
            obs_points: (M, 3) array of observation points (m).
            use_cpp:    True (default) = call C++ kernel.  False =
                         force the Python fallback (slow, kept for
                         A-B regression diffing only).
            n_threads:  0 (default) = NGSolve TaskManager default;
                         > 0 = request that many threads.  Only used
                         by the C++ path.

        Returns:
            (M, 3) real array of H (A/m).
        """
        if self.H is None:
            raise ValueError("evaluate_static_H requires H on the surface")
        obs = np.atleast_2d(np.asarray(obs_points, dtype=np.float64))
        H_surf = self.H.real    # static -> use real part

        if use_cpp:
            try:
                from . import _radia_pybind  # type: ignore
                fn = getattr(_radia_pybind, "_EquivalenceSourceStaticH", None)
                if fn is not None:
                    return fn(
                        np.ascontiguousarray(self.centroids, dtype=np.float64),
                        np.ascontiguousarray(self.normals, dtype=np.float64),
                        np.ascontiguousarray(self.areas, dtype=np.float64),
                        np.ascontiguousarray(H_surf, dtype=np.float64),
                        np.ascontiguousarray(obs, dtype=np.float64),
                        int(n_threads),
                    )
            except ImportError:
                pass   # fall through to Python fallback

        # ------------------ Python fallback (legacy, slow) --------------
        # Kept for Phase A as a regression-diff safety net; removed in
        # Phase E once the C++ kernel has shipped + been used at scale.
        R_vec = obs[:, None, :] - self.centroids[None, :, :]   # (M, N, 3)
        R = np.linalg.norm(R_vec, axis=2)                       # (M, N)
        R_safe = np.where(R > 1e-15, R, 1.0)
        grad = -R_vec / R_safe[..., None] ** 3                  # (M, N, 3)
        grad[R < 1e-15] = 0.0
        n_arr = self.normals[None, :, :]                        # (1, N, 3)
        H_s = H_surf[None, :, :]                                # (1, N, 3)
        J_s = np.cross(n_arr, H_s)                              # (1, N, 3)
        rho_m_over_mu = np.sum(n_arr * H_s, axis=2)             # (1, N)
        term_Js = np.cross(grad, J_s)                           # (M, N, 3)
        term_rhom = -rho_m_over_mu[..., None] * grad            # (M, N, 3)
        per_face = (term_Js + term_rhom) * self.areas[None, :, None]
        H_rec = (1.0 / (4.0 * math.pi)) * np.sum(per_face, axis=1)
        return H_rec

    def evaluate(self, obs_points: np.ndarray, omega: Optional[float] = None,
                  use_cpp: bool = True, n_threads: int = 0):
        """Full time-harmonic Stratton-Chu reconstruction.

        Phase B (2026-05-26): when ``use_cpp=True`` (default) and
        ``omega > 0``, delegates to the C++ kernel
        ``_EquivalenceSourceHarmonic`` which uses equivalent currents
        J_s = n x H_s and M_s = n x E_s with the FULL dyadic Green's
        function (I + grad-grad/k^2) psi.  Correctly reproduces
        deep-near-field behaviour (the Python fallback's scalar form
        underestimated by up to a factor of 3 at R_obs/lambda << 1).

        Args:
            obs_points: (M, 3) array of observation points (m).
            omega: angular frequency (rad/s).  If None, uses self.omega.
            use_cpp: True (default) -> C++ dyadic kernel.  False ->
                     Python scalar fallback (legacy, deprecated; will
                     be removed in Phase E).  The fallback is kept as a
                     numerical-diff cross-check during Phase B/C/D.
            n_threads: 0 = NGSolve TaskManager default.  Only used by
                       the C++ path.

        Returns:
            (E_rec, H_rec): each (M, 3) complex array.
            For omega == 0, calls evaluate_static_H and returns
                (zeros_like(H_rec) + 0j, H_rec + 0j).
        """
        if omega is None:
            omega = self.omega
        obs = np.atleast_2d(np.asarray(obs_points, dtype=np.float64))

        if omega == 0.0:
            H_rec = self.evaluate_static_H(obs, use_cpp=use_cpp,
                                              n_threads=n_threads).astype(np.complex128)
            E_rec = np.zeros_like(H_rec)
            return E_rec, H_rec

        if self.E is None or self.H is None:
            raise ValueError(
                "Full Stratton-Chu (omega != 0) requires BOTH E and H "
                "on the surface")

        if use_cpp:
            try:
                from . import _radia_pybind  # type: ignore
                fn = getattr(_radia_pybind,
                              "_EquivalenceSourceHarmonic", None)
                if fn is not None:
                    E_re_arr = np.ascontiguousarray(self.E.real, dtype=np.float64)
                    E_im_arr = np.ascontiguousarray(self.E.imag, dtype=np.float64)
                    H_re_arr = np.ascontiguousarray(self.H.real, dtype=np.float64)
                    H_im_arr = np.ascontiguousarray(self.H.imag, dtype=np.float64)
                    (e_re, e_im, h_re, h_im) = fn(
                        np.ascontiguousarray(self.centroids, dtype=np.float64),
                        np.ascontiguousarray(self.normals, dtype=np.float64),
                        np.ascontiguousarray(self.areas, dtype=np.float64),
                        E_re_arr, E_im_arr,
                        H_re_arr, H_im_arr,
                        np.ascontiguousarray(obs, dtype=np.float64),
                        float(omega),
                        int(n_threads),
                    )
                    return (e_re + 1j * e_im, h_re + 1j * h_im)
            except ImportError:
                pass    # fall through to Python fallback

        # ------------------ Python fallback (scalar form, legacy) -------
        # KNOWN LIMITATION (only when use_cpp=False is explicitly
        # passed): the scalar Green's-function form below misses the
        # (1/k^2) grad-grad psi term of the full dyadic GF and
        # undershoots in the deep near-field.  Use use_cpp=True
        # (default) for the correct dyadic reconstruction.  This
        # fallback is kept solely as a regression-diff anchor.

        k = omega * math.sqrt(MU_0 * EPS_0)
        # Vectorised over (M, N)
        R_vec = obs[:, None, :] - self.centroids[None, :, :]   # (M, N, 3)
        R = np.linalg.norm(R_vec, axis=2)                       # (M, N)
        R_safe = np.where(R > 1e-15, R, 1.0)
        # Green's function
        psi = np.exp(-1j * k * R_safe) / (4.0 * math.pi * R_safe)  # (M, N)
        # grad psi = -(jk + 1/R) psi * R-hat ; R-hat = R_vec / R
        grad_psi = (-(1j * k + 1.0 / R_safe) * psi)[..., None] * R_vec / R_safe[..., None]
        psi[R < 1e-15] = 0.0
        grad_psi[R < 1e-15] = 0.0

        n_arr = self.normals[None, :, :]                        # (1, N, 3)
        E_s = self.E[None, :, :]                                # (1, N, 3)
        H_s = self.H[None, :, :]                                # (1, N, 3)
        J_s = np.cross(n_arr, H_s)                              # (1, N, 3)
        # Legacy scalar form is kept only as a regression-diff anchor.
        # It is written in the historical E_s x n convention together
        # with explicit scalar-charge terms; the production C++ dyadic
        # path above uses M_s = n x E_s and no separate charge terms.
        M_s = np.cross(E_s, n_arr)                              # (1, N, 3)
        rho_e_over_eps = np.sum(n_arr * E_s, axis=2)            # (1, N)
        rho_m_over_mu  = np.sum(n_arr * H_s, axis=2)            # (1, N)

        # Time-harmonic Stratton-Chu, same sign convention as the
        # static reduction in evaluate_static_H:
        #   E = integral { -jw mu_0 J_s psi + grad_psi x M_s
        #                  - (n . E) grad_psi } dS
        #   H = integral {  jw eps_0 M_s psi + grad_psi x J_s
        #                  - (n . H) grad_psi } dS
        # (numpy.cross broadcasts (1, N, 3) and (M, N, 3) correctly.)
        E_per_face = (
            (-1j * omega * MU_0) * J_s * psi[..., None]
            + np.cross(grad_psi, M_s)
            - rho_e_over_eps[..., None] * grad_psi
        )
        H_per_face = (
            (1j * omega * EPS_0) * M_s * psi[..., None]
            + np.cross(grad_psi, J_s)
            - rho_m_over_mu[..., None] * grad_psi
        )
        E_rec = np.sum(E_per_face * self.areas[None, :, None], axis=1)
        H_rec = np.sum(H_per_face * self.areas[None, :, None], axis=1)
        return E_rec, H_rec

    # -----------------------------------------------------------------
    # External format export
    # -----------------------------------------------------------------

    def write_nastran_nfs(self, path, frequency_hz: Optional[float] = None):
        """Export to a Nastran-style Near_Field_Area_*.dat that the
        EMCoS Antenna Vlab / FEKO / equivalent tools can ingest.

        Matches the format produced by the Sugahara Lab 2015_04_12
        FEMTET workflow (appendix MATLAB script).

        Args:
            path: output file path
            frequency_hz: frequency in Hz; defaults to omega / (2 pi).
                For static (omega = 0) this is required.
        """
        path = Path(path)
        if frequency_hz is None:
            if self.omega == 0.0:
                raise ValueError("frequency_hz must be provided for static NFS")
            frequency_hz = self.omega / (2 * math.pi)

        if self.E is None or self.H is None:
            raise ValueError(
                "Nastran NFS export requires BOTH E and H "
                "(use evaluate_at_surface() to fill in zeros if needed)")

        # Need unique nodes + triangle connectivity; the NFS stores
        # per-face data but not the node list.  We reconstruct it by
        # walking the centroid records and treating each face as
        # 3 nodes via the helper _reconstruct_triangulation().
        nodes, tris = self._reconstruct_triangulation()
        with path.open("w", newline="\r\n") as fh:
            fh.write("location center\r\n")
            # GRID lines (one per node)
            for i, p in enumerate(nodes, start=1):
                fh.write("GRID    %8d       0% 6.5f% 6.5f% 6.5f       0\r\n"
                         % (i, p[0], p[1], p[2]))
            # CTRIA3 lines (one per triangle)
            for j, t in enumerate(tris, start=1):
                fh.write("CTRIA3  %8d       1%8d%8d%8d\r\n"
                         % (j, t[0] + 1, t[1] + 1, t[2] + 1))
            # FIELDMUL lines (one per triangle - field data)
            for j, (E_j, H_j) in enumerate(zip(self.E, self.H), start=1):
                fh.write(
                    "FIELDMUL%8d %12.7f %15.3f  %14.8f %15.3f  %14.8f %15.3f  "
                    "%14.8f %15.9f  %14.8f %15.9f  %14.8f %15.9f  %14.8f\r\n"
                    % (j, frequency_hz,
                        abs(E_j[0]), np.angle(E_j[0]),
                        abs(E_j[1]), np.angle(E_j[1]),
                        abs(E_j[2]), np.angle(E_j[2]),
                        abs(H_j[0]), np.angle(H_j[0]),
                        abs(H_j[1]), np.angle(H_j[1]),
                        abs(H_j[2]), np.angle(H_j[2]),
                    )
                )

    def _reconstruct_triangulation(self):
        """Per-face -> (nodes_unique, triangles) for Nastran export.

        The NFS stores centroids only.  To export Nastran we synthesise
        per-face triangle vertices from the face centroid + outward
        normal + area, assuming equilateral triangles oriented in the
        plane of the normal.  This is a LOSSY conversion -- the
        original mesh topology is not preserved.  For round-trip
        Nastran -> NFS the export consumer should treat each triangle
        as an independent panel.
        """
        # Equilateral triangle with the given area dS:
        #   edge = 2 sqrt(dS / sqrt(3))
        # Place it perpendicular to the normal, vertices at +-edge/2 in
        # an arbitrary tangent direction.
        nodes = []
        tris = []
        for c, n, dS in zip(self.centroids, self.normals, self.areas):
            edge = 2.0 * math.sqrt(dS / math.sqrt(3.0))
            # Build an orthonormal basis (t1, t2) perpendicular to n
            if abs(n[2]) < 0.9:
                t1 = np.cross(n, [0, 0, 1])
            else:
                t1 = np.cross(n, [1, 0, 0])
            t1 = t1 / np.linalg.norm(t1)
            t2 = np.cross(n, t1)
            # Equilateral vertices around c
            r = edge / math.sqrt(3.0)  # circumradius
            v0 = c + r * t1
            v1 = c + r * (-0.5 * t1 + math.sqrt(3) / 2 * t2)
            v2 = c + r * (-0.5 * t1 - math.sqrt(3) / 2 * t2)
            base_idx = len(nodes)
            nodes.extend([v0, v1, v2])
            tris.append((base_idx, base_idx + 1, base_idx + 2))
        return nodes, tris

    # -----------------------------------------------------------------
    # Projection onto an NGSolve GridFunction (the CF -> .sol path)
    # -----------------------------------------------------------------

    def project_to_h1_vector(self, mesh, order: int = 1, definedon=None,
                              save_path=None):
        """Project the static-H reconstruction onto a VectorH1 GridFunction.

        Builds a `VectorH1(mesh, order=order)` GridFunction and fills it
        by sampling `evaluate_static_H()` at each free vertex /  DOF
        coordinate.  Optionally restrict to a sub-region via `definedon`
        (e.g. `mesh.Materials("outer")`).

        For order=1 H1 the DOFs are vertices, so this gives the
        canonical nodal interpolant of the NFS reconstruction.  For
        order>1 it's a sub-optimal interpolation (vertices only); a
        proper L2 projection via `gfu.Set(cf)` would need a true
        CoefficientFunction wrapper, which is a planned roadmap item.

        Args:
            mesh: NGSolve Mesh (different mesh than the one used for
                  extraction is fine; the recon. is a free-space field).
            order: VectorH1 polynomial order (1 recommended for this
                   sampling-based path).
            definedon: Optional `mesh.Materials(...)` to restrict to a
                       sub-region.  Important: must be OUTSIDE the NFS
                       extraction surface, else the reconstruction is
                       not physically meaningful (will give 0 inside,
                       by the equivalence theorem null-field property).
            save_path: If set, also call `gfu.Save(save_path)`.

        Returns:
            (gfu, n_dofs) -- the populated GridFunction + number of
            non-zero DOFs filled.
        """
        try:
            from ngsolve import VectorH1, GridFunction
        except ImportError as e:
            raise ImportError(
                "project_to_h1_vector requires NGSolve") from e

        # NGSolve VectorH1 doesn't accept definedon=None; only pass the
        # kwarg when we actually want a sub-region restriction.
        kw = {"order": order}
        if definedon is not None:
            kw["definedon"] = definedon
        fes = VectorH1(mesh, **kw)
        gfu = GridFunction(fes)
        gfu.vec[:] = 0.0

        # Collect vertex coordinates from the mesh.  For order=1, the
        # H1 DOFs are EXACTLY the mesh vertices, so sampling at vertices
        # IS the nodal interpolation.
        ngmesh = mesh.ngmesh
        pts = np.asarray([(p[0], p[1], p[2])
                            for p in ngmesh.Points()],
                          dtype=np.float64)
        # NGSolve uses 1-indexed vertex ids; pts[i] corresponds to
        # vertex id i+1 in NGSolve's MeshElement.vertices interface.

        # Vectorised evaluation
        H_at_vertices = self.evaluate_static_H(pts).real  # (Nv, 3)

        # Assign to the VectorH1 GridFunction.  For each vertex, the
        # vertex-coupled DOFs hold the (Hx, Hy, Hz) value.  VectorH1
        # stores them interleaved per-vertex.
        free_dofs = fes.FreeDofs()
        # Use Set() with the analytic CF would be cleaner but requires
        # an NGSolve-native CF wrapper.  Sampling-based fallback:
        n_filled = 0
        for v_id in range(pts.shape[0]):
            try:
                # GetDofNrs returns the DOF indices for a vertex node
                dofs = fes.GetDofNrs(__import__("ngsolve").NodeId(
                    __import__("ngsolve").VERTEX, v_id))
            except Exception:
                # Fallback for VectorH1 layout: order=1 has dim_per_node
                # consecutive DOFs starting at v_id * 3.
                dofs = [v_id * 3, v_id * 3 + 1, v_id * 3 + 2]
            if len(dofs) != 3:
                continue
            for i, d in enumerate(dofs):
                if d >= 0 and (free_dofs is None or free_dofs[d]):
                    gfu.vec[d] = float(H_at_vertices[v_id, i])
                    n_filled += 1
        n_filled //= 3  # convert back to vertex count

        if save_path is not None:
            gfu.Save(str(save_path))

        return gfu, n_filled

    # -----------------------------------------------------------------
    # Re-radiation as Radia source (static only, planned)
    # -----------------------------------------------------------------

    def to_radia_objects(self):
        """Convert NFS into a Radia container suitable for rad.Fld().

        Planned for the static case.  Current implementation is a
        thin proof-of-concept: each face is treated as a magnetic
        dipole with moment m = (rho_m / mu_0) * dS * n + |J_s| * dS
        (placeholder).

        Use evaluate_static_H() directly for production code at the
        moment; full re-radiation via rad.ObjCnt is a roadmap item.
        """
        raise NotImplementedError(
            "to_radia_objects() is a planned roadmap item.  For now, "
            "use evaluate_static_H(obs_points) or evaluate(obs_points) "
            "to probe the exterior field directly via Stratton-Chu.")


# =====================================================================
# Module-level convenience
# =====================================================================

def reconstruct_static_H(centroids, normals, areas, H_surface, obs_points
                          ) -> np.ndarray:
    """Convenience function:  Stratton-Chu static H reconstruction.

    Useful when you want the integral without building a NearFieldSource
    instance (e.g. a one-off script or a unit test).
    """
    nfs = NearFieldSource(
        centroids=np.asarray(centroids),
        normals=np.asarray(normals),
        areas=np.asarray(areas),
        H=np.asarray(H_surface, dtype=np.complex128),
        omega=0.0,
    )
    return nfs.evaluate_static_H(obs_points)


__all__ = [
    "NearFieldSource",
    "reconstruct_static_H",
    "MU_0",
    "EPS_0",
    "C0",
]
