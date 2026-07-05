"""Metamaterial knowledge: physics, design and lab cross-links.

Distilled from PDFs in
  public-safe curated corpus
  public-safe curated corpus
  public-safe curated corpus
  public-safe curated corpus

Source papers
-------------
- Sato (2005), "Research Trends in Metamaterials", R&D Review of Toyota
  CRDL, Vol. 40, No. 2 -- LH materials, SRR, CRLH transmission lines,
  EBG, topology-optimization of dielectric microstructure.
- Sadatgol, Rahman, Forati, Levy, Guney (2016), "Enhanced Faraday
  rotation in hybrid magneto-optical metamaterial structure of
  bismuth-substituted-iron-garnet with embedded-gold-wires",
  J. Appl. Phys. 119, 103105.
- Bergamin (2009), "Negative index of refraction, spacetime folding
  and perfect imaging in transformation optics", arXiv:0912.2649.
- Vesalago (1968), Pendry (1996, 1999, 2000, 2006), Smith (2000),
  Shelby (2001), Leonhardt (2006) -- foundational references.
- "ジャイレータを用いた左手系伝送線路に関する研究" (lab thesis,
  password-protected -- summary only).
- Kiziltas (2003), Nomura (2005), Sigmund and Jensen (2004) --
  topology optimization for EM materials.

This file is the single source of truth for the dispatcher's accepted
topic enum.  Wired into the `metamaterial(topic)` MCP tool and the
companion `metamaterial_topics()` tool via common.register_topics_tool.
"""
from __future__ import annotations


TOPICS: dict[str, str] = {
    "overview":
        "Field map: DNG / LH / NIR materials, EBG, CRLH transmission "
        "line, SRR, transformation optics. Lab connection to Kelvin "
        "transform (radia_mcp.electromagnet).",
    "lh_materials_neg_index":
        "Left-handed (Veselago 1968) materials with simultaneously "
        "negative epsilon and mu. Inverted Snell, reversed Doppler, "
        "reversed Cherenkov radiation.",
    "srr_split_ring":
        "Pendry split-ring resonator (1999): equivalent LC circuit, "
        "negative-mu band, mu_eff(omega) Lorentzian, wire grid for "
        "negative-epsilon, Smith composite (2000).",
    "effective_medium":
        "Homogenization of periodic sub-wavelength structures. "
        "Wire-medium plasma frequency, Drude-like epsilon_eff, "
        "MOM permittivity tensor, gyration g, retrieval procedure.",
    "transformation_optics_kelvin":
        "Pendry-Leonhardt transformation optics: coordinate map -> "
        "anisotropic epsilon, mu. Cloaking. Kelvin inversion as a "
        "special conformal map -- direct link to radia_mcp.electromagnet "
        "open-boundary Kelvin transform for accelerator magnets.",
    "magneto_optical":
        "Faraday rotation, bismuth-substituted iron garnet (Bi:YIG), "
        "MOM with embedded gold wires (Sadatgol 2016): 9x enhancement, "
        "permittivity tensor with off-diagonal gyration g.",
    "topology_optimization_design":
        "Density-method topology optimization for periodic EM "
        "microstructure (Toyota CRDL, Kiziltas, Sigmund-Jensen). "
        "Inverse design of band-gap, NIR, antenna substrate.",
    "crlh_transmission_line":
        "Composite Right/Left-Handed (CRLH) transmission line "
        "(Caloz-Itoh, Iyer-Eleftheriades, Sanada, Sato). Series C_L "
        "+ shunt L_L gives LH band; parasitic L_R, C_R add RH band. "
        "Balanced CRLH -> seamless LH<->RH transition without bandgap.",
    "gyrator_lh_line":
        "Gyrator-based left-handed transmission line (lab thesis, "
        "password-protected -- summary). Non-reciprocal LH line "
        "constructed from active gyrator unit cells.",
    "perfect_lens_evanescent":
        "Pendry perfect lens (2000): n=-1 slab amplifies evanescent "
        "modes -> sub-wavelength imaging. Bergamin (2009) caution: "
        "transformation-optics lens is NOT the same -- it triples "
        "sources, does NOT amplify evanescent waves, restricted to "
        "stationary regime.",
    "applications_automotive":
        "Toyota CRDL roadmap (Sato 2005): mm-wave LH leaky-wave "
        "antenna for 76-77 GHz ACC radar, EBG substrate antennas, "
        "absorbers, beam-steering, EMC shielding, motor magnetic "
        "materials.",
    "acoustic_sonic_crystal":
        "Acoustic metamaterials / sonic crystals: Bloch band gaps need "
        "confinement or transverse periodicity (a sparse free-space "
        "chain has NONE - verified 4 ways); unit-cell + duct-transmission "
        "methodology with empty-lattice and port gates; FEM/BEM "
        "validation discipline (invisibility null gate, Anderson "
        "fluid-sphere series, irregular frequencies, CHIEF vs "
        "Burton-Miller).",
    "all":
        "Concatenated dump of every topic above.",
}


# ============================================================
# Individual topic strings
# ============================================================


OVERVIEW = r"""
# Metamaterials: physics, design and lab cross-links

A *metamaterial* is an engineered periodic composite whose
sub-wavelength unit cells produce an effective electromagnetic
response not available in natural media.  The most striking examples
exhibit simultaneously negative permittivity and permeability
(left-handed / negative-index materials), but the family also
includes electromagnetic bandgap (EBG) structures, complex
high-impedance surfaces, magneto-optical (MO) composites, and
transformation-optics media (cloaks, lenses).

## Field map

| Class | Defining feature | Lab connection |
|-------|------------------|----------------|
| DNG / LH (Vesalago 1968) | epsilon < 0 AND mu < 0 | radia.peec for the equivalent-circuit view |
| NIR (negative index of refraction) | n = -sqrt(eps mu) < 0 | analysis with radia.bem / FDTD |
| EBG (electromagnetic band-gap) | photonic bandgap, surface-wave suppression | radia_mcp.radia_ngsolve periodic FES |
| CRLH transmission line | series C_L + shunt L_L (LH) + parasitics (RH) | radia_mcp.peec circuit extraction |
| SRR (split-ring resonator) | resonant negative-mu | radia_mcp.electromagnet (analytical formulas) |
| Magneto-optical metamaterial (MOM) | tailored Faraday rotation via wire mesh | radia.bem for full-wave check |
| Transformation-optics media | coordinate map -> anisotropic eps, mu | **radia_mcp.electromagnet Kelvin** (see transformation_optics_kelvin) |

## Historical timeline

| Year | Milestone | Reference |
|------|-----------|-----------|
| 1968 | Vesalago predicts DNG materials | Sov. Phys. Usp. 10 (4), 509-514 |
| 1996 | Pendry: low-frequency plasmons in wire array | Phys. Rev. Lett. 76, 4773 |
| 1999 | Pendry: SRR magnetic resonance | IEEE Trans. MTT 47, 2075 |
| 2000 | Smith: first composite DNG experiment | Phys. Rev. Lett. 84, 4184 |
| 2000 | Pendry: perfect lens (evanescent amplification) | Phys. Rev. Lett. 85, 3966 |
| 2001 | Shelby: experimental negative refraction | Science 292, 77 |
| 2002 | Caloz-Itoh, Iyer-Eleftheriades: CRLH TL | USNC/URSI, IEEE MTT-S |
| 2003 | Kiziltas: topology optimization for EM substrates | IEEE Trans. AP 51, 2732 |
| 2006 | Pendry-Schurig-Smith: cloaking via transformation optics | Science 312, 1780 |
| 2006 | Leonhardt: optical conformal mapping | Science 312, 1777 |
| 2009 | Bergamin: TO lens triples sources, does NOT amplify evanescent | arXiv:0912.2649 |
| 2016 | Sadatgol: 9x Faraday rotation in Bi:YIG + Au-wire MOM | J. Appl. Phys. 119, 103105 |

## Why metamaterials matter for Sugahara Lab

Three of the lab's core competencies map directly onto metamaterial
design:

1. **Kelvin transformation** for unbounded-domain accelerator
   magnets (`radia_mcp.electromagnet`) is the **same family of
   coordinate transformations** that Pendry/Leonhardt use for
   cloaks and transformation-optics lenses.  See the
   `transformation_optics_kelvin` topic.

2. **Topology optimization** for iron-yoke shaping and PM
   placement is the same algorithm Toyota CRDL applies to EBG
   substrates and Sigmund/Jensen apply to photonic-crystal
   waveguides (`topology_optimization_design`).

3. **PEEC + SIBC** for induction heating gives equivalent-circuit
   parameters that map directly onto CRLH TL unit-cell parameters
   (C_L, L_L, C_R, L_R) -- the metamaterial leaky-wave antenna
   is just a PEEC ladder network with frequency dispersion
   (`crlh_transmission_line`).

## Cross-references

- `radia_mcp.electromagnet` -- Kelvin transform, CoilBuilder
- `radia_mcp.radia_ngsolve` -- Periodic FES, BEM, analytical formulas
- `radia_mcp.peec` -- PEEC ladder networks (CRLH cell)
- `radia_mcp.ih` -- ESIM cell problem (effective surface impedance)
- `radia_mcp.fem` -- Periodic and Floquet conditions
"""


LH_MATERIALS_NEG_INDEX = r"""
# Left-handed materials and negative index of refraction

## Vesalago's prediction (1968)

V. G. Vesalago, "The Electrodynamics of Substances with
Simultaneously Negative Values of epsilon and mu", Sov. Phys. Usp.
10 (4), 509-514, 1968.  The original DNG (double-negative) paper.

Starting from Maxwell's equations in a homogeneous medium with
simultaneously negative permittivity and permeability:

    eps < 0,  mu < 0   =>   n = -sqrt(eps mu)

A medium with negative n is *left-handed*: the triad (E, H, k) forms
a left-handed orthogonal system, instead of the usual right-handed
one.  This produces several counter-intuitive effects:

| Effect | Sign in DNG medium |
|--------|---------------------|
| Snell refraction | **inverted** (ray bends to the same side of the normal as the incident ray) |
| Doppler shift | **reversed** (source moving toward observer -> red shift) |
| Cherenkov radiation | **reversed** cone (backward) |
| Energy flow S = E x H | **opposite** to phase velocity k |
| Group velocity v_g | **anti-parallel** to phase velocity v_p |

Vesalago noted that no natural material exhibits eps < 0 AND mu < 0
simultaneously over a useful frequency band; the prediction sat
dormant for 30 years until Pendry and Smith engineered artificial
composites.

## Pendry's recipes

### Negative epsilon via wire array (Pendry 1996)

J. B. Pendry et al., "Extremely Low Frequency Plasmons in Metallic
Microstructures", Phys. Rev. Lett. 76 (25), 4773-4776, 1996.

A periodic array of thin conducting wires with spacing a much
smaller than wavelength behaves as an effective plasma with:

    eps_eff(omega) = 1 - omega_p^2 / (omega (omega + i gamma))

The plasma frequency is dramatically reduced relative to a bulk metal
because only a tiny fraction of the unit cell volume is conductor:

    omega_p^2 = (2 pi c^2) / (a^2 ln(a / r))

where r is wire radius and a is the lattice spacing.  For a = 5 mm
and r = 50 micron, omega_p falls into the GHz range -- entirely
classical, no quantum effects.

### Negative mu via split-ring resonator (Pendry 1999)

J. B. Pendry et al., "Magnetism from Conductors and Enhanced
Nonlinear Phenomena", IEEE Trans. MTT 47, 2075-2084, 1999.

A non-magnetic metallic split-ring resonator exhibits an LC
resonance.  Time-varying H along the ring axis induces a current
that produces a magnetic dipole moment -- with a Lorentzian
frequency response that yields mu_eff < 0 in a band just above
resonance.  See the `srr_split_ring` topic.

### Composite DNG (Smith 2000)

D. R. Smith et al., "Composite Medium with simultaneously Negative
Permeability and Permittivity", Phys. Rev. Lett. 84 (18), 4184-4187,
2000.

Combine the wire array (negative eps) and SRR (negative mu) at
overlapping frequencies -> first experimental DNG medium.
Confirmed by Shelby et al. (Science 292, 77, 2001) using a
microwave prism geometry that directly measures n < 0.

## Loss and bandwidth limitations

Resonant SRR + wire media:
- Narrow bandwidth (set by SRR Q-factor)
- Significant ohmic loss near resonance
- Difficult to scale beyond microwave

These limitations motivated:
- CRLH transmission lines (Caloz, Iyer-Eleftheriades, 2002):
  non-resonant, broader band, lower loss -- see
  `crlh_transmission_line`.
- Transformation-optics cloaks (Pendry 2006, Leonhardt 2006):
  use geometric coordinate maps instead of resonant unit cells --
  see `transformation_optics_kelvin`.

## Lab connection

The PEEC ladder network in `radia_mcp.peec` produces the same
equivalent-circuit description that Vesalago postulated for DNG
media.  Sub-wavelength PEEC unit cells with negative effective
inductance / capacitance behave as Vesalago slabs in the LF limit.
"""


SRR_SPLIT_RING = r"""
# Split-ring resonator (SRR)

## Geometry and equivalent circuit

The classic SRR (Pendry 1999, IEEE Trans. MTT 47, 2075) is a pair
of concentric metal rings, each with a small gap, etched on a
dielectric substrate.  The double-ring configuration with
diametrically opposite gaps produces strong capacitive coupling
between the rings and shifts the resonance into the GHz / mm-wave
band for sub-millimetre dimensions.

Equivalent LC circuit:
- L from current loop of one ring
- C from gap and ring-to-ring capacitance
- Resonance: omega_0 = 1 / sqrt(L C)

The unit cell's effective permeability (Pendry 1999):

    mu_eff(omega) = 1 - F omega^2 / (omega^2 - omega_0^2 + i Gamma omega)

with F the geometric filling factor and Gamma the loss term.
Negative mu band: omega_0 < omega < omega_0 / sqrt(1 - F).

## Polarization sensitivity

SRR responds only to H normal to the ring plane.  For 3D isotropic
DNG, three orthogonal SRR orientations per unit cell are needed.

## Variants

| Variant | Lab purpose |
|---------|-------------|
| Double-ring SRR | original Pendry; resonance at 5-10 GHz |
| Edge-coupled SRR (EC-SRR) | broadside coupling, lower omega_0 |
| Broadside-coupled SRR (BC-SRR) | strong capacitive coupling -> smaller cell |
| Spiral resonator | multi-turn, even lower omega_0 |
| Complementary SRR (CSRR) | Babinet dual on ground plane; negative epsilon |

## Lab connection: surface impedance

The SRR resonance produces an effective frequency-dependent surface
impedance similar to what `radia_mcp.ih` ESIM computes for
nonlinear iron workpieces.  The skin-depth + susceptibility
formalism is the same:

    Z_s(omega) = (1 + j) / (sigma delta)   (linear case)

For the magneto-quasi-static (MQS) regime, the SRR equivalent
circuit can be wired into a `radia.peec.PEECBuilder` as a
single-port resonator and used as a metamaterial-aware load in
WPT or IH ladder networks.

## Smith composite (negative epsilon + negative mu)

D. R. Smith et al., "Composite Medium with simultaneously Negative
Permeability and Permittivity", Phys. Rev. Lett. 84, 4184, 2000.

The first experimental DNG medium combined Pendry wire arrays
(eps < 0 below diluted plasma frequency) with SRRs (mu < 0 in a
narrow band above SRR resonance) tuned to overlap in the same band.
Shelby et al. (Science 292, 77, 2001) verified negative refraction
in a prism geometry.

## References

- J. B. Pendry et al., IEEE Trans. MTT 47 (1999) 2075-2084.
- D. R. Smith et al., Phys. Rev. Lett. 84 (2000) 4184-4187.
- D. R. Smith, N. Kroll, Phys. Rev. Lett. 85 (2000) 2933-2936.
- R. A. Shelby et al., Science 292 (2001) 77-79.
- K. Sarabandi, H. Mosallaei, IEEE AP-S Int. Symp. 3 (2003) 355.

## Cross-references

- `crlh_transmission_line` -- non-resonant alternative to SRR
- `radia_mcp.ih` -- ESIM cell problem (related surface impedance)
- `radia_mcp.peec` -- PEEC ladder + equivalent circuit
"""


EFFECTIVE_MEDIUM = r"""
# Effective medium / homogenization

## When the EM approximation is valid

The effective medium (homogenization) treatment of a periodic
sub-wavelength composite is valid when the lattice constant a is
much smaller than the wavelength inside the medium:

    a << lambda / n_eff

In practice a / lambda < 0.1 is "safe", a / lambda > 0.25 starts to
show spatial-dispersion (non-local) effects and the homogenized
parameters depend on direction of propagation.

## Wire-medium plasma frequency

For a square array of thin parallel wires (radius r, spacing a) in
a dielectric host of permittivity eps_h, embedded in vacuum:

    (k_p a)^2 = 2 pi / ln( a^2 / (4 r (a - r)) )

This gives the diluted plasma wavenumber.  Above k_p, the wire
medium is dielectric (Re(eps_eff) > 0); below, it is metallic
(Re(eps_eff) < 0) and TE waves polarized along the wires are
evanescent.

The Drude form for the wire medium's effective permittivity
(neglecting losses):

    eps_eff(omega) = eps_h - eps_h k_p^2 / k_0^2

where k_0 = omega / c.

## MOM permittivity tensor (Sadatgol 2016)

For a magneto-optical host with permittivity tensor

    eps_MO = [[eps_perp, i g, 0],
              [-i g, eps_perp, 0],
              [0, 0, eps_par]]

the addition of an embedded wire mesh (in X and Y directions only,
preserving the Z-axis magnetization) gives the homogenized MOM
permittivity tensor:

    eps_MOM ~= eps_MO - kappa * I_xy

where kappa is the wire-mesh contribution

    kappa = ( k_0^2 / k_p^2  -  1 / (f_v (eps_m - eps_perp)) )^(-1)

f_v = pi r^2 / a^2 is the volumetric wire-filling fraction, and
eps_m is the metal permittivity (typically Drude gold or silver).
I_xy = diag(1, 1, 0).

Key consequence: the diagonal Re(eps_perp) - kappa can be tuned
small while gyration g is preserved.  Faraday rotation per length

    d theta / dz = (omega sqrt(mu)) / (2 c_0) * (sqrt(eps_perp + g) - sqrt(eps_perp - g))

scales as ~g / sqrt(eps_perp) and is maximized when eps_perp - kappa
approaches g -- giving the 9x rotation enhancement reported in
Bi:YIG + Au-wire MOM at telecom wavelength.

## Retrieval procedure (Smith-Vier-Koschny-Soukoulis)

For a finite metamaterial slab the effective n and Z can be
retrieved from simulated or measured S parameters:

1. From transmission T and reflection R, compute
       cos(n k_0 d) = (1 - R^2 + T^2) / (2 T)
   and
       Z = sqrt( ((1+R)^2 - T^2) / ((1-R)^2 - T^2) )
2. Choose branch of arccos using passivity (Im(n) >= 0).
3. eps = n/Z,  mu = n Z.

For chiral / gyrotropic media (MOM, MO composites) the procedure
extends to LCP and RCP eigenmodes:

    n_+ = sqrt(mu (eps_perp + g)),   n_- = sqrt(mu (eps_perp - g))

retrieved separately from LCP/RCP T, R -- see Sadatgol (2016)
Section III.

## Spatial dispersion warning

When the wire spacing approaches a /lambda ~ 0.25, the wire medium
shows non-local (k-dependent) eps -- the simple Drude formula
breaks down and a hydrodynamic correction is needed.  Silveirinha
(Phys. Rev. B 79, 035118, 2009) and Hanson et al. (IEEE Trans. AP
60 (9), 4219, 2012) give the rigorous spatially-dispersive
homogenization used in MOM derivation.

## References

- M. G. Silveirinha, Phys. Rev. B 79 (3) (2009) 035118.
- G. W. Hanson, E. Forati, M. G. Silveirinha, IEEE Trans. AP 60 (9)
  (2012) 4219-4231.
- D. R. Smith, D. C. Vier, Th. Koschny, C. M. Soukoulis,
  Phys. Rev. E 71 (3) (2005) 036617.
- Th. Koschny et al., Phys. Rev. B 71 (24) (2005) 245105.

## Cross-references

- `srr_split_ring` -- resonant unit cell case
- `magneto_optical` -- application of MOM homogenization
- `radia_mcp.radia_ngsolve` -- periodic FES (microscopic check)
"""


TRANSFORMATION_OPTICS_KELVIN = r"""
# Transformation optics and the Kelvin transform connection

## The big idea

Transformation optics (Pendry-Schurig-Smith 2006; Leonhardt 2006)
recasts a coordinate diffeomorphism in vacuum as a position-
dependent anisotropic constitutive relation for an equivalent
medium in flat space.  This lets us *design* a medium that bends,
focuses, or hides light along any prescribed coordinate map.

Starting from Maxwell's equations in vacuum and a coordinate map
x^i -> x_bar^i (x), the medium that reproduces the transformed
solution in the laboratory frame has

    eps^ij = +/- ( det(L)^{-1} ) L^i_k L^j_l delta^kl  (and similarly mu)

where L^i_k = partial x_bar^i / partial x^k.  The +/- sign comes
from the orientation of the map (Bergamin 2009; see below).

This is the *same* mathematical operation that the lab's Kelvin
transformation uses to map an open-domain accelerator-magnet
problem to a bounded computational domain -- the two are different
applications of a single coordinate-transformation framework.

## Pendry cloaking (2006)

J. B. Pendry, D. Schurig, D. R. Smith, "Controlling Electromagnetic
Fields", Science 312 (2006) 1780-1782.

Map a spherical hole (radius R_1 < r < R_2) in the laboratory to a
single point r_bar = 0 in electromagnetic space.  Any ray emanating
from the origin in electromagnetic space wraps smoothly around the
hole in laboratory space -- the hole and its contents are invisible
to external observers.

The required medium parameters (spherical cloak):
    eps_r = mu_r = ((r - R_1)/r)^2 * R_2 / (R_2 - R_1)
    eps_theta = mu_theta = R_2 / (R_2 - R_1)
    eps_phi = mu_phi = R_2 / (R_2 - R_1)

These are extreme anisotropic profiles realized in microwave bands
with SRR + wire metamaterials.

## Leonhardt conformal mapping (2006)

U. Leonhardt, "Optical Conformal Mapping", Science 312 (2006) 1777.

Uses 2D conformal maps (w = f(z) in complex coordinates) to design
isotropic-index media that route light around objects.  In 2D, a
conformal map preserves angles and gives an *isotropic* n(x, y)
profile, much easier to fabricate than the anisotropic 3D
Pendry-Schurig-Smith cloak.

## Kelvin inversion as a special transformation-optics map

The Kelvin (sphere) inversion map

    x -> x_bar = R^2 x / |x|^2

is a *special conformal transformation* of 3-space: it maps
spheres to spheres and preserves angles.  Applied to a Laplace /
magnetostatic problem, it has the remarkable property that **a
harmonic function in the original domain transforms into a
harmonic function in the inverted domain** -- with a known
multiplicative factor (1/|x|).

This is exactly the property the lab exploits in
`radia_mcp.electromagnet`:

> Map the unbounded exterior r > R of an accelerator magnet to the
> compact ball r_bar < R via Kelvin inversion.  Solve the
> Laplace / Poisson equation for the magnetic scalar potential
> Omega on the compact mesh.  The Kelvin-image boundary condition
> at the sphere r = R guarantees the correct far-field behaviour
> (Omega -> 0 at infinity in the original domain).

Mathematically, this is the same coordinate-transformation
technology Pendry / Leonhardt use -- it just happens that the
target is a *vacuum Laplace problem* (no constitutive jumps)
rather than a *medium-design* problem.

Practical consequence: the discrete Kelvin mesh identification
infrastructure in the lab (sphere-paired vertices, OCC Identify,
post-hoc Kelvin via `kelvin-identify-post-hoc` skill) is reusable
for transformation-optics cloak meshing.

## Bergamin's caveat on negative index

L. Bergamin, "Negative Index of Refraction, Spacetime Folding and
Perfect Imaging in Transformation Optics", arXiv:0912.2649, 2009.

Negative-index emerges in transformation optics only when the
coordinate map *changes the orientation of space* (e.g. z_bar = -z
or a swap x_bar = y, y_bar = x).  The +/- sign is then a deliberate
*choice* in the re-interpretation of the transformed Maxwell
solution -- not an inherent geometric feature.  Three competing
prescriptions (A, B, C in Table 1 of Bergamin) distribute the
signs differently:

| Prescription | Negative index when |
|--------------|---------------------|
| (A) Leonhardt-Philbin | spatial orientation flips |
| (B) Bergamin minimal-sign | spacetime orientation flips |
| (C) "no sign change anywhere" | NEVER (no negative index) |

Implication for cloak design: the spherical Pendry cloak does NOT
require a negative-index region (orientation-preserving map);
**flat** perfect lenses (Pendry-Veselago) do require it -- and even
then, Bergamin shows the transformation-designed lens does not
amplify evanescent modes (it triples sources instead).

This is a subtle point worth flagging when reading the literature:
"transformation-optics flat lens" and "Pendry-Veselago flat lens"
are NOT the same physical device, even though both have n = -1 in
the slab.

## Cross-references

- `radia_mcp.electromagnet` -- Kelvin transform for open boundaries
- `radia_mcp.radia_ngsolve.kelvin_periodic_identifications`
- `lh_materials_neg_index` -- the n < 0 story
- `perfect_lens_evanescent` -- the Bergamin caveat in detail

## References

- J. B. Pendry, D. Schurig, D. R. Smith, Science 312 (2006) 1780.
- U. Leonhardt, Science 312 (2006) 1777.
- U. Leonhardt, T. G. Philbin, New J. Phys. 8 (2006) 247.
- L. Bergamin, arXiv:0912.2649 (2009).
- L. Bergamin, Phys. Rev. A 78 (2008) 043825.
- L. Bergamin, Phys. Rev. A (2009), arXiv:09031821.
- L. D. Landau, E. M. Lifshitz, "The Classical Theory of Fields",
  4th ed., Butterworth-Heinemann (2006) -- covariant formulation.
"""


MAGNETO_OPTICAL = r"""
# Magneto-optical metamaterials (MOM)

## The Faraday effect in natural media

When a linearly polarized wave propagates through a magnetic medium
along the direction of magnetization, Zeeman / spin-orbit splitting
of electronic states makes the left- and right-circular eigenmodes
travel at different phase velocities.  The plane of polarization
rotates by an angle proportional to length:

    delta theta = V * B * L     (V is Verdet constant)

For bismuth-substituted yttrium iron garnet (Bi:YIG) at telecom
1550 nm, the specific Faraday rotation at saturation magnetization
is 0.02 -- 1 deg / micron, requiring tens of microns of material
for the ~45 deg rotation needed in an optical isolator -- too thick
for integrated photonics.

## Enhancing Faraday rotation: the MOM approach (Sadatgol 2016)

M. Sadatgol, M. Rahman, E. Forati, M. Levy, D. O. Guney,
"Enhanced Faraday rotation in hybrid magneto-optical metamaterial
structure of bismuth-substituted-iron-garnet with embedded-gold-
wires", J. Appl. Phys. 119 (2016) 103105.

**Key idea**: embed a non-MO wire mesh in the MO host.  The wire
mesh tunes only the *diagonal* elements of the MO permittivity
tensor (not the off-diagonal gyration g).  By bringing the diagonal
eps_perp close to the gyration g (just above the diluted plasma
frequency), the Faraday rotation per unit length is enhanced
**9x** without loss of MO performance.

### Permittivity tensor

Bi:YIG host, magnetized in Z:

    eps_MO = [[eps_perp, i g, 0],
              [-i g, eps_perp, 0],
              [0, 0, eps_par]]

At 1550 nm: eps_perp ~= 5, g ~= 0.01, mu ~= 1, near-zero absorption.

### MOM with embedded wires

Gold wires (radius r = 17 nm) in X and Y directions, unit cell
a = 130 nm.  After homogenization (Sec. effective_medium for the
formula):

    eps_MOM ~= eps_MO - kappa * diag(1, 1, 0)

Faraday rotation per unit length

    d theta / dz = omega sqrt(mu) / (2 c_0) *
                   ( sqrt(eps_perp - kappa + g) - sqrt(eps_perp - kappa - g) )

Maximized when eps_perp - kappa ~= g (just above plasma frequency
where both LCP and RCP modes propagate but eps - kappa is small).

### Results

For 16 unit cells (~ 2.1 micron) at 190 THz (1.58 micron):
- 9x Faraday rotation enhancement vs pure Bi:YIG
- Polarization extinction ratio < -20 dB over > 2 THz bandwidth
- Attenuation ~ 1.5 dB / micron (gold ohmic loss)
- Magneto-optical Kerr effect (MOKE) enhancement ~ 10x at 190.2 THz

### Practical alternative: planar multilayer

The wire-mesh feature size (17 nm radius) is challenging.  Sadatgol
(2016) shows that a periodic stack of 6-nm gold films + 124-nm
Bi:YIG layers gives essentially the same enhancement, fabricable by
standard RF sputtering.

## Other MO metamaterials in the literature

| Approach | Enhancement | Bandwidth |
|----------|-------------|-----------|
| Fabry-Perot cavity + MO (Inoue, Fujii 1997) | high | narrow |
| Plasmonic nanostructures + MO (Chin 2013) | moderate | narrow |
| MOM wire-mesh (Sadatgol 2016) | 9x | > 2 THz |
| MOM multilayer (Sadatgol 2016) | 9x | > 2 THz |
| Gold grating + graphene + Si (APL Photonics 3, 016103, 2018) | high | narrow |

## Lab connection

- `radia_mcp.peec` for the equivalent electromagnetic model of the
  wire mesh (LF / quasi-static limit).
- Faraday-rotation MO sensors and isolators are standard components
  in cryogenic accelerator-magnet diagnostics; the MOM design route
  could shrink isolator size by 10x.
- The off-diagonal gyration g formalism is analogous to the chiral
  media handled in `radia_mcp.radia_ngsolve.bem`.

## References

- M. Sadatgol et al., J. Appl. Phys. 119 (2016) 103105.
- J. Y. Chin et al., Nat. Commun. 4 (2013) 1599.
- V. I. Belotelov et al., Nat. Nanotechnol. 6 (2011) 370.
- V. V. Temnov et al., Nat. Photonics 4 (2010) 107.
- M. Inoue, T. Fujii, J. Appl. Phys. 81 (1997) 5659.
- V. J. Fratello, R. Wolfe, "Epitaxial garnet films for non-
  reciprocal magneto optic devices", in *Handbook of Thin Film
  Devices*, vol. 4, Academic Press (2000) 93-141.
- L. Bi et al., Nat. Photonics 5 (2011) 758 (telecom isolator).
"""


TOPOLOGY_OPTIMIZATION_DESIGN = r"""
# Topology optimization of electromagnetic materials

Topology optimization (TO) -- the algorithmic search for an
optimal material distribution within a design domain -- is one of
the most general design tools available for metamaterials, because
it does not presuppose any specific unit-cell topology.  TO is the
recommended route when the desired electromagnetic response cannot
be matched to a known canonical structure (SRR, wire mesh, mushroom
EBG, ...).

## Density-method TO for periodic EM microstructure

Following Bendsoe and Sigmund (structural mechanics) and adapted to
EM by Toyota CRDL (Sato 2005; Nomura 2005):

1. Discretize the unit cell into a fine FE mesh.
2. Assign a *design density* rho_e in [0, 1] to each element.
   rho_e = 0 -> void (air, eps_r0 = 1).
   rho_e = 1 -> bulk material (e.g. eps_r = 10).
   Intermediate values represent porous / mixed regions.
3. Set up a forward FE solve for the figure of merit (band-gap
   width, |S21|^2 in a pass-band, antenna bandwidth, ...).
4. Compute design-variable sensitivities by the **adjoint method**
   (Hang and Arora 1978) -- a single additional FE solve gives
   d(FOM)/d(rho_e) for *all* elements simultaneously.
5. Update rho via a nonlinear-programming step (MMA, OC, SLP).
6. Penalize intermediate densities via SIMP:
       eps(rho) = 1 + (eps_max - 1) * rho^p,   p = 3
   to drive the final design toward 0/1.

The microstructure emerges from random / homogeneous initial
conditions; the optimizer discovers shapes that often look nothing
like SRR or wire-mesh geometries.

## Examples from the literature

| Application | Authors | Reference |
|-------------|---------|-----------|
| Microstrip antenna dielectric substrate | Kiziltas et al. | IEEE TAP 51 (2003) 2732 |
| EBG dielectric periodic substrate | Nomura et al. | 6th WCSMO (2005) |
| Photonic-crystal waveguide bends | Sigmund, Jensen | APL 84 (2004) 2022 |
| Magneto-static topology | Yoo, Hong | IJSS 41 (2004) 2461 |
| EBG for surface-wave suppression | Sievenpiper et al. | IEEE TMT 47 (1999) 2059 |

## Multi-physics extension (Toyota CRDL)

The Toyota CRDL flowchart (Sato 2005, Fig. 3) couples three FE
solves -- electromagnetic + structural + thermal -- within a single
density-method TO loop.  This is the natural way to design
multi-functional metamaterials where the EM band-gap must coexist
with mechanical stiffness or heat dissipation constraints.

## Lab connection: TO for accelerator magnets and IH workpieces

The same density-method TO machinery applies to:

1. **Iron yoke shape optimization** for accelerator magnets
   (`radia_mcp.electromagnet`) -- FOM = field uniformity or
   multipole purity in the good-field region.  The lab's Hantila
   polarization method (`radia_mcp.electromagnet.hantila`) gives
   the LU-factored-once forward solve that TO needs (~100x speedup
   over Picard).
2. **PM placement** for IH coil array design
   (`radia_mcp.ih`) -- FOM = uniform workpiece heating.
3. **Cubit hex mesh + density variable** -- the hex-mesh sweep
   topology in `radia_mcp.cubit` provides exactly the structured
   unit-cell grid that TO requires.

Cross-references:

- `radia_mcp.electromagnet.hantila` -- fast inverse for nonlinear
  iron (TO inner loop)
- `radia_mcp.radia_ngsolve.adjoint` -- adjoint formulation for FEM
- `crlh_transmission_line` -- the canonical "non-TO" alternative

## References

- G. Kiziltas et al., IEEE Trans. AP 51 (2003) 2732-2743.
- T. Nomura et al., Proc. 6th WCSMO (2005).
- O. Sigmund, J. Jensen, Appl. Phys. Lett. 84 (12) (2004) 2022-2024.
- J. Yoo, H. Hong, Int. J. Solids Struct. 41 (2004) 2461-2477.
- E. Hang, S. Arora, Comput. Methods Appl. Mech. Eng. 15 (1978)
  35-62.
- K. Sato, "Research Trends in Metamaterials", R&D Review of Toyota
  CRDL 40 (2) (2005) 24-30.
- M. P. Bendsoe, O. Sigmund, "Topology Optimization: Theory,
  Methods, and Applications", Springer (2004).
"""


CRLH_TRANSMISSION_LINE = r"""
# Composite Right/Left-Handed (CRLH) transmission line

## Why an alternative to SRR

The Pendry-Smith DNG composite (SRR + wires) is *resonant* -- it
exhibits negative parameters only in a narrow band near the SRR
resonance and is intrinsically lossy.  In 2002 two groups
independently proposed a *non-resonant*, *circuit-oriented*
alternative: the composite right/left-handed transmission line.

- C. Caloz, A. Sanada, T. Itoh -- "Transmission Line Approach of
  Left-handed (LH) Materials", USNC/URSI 1 (2002) 39.
- A. K. Iyer, G. V. Eleftheriades -- "Negative Refractive Index
  Metamaterials Supporting 2-D Waves", IEEE MTT-S 2 (2002) 1067.

## Unit cell

The CRLH unit cell augments a conventional (right-handed) TL by
adding a *series capacitance C_L* and a *shunt inductance L_L*.

```
              L_R/2          L_R/2
   ----[||]--[~~~~]--+--[~~~~]--[||]----
       2C_L            2C_L
                  |
                 [~~~~] L_L
                  |
                === C_R
                  |
                 GND
```

- Series C_L + parasitic series inductance L_R -> LH mode at low f
- Shunt L_L + parasitic shunt capacitance C_R -> LH mode at low f
- Parasitics L_R, C_R dominate at high f -> RH mode

## Dispersion (balanced CRLH)

For a unit cell of length a, the propagation constant beta(omega)
satisfies (after Caloz-Itoh):

    cos(beta a) = 1 - (1/2) (omega L_L - 1/(omega C_R))(omega C_L - 1/(omega L_R))

LH region (omega < omega_se = 1/sqrt(L_L C_L)):
- beta < 0 (negative phase velocity)
- d omega / d beta > 0 (positive group velocity)
- backward wave (v_p . v_g < 0)

RH region (omega > omega_sh = 1/sqrt(L_R C_R)):
- standard forward wave

**Balanced CRLH**: omega_se = omega_sh -> seamless LH<->RH
transition with NO bandgap.  Practical engineering design hits
this exactly by tuning gap and stub dimensions.

## Applications (Sato 2005, Toyota CRDL)

| Device | Reference |
|--------|-----------|
| 0-dB tightly coupled directional coupler (35% BW) | Caloz et al., IEEE TMT 52 (2004) 980 |
| Zero-order resonator (size-independent f_0) | Sanada, APMC 3 (2003) 1588 |
| Left-handed leaky-wave antenna (LH-LWA) | Grbic, Eleftheriades, IEEE AP-S 4 (2003) 340 |
| 76 GHz mm-wave LH-LWA for ACC radar | Matsuzawa et al., IEEE Int. Workshop AT (2005) 183 |
| Via-free LH microstrip TL | Sanada et al., IEEE MTT-S (2004) 301 |
| Electronically-scanned LWA (varactor-loaded) | Lim et al., IEEE MTT-S 1 (2004) 313 |

## Backfire-to-endfire scanning

The LH-LWA scans its beam continuously from backfire (-90 deg)
through broadside (0 deg) to endfire (+90 deg) as frequency sweeps
through the LH-RH transition (balanced CRLH).  No bandgap means no
"blind angle" at broadside.  Toyota CRDL prototype achieved
70 deg -> 100 deg scan across 76-78 GHz for the 16-cell antenna.

## Lab connection

The CRLH unit cell maps directly onto a 2-port
`radia.peec.PEECBuilder` cell with series C_L and shunt L_L.  This
gives a Loop-Star ladder representation of the LH band and lets the
lab's PRIMA model-order-reduction tool extract a SPICE netlist for
system-level simulation -- something not possible with the resonant
SRR + wire approach.

The dispersion analysis is identical to the standard PEEC
filament-panel two-port impedance solved by
`radia.peec_topology.PEECCircuitSolver` at swept frequencies.

## References

- C. Caloz, T. Itoh, IEEE Trans. AP 52 (2004) 1159-1166.
- A. Sanada et al., IEEE Microwave Wireless Comp. Lett. 14 (2004)
  68-70.
- A. Sanada, IEEE Trans. MTT 52 (4) (2004) 1252-1263.
- K. Sato, "Research Trends in Metamaterials", R&D Review of Toyota
  CRDL 40 (2) (2005) 24-30.
- S. Matsuzawa et al., IEEE Int. Workshop Antenna Tech. (2005)
  183-186.
- C. Caloz, T. Itoh, *Electromagnetic Metamaterials: Transmission
  Line Theory and Microwave Applications*, Wiley-IEEE (2005).

## Cross-references

- `gyrator_lh_line` -- active LH TL (lab thesis)
- `radia_mcp.peec` -- PEEC ladder network
- `lh_materials_neg_index` -- the underlying physics
"""


GYRATOR_LH_LINE = r"""
# Gyrator-based left-handed transmission line (lab thesis)

## Background

Source: "ジャイレータを用いた左手系伝送線路に関する研究"
(Study of left-handed transmission lines using gyrators), lab
thesis archived at
public-safe curated corpus

The PDF is password-protected (lab-internal); this topic summarizes
the publicly known approach to gyrator-based LH TL design.

## Gyrator basics

A gyrator is a two-port linear element defined by its scattering
relation

    V_1 = + R_g * I_2,   V_2 = - R_g * I_1

It is non-reciprocal and "rotates" the port impedance:
attaching an impedance Z at port 2 makes port 1 appear as
R_g^2 / Z.  A gyrator therefore turns a *capacitance* at one port
into an *inductance* at the other port (and vice versa) -- which
is exactly what is needed to synthesize the LH unit cell (series
C_L, shunt L_L) entirely from passive C plus gyrators, without
discrete inductors.

## Construction of LH TL with gyrators

Tellegen's original gyrator proposal (1948) and the modern
realizations using operational-amplifier circuits (Antoniou, 1969;
Riordan, 1967) make it possible to build LH transmission lines
that:

- avoid lumped-inductor parasitics (no winding loss / no DCR)
- are reciprocal at the *line* level even when each gyrator is
  individually non-reciprocal (paired-gyrator design)
- can be tuned electronically (gyrator R_g set by op-amp bias)

The lab thesis investigates such an active LH TL implementation
and characterizes its dispersion, loss, and noise figure.

## Why this is interesting (vs CRLH)

| Property | CRLH | Gyrator LH |
|----------|------|------------|
| Passive | Yes | No (active op-amps) |
| Bandwidth | Wide | Tunable per gyrator R_g |
| Loss | Conductor + dielectric | Op-amp noise / supply |
| Non-reciprocity | Possible with ferrite | Built in (gyrator) |
| Integration | RF PCB / mm-wave | RF IC / CMOS |
| Quasi-static limit | OK | Excellent |

Active gyrator LH TL is particularly attractive in the kHz - MHz
range where physical inductors with high Q are large and lossy.
This overlaps with the lab's IH (1-500 kHz) and WPT (kHz-MHz)
target bands.

## Lab connection

The gyrator LH TL can be analyzed in `radia.peec` by including
controlled-source elements (VCCS, CCVS) representing the gyrator
in the PEEC ladder.  The same Loop-Star MNA assembly used for
passive PEEC accepts the gyrator stamp directly.

For papers / talks involving this lab thesis, please obtain the
password from the supervisor.

## References (open literature)

- B. D. H. Tellegen, "The Gyrator, a New Electric Network Element",
  Philips Res. Rep. 3 (1948) 81-101.
- A. Antoniou, "Realisation of gyrators using operational
  amplifiers, and their use in RC-active-network synthesis",
  Proc. IEE 116 (11) (1969) 1838-1850.
- R. H. S. Riordan, "Simulated inductors using differential
  amplifiers", Electron. Lett. 3 (2) (1967) 50-51.

## Cross-references

- `crlh_transmission_line` -- the passive alternative
- `radia_mcp.peec` -- PEEC ladder with controlled sources
"""


PERFECT_LENS_EVANESCENT = r"""
# Perfect lens and the evanescent-mode subtlety

## Pendry's perfect lens (2000)

J. B. Pendry, "Negative Refraction Makes a Perfect Lens",
Phys. Rev. Lett. 85 (2000) 3966-3969.

A slab of n = -1 material (eps = -1, mu = -1) of thickness d
placed between source and observation plane focuses BOTH the
propagating waves AND the evanescent waves -- the latter by
*resonant amplification* in the slab.  Result: sub-wavelength
imaging, theoretically without resolution limit.

Mechanism:
- Propagating wave: Snell-inverted refraction, two-pass refocus
- Evanescent wave: surface-plasmon-like amplification in the
  n = -1 slab restores the decaying field amplitude on the image
  side

## Transformation-optics lens (Leonhardt-Philbin 2006)

U. Leonhardt, T. G. Philbin, "General Relativity in Electrical
Engineering", New J. Phys. 8 (2006) 247.

A flat lens can also be derived from a folding coordinate map

    z(z_bar) = { z_bar         if z_bar < 0,
                 -alpha z_bar  if 0 < z_bar < D,
                 z_bar - (alpha + 1) D  if z_bar > D }

This three-valued map makes any source point map to three points
in laboratory space: source, mirror image inside the slab, image
behind the slab.  The medium inside the slab has eps = mu = -alpha
(negative-index region from orientation flip of the map).

## Bergamin's caveat (2009)

L. Bergamin, "Negative Index of Refraction, Spacetime Folding and
Perfect Imaging in Transformation Optics", arXiv:0912.2649.

Bergamin's three observations:

1. **Causality**: the transformation-optics lens forces an
   *instantaneous* update of the mirror image when the source
   changes -- so the device is strictly *stationary*.  Any
   time-varying signal violates causality.

2. **No imaging, just duplication**: a source at the source point
   AUTOMATICALLY creates a mirror sink inside the slab and a
   second source behind the slab.  The medium does not "image"
   the source -- it *triples* it.

3. **No evanescent amplification**: all evanescent fields on the
   image side of the slab are explained as the standard
   exponentially-decaying near-field of the *duplicated* source
   behind the slab.  There is *no* resonant amplification in the
   slab.  See Fig. 2 of Bergamin (2009) for the schematic.

Therefore the transformation-optics "perfect lens" and the
Pendry-Veselago "perfect lens" are NOT the same physical device,
despite having the same n = -1 slab description.  The former is a
clever stationary boundary-value construction; the latter is a
genuine evanescent-wave amplifier.

## Practical implication for the lab

If a lab project ever requires sub-wavelength imaging:

- Use Pendry-Veselago (resonant SRR + wire DNG) for the actual
  evanescent amplification.  Accept the narrow band and the loss.
- Do NOT design a transformation-optics flat lens and expect
  sub-wavelength imaging -- it will only work for static / slowly
  varying sources, and even then it triples sources rather than
  imaging.
- For sub-wavelength imaging of *static* magnetic fields (e.g. in
  scanning Hall microscopy, MFM), the Kelvin transform plus
  high-order BEM in `radia_mcp.radia_ngsolve` is a more honest
  route than a Pendry perfect lens.

## References

- J. B. Pendry, Phys. Rev. Lett. 85 (2000) 3966-3969.
- V. G. Veselago, Sov. Phys. Usp. 10 (1968) 509-514.
- U. Leonhardt, T. G. Philbin, New J. Phys. 8 (2006) 247.
- L. Bergamin, arXiv:0912.2649 (2009).
- U. Leonhardt, T. Philbin, "Perfect imaging with positive
  refraction in three dimensions", arXiv:0911.0522 (2009).
- U. Leonhardt, "Perfect imaging without negative refraction",
  New J. Phys. 11 (2009) 093040.
- S. A. Ramakrishna, "Physics of negative refractive index
  materials", Rep. Prog. Phys. 68 (2005) 449.

## Cross-references

- `lh_materials_neg_index` -- the physics
- `transformation_optics_kelvin` -- the related Kelvin map
"""


APPLICATIONS_AUTOMOTIVE = r"""
# Automotive and engineering applications

## Toyota CRDL roadmap (Sato 2005)

K. Sato, "Research Trends in Metamaterials", R&D Review of Toyota
CRDL Vol. 40 No. 2 (2005) 24-30.

The Toyota CRDL review identified three automotive use-cases that
remain the most cited industrial drivers:

| Use case | Device | Frequency |
|----------|--------|-----------|
| Adaptive Cruise Control (ACC) radar | mm-wave LH-LWA | 76-77 GHz |
| Mobile communication (5G mm-wave) | beam-steering antenna | 28 / 60 GHz |
| Electric vehicle motor | high-mu_r magnetic metamaterial | DC, kHz |
| EV electromagnetic compatibility | EBG absorber / shield | wideband |

## mm-wave LH-LWA for ACC radar (76-77 GHz)

Toyota CRDL prototype (Matsuzawa et al., IEEE Int. Workshop AT
2005, 183-186):

- 16-cell balanced CRLH unit cell
- Microstrip via-free implementation (gap capacitor + virtual
  ground patch)
- Unit cell parameters: L_L = 0.1 nH, L_R = 0.06 nH, C_L =
  0.07 pF, C_R = 0.8 pF, C_g = 0.045 pF
- Beam scan: 70 deg (backfire) to 100 deg (forward) across
  76 -> 78 GHz
- No broadside blind angle (balanced CRLH)
- Etching tolerance much better than interdigital-capacitor
  microwave designs

Extension: liquid-crystal-controlled CRLH for fixed-frequency
beam scanning (electronic ACC sweep without VCO).

## EV motor magnetic materials

The Toyota CRDL roadmap (Sato 2005) flagged "novel high-performance
magnetic materials for electric motors" as a strategic
metamaterial application.  Concept:

- Sub-mm periodic composite of soft magnetic + insulating phases
- Eddy-current loss reduced by sub-skin-depth phase size
- Effective mu_r tunable by composite microstructure
- TO-designed lamination patterns

This overlaps directly with the lab's IH (induction-heating) and
ESIM workflow (`radia_mcp.ih`):

- IH uses ESIM to compute the *effective surface impedance* of a
  nonlinear iron workpiece
- Motor metamaterial uses the *effective bulk* mu and sigma of a
  laminated stack
- Same homogenization mathematics, different boundary conditions

## EMC absorber / shield

Metamaterial absorbers (Engheta 2002; Landy et al. PRL 100, 207402,
2008) provide near-unity absorption at chosen frequencies with
sub-wavelength thickness, replacing the bulky carbon-loaded foams
used for EV EMC compliance.  Adjacent topic: high-impedance
surfaces (Sievenpiper, IEEE TMT 47 (1999) 2059) for antenna
ground planes.

## EBG substrate antennas

Sanada (IEEE Trans. MTT 52 (4) (2004) 1252) showed planar EBG
structures with negative refractive index for sub-wavelength
focusing lenses, demonstrated at C-band.  Production EV roof
antennas use mushroom-EBG ground planes for low-profile vertical
patches.

## Lab connection

- `radia_mcp.ih` ESIM is the closed-form workpiece-side analogue
  of the EV-motor metamaterial homogenization
- `radia_mcp.peec` directly synthesizes the CRLH cell as a 2-port
  PEEC ladder (model-order-reducible to a SPICE-friendly netlist
  via PRIMA)
- `radia_mcp.radia_ngsolve.kelvin_periodic_identifications`
  provides the open-boundary FEM tool that EV antenna designers
  need for near-field radiation analysis without PML

## References

- K. Sato, R&D Review of Toyota CRDL 40 (2) (2005) 24-30.
- S. Matsuzawa et al., IEEE Int. Workshop Antenna Tech. (2005)
  183-186.
- N. I. Landy et al., Phys. Rev. Lett. 100 (2008) 207402.
- N. Engheta, IEEE AP-S Int. Symp. Dig. 2 (2002) 392-395.
- D. Sievenpiper, IEEE Trans. MTT 47 (1999) 2059-2074.
- A. Sanada, IEEE Trans. MTT 52 (4) (2004) 1252-1263.
- F. Yang, Y. Rahmat-Samii, Microwave Opt. Technol. Lett. 31 (3)
  (2001) 165-168 (low-profile EBG antenna).
"""


ACOUSTIC_SONIC_CRYSTAL = r"""
# Acoustic metamaterials / sonic crystals: band gaps and validation discipline

The acoustic counterpart of the EM metamaterial family: a *sonic
crystal* is a periodic array of scatterers in a fluid (air/water)
engineered for stop bands, negative effective density/modulus, lensing
or absorption; *phononic crystals* are the elastic-wave siblings. The
same homogenization mathematics as `effective_medium` applies (cell
problems -> effective parameters), and the same unit-cell discipline as
EM band-gap design.

## The load-bearing finding: confinement is NOT optional

A sparse FREE-SPACE chain of strong scatterers does NOT develop a Bragg
stop band at k d = pi. Verified four independent ways on a chain of
sound-soft spheres (R/d = 0.2, five periods): Foldy monopole multiple
scattering, a Galerkin BEM, a second independent BEM implementation, and
a volume FEM with a radiation ABC all agree that the insertion loss is a
FLAT broadband plateau (sub-wavelength attenuation: a sound-soft sphere
keeps scattering length ~ R down to k -> 0) with no feature at the Bragg
wavenumber. In open 3D the energy leaks sideways; the 1D interference
mechanism never accumulates.

The band gap appears the moment the SAME scatterer chain is CONFINED:
in a rigid duct (cross-section a x a, single-mode below the transverse
cutoff pi/a) the Bloch band structure of the unit cell and the
transmission through a finite N-cell crystal show the full
stop / pass / stop structure, aligned quantitatively:

- below band 1: soft (pressure-release) inclusions cut off the
  long-wavelength plane mode entirely (T ~ 1e-4),
- inside band 1: a passband (finite-crystal transmission with
  Fabry-Perot ripple),
- inside the Bragg gap above band 1: T drops by ~2 orders of magnitude.

Locally resonant variants (concentric C-shaped shells) push the gap far
below the Bragg frequency: Elford, Chalmers, Kusmartsev, Sambles,
"Matryoshka locally resonant sonic crystal", J. Acoust. Soc. Am. 130,
2746 (2011).

## Methodology (tool-agnostic recipe)

1. **Unit cell (Bloch)**: Helmholtz eigenproblem on one cell, Floquet
   phase-periodic BC exp(-i q d) on the lattice faces, Dirichlet (soft)
   or Neumann (rigid) inclusion, complex FE space. GATE before physics:
   the EMPTY cell's bands must reproduce the analytic duct dispersion
   k = sqrt((m1^2+m2^2) pi^2/a^2 + ((q d + 2 pi n)/d)^2)
   exactly (measured 2e-4 with order-2 FEM) - this catches
   phase-convention errors AND boundary-classification bugs (a
   misclassified face silently turns the duct into a Dirichlet cavity;
   classify faces geometrically, e.g. by whether the face center lies on
   a wall plane, never by area thresholds).
2. **Finite crystal (transmission)**: duct of N cells + pads, one-way
   port Robin conditions (inlet dp/dn - i k p = -2 i k for a unit
   incident plane mode, outlet dp/dn - i k p = 0), exact below the
   transverse cutoff. GATE: the EMPTY duct must be transparent
   (|T - 1| < 1e-2 measured at order 2). Transmission dips must sit
   INSIDE the Bloch gap.
3. **Effective medium**: for sub-wavelength cells the chain/crystal
   behaves as a modified compressibility/density medium - the cell
   problem is the same mathematics as the EM `effective_medium` topic
   and the eddy-current lamination homogenization used in motor
   modeling.

## Multiple-scattering reference (analytic class)

Foldy (1945): self-consistent point-monopole system. For sound-soft
spheres the exact monopole coefficient is f = -j_0(kR)/h_0(kR); exciting
fields solve (I - f H) E = p_inc with H_ij = h_0(k |x_i - x_j|). Its
validity envelope is LOW-k only (the neglected l >= 1 single-sphere
terms grow from ~4% at kR = 0.18 to ~47% at kR = 0.9): MEASURE the
envelope against the exact single-sphere series and lock the degradation
itself, never assume the model.

## BEM validation discipline for scatterer-level acoustic solvers

- **Kernel conventions are measured, never assumed.** Unit-sphere
  eigenvalues pin everything: V[Y_l] = i k j_l(k) h_l(k),
  K[Y_l] = 1/2 + i k^2 j_l(k) h_l'(k) (outward-normal PV convention);
  the imaginary signs distinguish e^{+ikr} from e^{-ikr}, and an
  operator must match a trusted reference far better than its
  conjugate.
- **The exterior Dirichlet EXACT gate**: boundary data from an interior
  point source must be reproduced at every exterior point (uniqueness +
  radiation condition) - zero series truncation, pure discretization
  error. Generalizes verbatim to multi-body (source inside any one
  scatterer).
- **Transmission (penetrable) solvers**: the acoustic INVISIBILITY null
  gate - set k1 = k0 and rho1 = rho0 and the scatterer must vanish;
  the residual is the pure discretization null error and must CONVERGE
  under refinement. Real contrast is gated by the Anderson (1950)
  fluid-sphere partial-wave series; solve its per-mode 2x2 transmission
  match by ANALYTIC elimination through the interior log-derivative
  (k1/rho1) j_l'(x1)/j_l(x1) - the naive 2x2 is catastrophically
  ill-conditioned at high l (h_l huge, j_l tiny) and quietly pollutes
  the invisible case. Size the series term count on the FARTHEST
  evaluation point (k r_max), not the scatterer radius.
- **Irregular frequencies**: the first-kind single-layer equation and
  the total-field double-layer (rigid) equation are singular at the
  interior Dirichlet eigenvalues of the surface (unit sphere: first at
  kR = pi). The discrete condition number can look BENIGN (the faceted
  eigenvalue shifts) while the solution is ~100% wrong - only analytic
  gates catch it. Classic fixes: CHIEF (Schenck 1968; append interior
  null-field rows at a few jittered interior points, least squares -
  restores regular-k accuracy) for teaching, Burton-Miller (1971;
  combined-field, needs the hypersingular operator) for production.
- **FEM/BEM coupling needs no ABC**: the BIE row IS the exact radiation
  condition. A volume-FEM leg with a first-order Sommerfeld ABC on a
  sphere at kR ~ 6-12 sits ~4-8% from the BEM answer at exterior
  probes - useful as a methods-diverse cross-check, with the ABC error
  quantified rather than argued.

## Cross-links

- `effective_medium` - the EM homogenization sibling (same cell-problem
  mathematics).
- `radia_mcp.radia_ngsolve` - periodic FE spaces (Bloch phases), BEM.
- `radia_mcp.ih` - ESIM 1D cell problems (effective surface impedance:
  the same "solve a cell, extract an effective parameter" pattern).
- `radia_mcp.motor` - Hollaus lamination homogenization.

## References

- L. L. Foldy, "The multiple scattering of waves", Phys. Rev. 67, 107
  (1945).
- V. C. Anderson, "Sound scattering from a fluid sphere", J. Acoust.
  Soc. Am. 22, 426 (1950).
- H. A. Schenck, "Improved integral formulation for acoustic radiation
  problems" (CHIEF), J. Acoust. Soc. Am. 44, 41 (1968).
- A. J. Burton, G. F. Miller, "The application of integral equation
  methods to the numerical solution of some exterior boundary-value
  problems", Proc. R. Soc. Lond. A 323, 201 (1971).
- D. P. Elford et al., "Matryoshka locally resonant sonic crystal",
  J. Acoust. Soc. Am. 130, 2746 (2011).
- P. A. Martin, "Multiple Scattering: Interaction of Time-Harmonic
  Waves with N Obstacles", Cambridge Univ. Press (2006).
"""


# ============================================================
# Dispatcher
# ============================================================


def get_knowledge(topic: str = "overview") -> str:
    """Dispatch metamaterial topics.

    See TOPICS for the authoritative enum and one-line summaries.
    """
    t = topic.lower().strip()
    if t in ("overview", "intro", ""):
        return OVERVIEW
    if t in ("lh_materials_neg_index", "lh", "lefthanded", "left_handed",
             "negative_index", "nir", "dng", "vesalago", "veselago"):
        return LH_MATERIALS_NEG_INDEX
    if t in ("srr_split_ring", "srr", "split_ring", "pendry_srr"):
        return SRR_SPLIT_RING
    if t in ("effective_medium", "homogenization", "homogenisation",
             "wire_medium", "retrieval"):
        return EFFECTIVE_MEDIUM
    if t in ("transformation_optics_kelvin", "transformation_optics",
             "to", "kelvin", "cloaking", "pendry_cloak", "leonhardt"):
        return TRANSFORMATION_OPTICS_KELVIN
    if t in ("magneto_optical", "mom", "faraday", "biyig", "bi_yig",
             "magneto_optic"):
        return MAGNETO_OPTICAL
    if t in ("topology_optimization_design", "topopt", "topology_opt",
             "topology_optimization", "inverse_design"):
        return TOPOLOGY_OPTIMIZATION_DESIGN
    if t in ("crlh_transmission_line", "crlh", "crlh_tl",
             "transmission_line", "leaky_wave_antenna", "lwa"):
        return CRLH_TRANSMISSION_LINE
    if t in ("gyrator_lh_line", "gyrator", "active_lh_line"):
        return GYRATOR_LH_LINE
    if t in ("perfect_lens_evanescent", "perfect_lens", "pendry_lens",
             "veselago_lens", "evanescent_amplification"):
        return PERFECT_LENS_EVANESCENT
    if t in ("applications_automotive", "automotive", "applications",
             "toyota", "ev_motor"):
        return APPLICATIONS_AUTOMOTIVE
    if t in ("acoustic_sonic_crystal", "acoustic", "acoustic_metamaterial",
             "sonic_crystal", "phononic", "phononic_crystal", "band_gap",
             "bloch", "chief", "burton_miller"):
        return ACOUSTIC_SONIC_CRYSTAL
    if t == "all":
        return "\n\n".join([
            OVERVIEW,
            LH_MATERIALS_NEG_INDEX,
            SRR_SPLIT_RING,
            EFFECTIVE_MEDIUM,
            TRANSFORMATION_OPTICS_KELVIN,
            MAGNETO_OPTICAL,
            TOPOLOGY_OPTIMIZATION_DESIGN,
            CRLH_TRANSMISSION_LINE,
            GYRATOR_LH_LINE,
            PERFECT_LENS_EVANESCENT,
            APPLICATIONS_AUTOMOTIVE,
            ACOUSTIC_SONIC_CRYSTAL,
        ])
    return (
        f"Unknown topic '{topic}'. Available topics: "
        + ", ".join(TOPICS.keys())
    )
