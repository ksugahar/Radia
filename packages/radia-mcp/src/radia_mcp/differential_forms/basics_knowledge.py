"""
Basic differential-form machinery: tangent vectors, k-forms, wedge product,
exterior derivative, Stokes' theorem.

Distilled from:
  - 微分形式.pdf  (sections 1-5, 8; tangent/cotangent spaces, tensor product,
                  k-forms via alternating tensor product, wedge product,
                  exterior derivative, Stokes' theorem)
  - Arnold-Falk-Winther 2006, sec. 2 (exterior algebra and exterior calculus
                                       on manifolds)
  - Bossavit 1998 App. A (mathematical background)
"""

TANGENT_VECTORS = """
# Tangent vectors and the dual (cotangent) space

## Curve and tangent vector
A curve on an n-dimensional space  c(t)  has tangent vector  v = dc/dt.
At a point  p, the set of all such tangent vectors forms the
**tangent space**  T_p (an n-dimensional vector space).

## Local coordinates and basis
Pick n curves  c_1(t), ..., c_n(t)  through  p, mark each with monotonically
increasing parameter  X_i(x_1, ..., x_n).  Then

    e_i := ∂X / ∂x_i          (contravariant basis, "ket": |∂/∂x_i⟩)
    ⟨e^i| := ∂x_i / ∂X        (covariant basis, "bra": ⟨dx_i|)

Duality:
    ⟨dx_i | ∂/∂x_j⟩ = δ_i^j
    ∑_i |∂/∂x_i⟩⟨dx_i| = 1

(Bra-ket notation used in 微分形式.pdf §1; this is the same as the
standard "contraction" pairing of T_p with T_p^*.)

## Gradient as a 1-form
For a scalar function f(x_1, ..., x_n):

    ⟨df| = ∑_i (∂f/∂x_i) ⟨dx_i|

— the differential of  f. Note that this is a covariant object (1-form),
not a vector. Identifying it with grad f requires a metric (Riesz isomorphism).

## Why this matters in EM
Vector calculus glosses over the distinction between vectors and covectors
because we use Cartesian coordinates with the Euclidean metric.  In a
general coordinate system (cylindrical/spherical/curvilinear), or on a
curved manifold, the distinction is unavoidable.  Differential forms
make the right object explicit:

- E (electric field) is a 1-form (acts on lines)
- B (magnetic flux density) is a 2-form (acts on surfaces)
- J (current density) is a 2-form
- ρ (charge density) is a 3-form

See `maxwell` topic for the explicit construction.
"""

WEDGE_PRODUCT = """
# Tensor product, alternating product, wedge product

## Tensor product of vectors
From n k-tuples of vectors  v_{i_1}, ..., v_{i_k}  build

    |v_{i_1} ⊗ v_{i_2} ⊗ ... ⊗ v_{i_k}⟩

Spans an  n^k-dimensional space.  The pairing
⟨u_{i_1} ⊗ ... ⊗ u_{i_k} | v_{j_1} ⊗ ... ⊗ v_{j_k}⟩
= ⟨u_{i_1}|v_{j_1}⟩ · ⟨u_{i_2}|v_{j_2}⟩ · ... · ⟨u_{i_k}|v_{j_k}⟩.

## Alternating (antisymmetric) part
Restrict to the **antisymmetric subspace** under permutation of indices:

    |v_{i_1} ∧ v_{i_2} ∧ ... ∧ v_{i_k}⟩ :=
        ∑_σ sign(σ) |v_{σ(i_1)} ⊗ ... ⊗ v_{σ(i_k)}⟩

This subspace has dimension  C(n,k) = n! / (k!(n-k)!).  Elements are
**k-vectors** (Grassmann algebra). Their duals (acting on k-tuples of
vectors and returning a real number) are **k-forms**.

## Wedge product (exterior product)
For ω a k-form and η an l-form:

    (ω ∧ η)(v_1, ..., v_{k+l})
        = (1 / (k! l!)) ∑_σ sign(σ) ω(v_{σ(1)}, ..., v_{σ(k)})
                                    · η(v_{σ(k+1)}, ..., v_{σ(k+l)})

(Arnold-Falk-Winther uses an alternative normalization that drops the
1/(k!l!) by restricting σ to ordered subsets.)

## Properties
- Bilinear:   (aω₁ + bω₂) ∧ η = a ω₁∧η + b ω₂∧η
- Associative: (ω ∧ η) ∧ γ = ω ∧ (η ∧ γ)
- Graded-anticommutative:  η ∧ ω = (-1)^{kl} ω ∧ η
- In particular   dx^i ∧ dx^i = 0  (so 1-form squared is zero).

## Dimension count
- 0-form: 1 component (scalar)
- 1-form: n components  (covector)
- 2-form: n(n-1)/2 components
- ...
- n-form: 1 component (top form / volume form)
- k-form with k > n: zero.

In E^3: 0-, 1-, 2-, 3-forms have 1, 3, 3, 1 components respectively.
The 1↔2-form and 0↔3-form symmetry is the **Hodge duality** — see `maxwell`.
"""

EXTERIOR_DERIVATIVE = """
# Exterior derivative d

## Definition
For a k-form  ω = ∑ f_{i_1,...,i_k} ⟨dx^{i_1} ∧ ... ∧ dx^{i_k}|  the
exterior derivative is a (k+1)-form:

    dω = ∑ ⟨df_{i_1,...,i_k}| ∧ ⟨dx^{i_1} ∧ ... ∧ dx^{i_k}|

where the gradient of each component coefficient is taken.

## Examples in 3D
**0-form  f**:
    df = (∂f/∂x_1) dx_1 + (∂f/∂x_2) dx_2 + (∂f/∂x_3) dx_3
                                                            ↔  grad f

**1-form  ω = f_1 dx_1 + f_2 dx_2 + f_3 dx_3**:
    dω = (∂f_3/∂x_2 − ∂f_2/∂x_3) dx_2∧dx_3
       + (∂f_1/∂x_3 − ∂f_3/∂x_1) dx_3∧dx_1
       + (∂f_2/∂x_1 − ∂f_1/∂x_2) dx_1∧dx_2          ↔  rot ω

**2-form  ω = f_1 dx_2∧dx_3 + f_2 dx_3∧dx_1 + f_3 dx_1∧dx_2**:
    dω = (∂f_1/∂x_1 + ∂f_2/∂x_2 + ∂f_3/∂x_3) dx_1∧dx_2∧dx_3   ↔  div ω

So in 3D the operators grad, rot, div are **the same operator d**, only
acting on forms of different degree.

## Properties
- Linearity:  d(aω + bη) = a dω + b dη
- Graded Leibniz:  d(ω ∧ η) = dω ∧ η + (−1)^k ω ∧ dη
- **Nilpotent**:  d(dω) = 0  for any ω.
   In 3D this becomes  rot(grad f) = 0  and  div(rot v) = 0.

## Poincaré lemma
On a **contractible** domain, the converse holds:
    dω = 0  ⇒  ∃ η  such that  ω = dη

so closed forms (dω = 0) are exact (ω = dη). In 3D:
- curl-free field is a gradient
- divergence-free field is a curl
- divergence-free scalar (degenerate) is the volume integral of something.

On non-contractible domains, the **de Rham cohomology**
H^k_dR := ker(d : Λ^k) / im(d : Λ^{k-1}) is non-trivial — its
dimension is the Betti number  b_k.  See `homology` topic.

## Stokes' theorem (the unifier)
For a k-form ω on an oriented (k+1)-chain σ with boundary ∂σ:

    ∫_σ dω = ∫_{∂σ} ω

In 3D:
- k=0: fundamental theorem of calculus
- k=1: classical Stokes ∫_S rot a · dS = ∮_∂S a · dl
- k=2: divergence theorem ∫_V div b dV = ∫_∂V b · dS
"""

STOKES_THEOREM = """
# Stokes' theorem — the master integration formula

## The single statement
For an oriented (k+1)-chain σ in an n-manifold and a k-form ω:

    ∫_σ dω = ∫_{∂σ} ω.

Discrete chain analogue (with Whitney elements):
    ⟨dω, c⟩ = ⟨ω, ∂c⟩
where ⟨·,·⟩ is the integration pairing (1-form · 1-chain → real).

## Recovery of classical theorems
**3D, k=0 (1-chain σ = curve from A to B, ω = f ∈ Λ^0):**
    ∫_A^B df = f(B) − f(A).
    Fundamental theorem of calculus.

**3D, k=1 (2-chain σ = surface S, ω = 1-form a corresponding to vector A):**
    ∫_S rot A · dS = ∮_∂S A · dl.
    Classical Stokes/Kelvin theorem.

**3D, k=2 (3-chain σ = volume V, ω = 2-form b corresponding to vector B):**
    ∫_V div B dV = ∫_∂V B · n dS.
    Gauss-Ostrogradsky divergence theorem.

## Why this matters for FEM
The discrete Whitney complex (next: `whitney`) has the matrix relations

    D R = 0,        R G = 0

(D, R, G the cell-face/face-edge/edge-node incidence matrices) as the
exact discrete analogue of  div(rot) = 0,  rot(grad) = 0.  This is **not
approximate** — it is the discrete Stokes theorem, structurally embedded
in the mesh.  See `whitney` and `homology` for details.

## Twisted vs straight integration
For "twisted" forms (densities, currents), the integration formula is
the same but the orientation of the integration domain is **outer**
(crossing-direction) rather than **inner** (tangential).
"""


def get_basics_documentation(topic: str = "all") -> str:
    """Dispatch by topic keyword.

    Topics:
      "all"                    — concatenate everything
      "tangent_vectors"        — bra-ket basis, gradient as 1-form
      "wedge_product"          — alternating tensor product, ∧
      "exterior_derivative"    — d, properties, Poincaré lemma
      "stokes"                 — Stokes' theorem and its classical instances
    """
    topic = topic.lower().strip()
    if topic == "tangent_vectors":
        return TANGENT_VECTORS
    if topic in ("wedge_product", "wedge"):
        return WEDGE_PRODUCT
    if topic in ("exterior_derivative", "d", "exterior"):
        return EXTERIOR_DERIVATIVE
    if topic in ("stokes", "stokes_theorem"):
        return STOKES_THEOREM
    if topic == "all":
        return "\n\n".join([
            TANGENT_VECTORS,
            WEDGE_PRODUCT,
            EXTERIOR_DERIVATIVE,
            STOKES_THEOREM,
        ])
    return (
        f"Unknown topic '{topic}'. Available: "
        "all, tangent_vectors, wedge_product, exterior_derivative, stokes."
    )
