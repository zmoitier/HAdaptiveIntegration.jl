---
title: '`HAdaptiveIntegration.jl`: Adaptive numerical integration over simplices and orthotopes'
tags:
  - Julia
  - Scientific computing
  - Numerical integration
authors:
  - name: Luiz Faria
    orcid: 0000-0002-8090-934X
    affiliation: 1
  - name: Zoïs Moitier
    corresponding: true
    orcid: 0000-0003-1736-5754
    affiliation: 2
affiliations:
  - name: POEMS, CNRS, Inria, ENSTA, Institut Polytechnique de Paris, 91120 Palaiseau, France
    index: 1
  - name: Inria, Unité de Mathématiques Appliquées, ENSTA, Institut Polytechnique de Paris, 91120 Palaiseau, France
    index: 2
date: 19 March 2026
bibliography: paper.bib
---

# Summary

<!--
A description of the high-level functionality and purpose of the software for a diverse, 
non-specialist audience.
-->

`HAdaptiveIntegration` is a Julia package for adaptive numerical integration on
multidimensional simplices and orthotopes. It approximates integrals of the form
$$
  I = \int_{\Omega} f(\boldsymbol{x}) \, \mathrm{d}\boldsymbol{x}
$$
where

- $f \colon \mathbb{R}^d \to \mathbb{F}$ is any Julia function mapping $d$-dimensional
  vectors to elements of type $\mathbb{F}$ supporting multiplication by a real scalar,
  addition and a norm (*i.e.* a normed real vector space);
- $\Omega \subset \mathbb{R}^d$ is the integration domain: simplices and orthotopes are
  supported.

The package targets adaptive cubature in low to moderate dimension, with a particular
emphasis on domains that arise naturally in mesh-based scientific computing. 
`HAdaptiveIntegration` combines uniform subdivision rules with embedded cubature, so that it
returns both an integral estimate and an a posteriori error estimate.

Its main features are:

- adaptive integration over simplices and orthotopes of arbitrary dimension;
- efficient tabulated cubature rules for low-dimensional simplices and orthotopes;
- support for user-defined embedded cubature rules and subdivision strategies;
- arbitrary-precision arithmetic.

## Usage

The package is centered around a single function, `integrate(f, domain)`, together with
constructors for common domains.

**Create a domain.**
There are seven domain constructors:

- `Segment`: defined by its two end points;
- `Triangle`, `Tetrahedron`, `Simplex`: defined by its vertices;
- `Rectangle`, `Cuboid`, `Orthotope`: defined by its low and high corner.

For example:
```julia
using HAdaptiveIntegration

segment     = Segment(0, 1)
triangle    = Triangle((0,0), (1,0), (0,1))
rectangle   = Rectangle((0,0), (1,1))
```

**Integrate.**
Once the domain is defined, the same interface is used for simplices and orthotopes:

```julia
f(x) = cis(sum(x)) / (sum(abs2, x) + 1e-2)

I, E = integrate(f, triangle)
I, E = integrate(f, rectangle)
```

The return value is a pair `(I, E)`, where `I` is the integral estimate and `E` is an a
posteriori error estimate.

# Statement of need

<!--
A section that clearly illustrates the research purpose of the software and places it in the context of related work.
This should clearly state what problems the software is designed to solve, who the target audience is, and its relation to other work.
-->

Adaptive numerical integration is a routine building block in scientific computing,
including finite-element and boundary-element methods, validation of semi-analytical
formulae, and parameter studies where repeated integrals must be computed over reference or
physical cells. In Julia, existing packages already cover important parts of this space:
`QuadGK.jl` is the natural tool for one-dimensional adaptive integration [@QuadGK],
`HCubature.jl` provides adaptive cubature on axis-aligned orthotopes [@HCubature], and
`Cuba.jl` exposes a family of multidimensional integration methods, including stochastic
ones, for rectangular domains [@Cuba].

The missing piece addressed by `HAdaptiveIntegration` is a Julia-native integrator that
handles simplices and orthotopes through a common interface, including triangles,
tetrahedra, higher-dimensional simplices, and low-dimensional orthotopes with specialized
tabulated rules. This is particularly useful for researchers working with simplicial meshes,
where mapping everything to boxes or maintaining separate code paths is inconvenient and can
obscure geometric intent. The package targets users who need deterministic adaptive
integration with explicit error estimates, access to custom cubature and subdivision
strategies, and the option to move beyond `Float64` arithmetic when validation or sensitive
applications require it.

# State of the field

<!--
A description of how this software compares to other commonly-used packages in the research area.
If related tools exist, provide a clear “build vs. contribute” justification explaining your unique scholarly contribution and why existing alternatives are insufficient.
-->

`HAdaptiveIntegration` is best viewed as complementary to the existing Julia ecosystem
rather than as a wholesale replacement. `QuadGK.jl` remains the appropriate choice in one
dimension because it is specialized for segments and offers very mature Gauss-Kronrod-based
integration [@QuadGK]. `HCubature.jl` is the closest comparison for multidimensional box
domains and is often the better option for medium- and high-dimensional orthotopes
[@HCubature]. `Cuba.jl` remains attractive for high-dimensional rectangular problems where
stochastic or importance-sampling methods are more suitable than deterministic refinement
[@Cuba].

The distinct contribution of `HAdaptiveIntegration` is its support for simplices of
arbitrary dimension together with orthotopes under one API, and its use of efficient
tabulated rules in low-dimensional cases where fixed-node formulas are especially effective.
That combination is not provided by the packages above. A standalone package is therefore
justified: domain geometry influences the reference mapping, admissible default cubature
rules, subdivision schemes, and the practical user workflow. Exposing these choices in a
single package makes it easier to write research code that works uniformly on triangles,
tetrahedra, rectangles, cuboids, and their higher-dimensional counterparts.

# Software design

<!--
An explanation of the trade-offs you weighed, the design/architecture you chose, and why it matters for your research application.
This should demonstrate meaningful design thinking beyond a superficial code structure description.
-->

## Embedded cubatures

Each local estimate is produced by an embedded cubature pair $(Q_h, Q_l)$ defined on a
reference domain. The higher-order rule $Q_h$ provides the integral estimate on the current
subdomain, while the difference between the high- and low-order rules yields a local error
indicator,
$$
E_K = \lVert Q_h(f; K) - Q_l(f; K) \rVert,
$$
after mapping the reference nodes to the physical subdomain $K$. This design has two
practical benefits: it provides an a posteriori error estimate without requiring derivative
information, and it lets the package switch rules depending on geometry and dimension.

The default rules are chosen to match the domain type. Segments use Gauss-Kronrod pairs,
triangles use a Radon-Laurie embedded rule, generic simplices use Grundmann-Moeller-based
pairs, generic orthotopes use Genz-Malik rules, and low-dimensional rectangles and cuboids
use tabulated formulas. Users can replace these defaults by supplying their own
`EmbeddedCubature`, which keeps the adaptive driver independent from any single cubature
family.

## Adaptive algorithm

The adaptive driver maintains a collection of active subdomains in a max-heap keyed by the
local error estimate. Starting from the initial domain, the algorithm repeatedly removes the
subdomain with the largest estimated error, subtracts its previous contribution from the
running totals, subdivides it, evaluates the children, and inserts them back into the heap.
This greedy strategy concentrates work where the current error indicator is largest while
keeping the global estimate available after every refinement.

Subdivision follows the geometry of the domain. Orthotopes are bisected along each axis,
producing $2^d$ children, while simplices are subdivided by midpoint-based simplicial
refinement, again producing $2^d$ children in dimension $d$. Integration stops when the
global error estimate satisfies either an absolute tolerance or a relative tolerance,
$$
E \leq \mathtt{atol} \quad \text{or} \quad E \leq \mathtt{rtol}\,\lVert I \rVert,
$$
or when a user-specified maximum number of subdivisions is reached. Two implementation
details matter in practice: reusable heap buffers allow allocation-free repeated calls, and
an optional callback exposes convergence information for diagnostics or benchmarking.

# Research impact statement

<!--
Evidence of realized impact (publications, external use, integrations) or credible near-term significance (benchmarks, reproducible materials, community-readiness signals).
The evidence should be compelling and specific, not aspirational.
-->

`HAdaptiveIntegration` is ready for direct use in research code. The repository includes
API documentation, advanced examples covering custom cubature and subdivision strategies, an
extension for arbitrary-precision workflows, and an automated test suite spanning segments,
triangles, rectangles, tetrahedra, cuboids, and four-dimensional simplices and orthotopes.
Continuous integration checks tests, documentation, and linting, which lowers the barrier
to adoption in reproducible computational projects.

Its near-term impact is strongest in workflows that already manipulate simplicial cells,
such as finite-element, boundary-element, and high-order discretization codes, because those
applications benefit from integrating directly on triangles and tetrahedra instead of first
rewriting problems on rectangular boxes. The package also serves validation and sensitivity
studies through its support for arbitrary-precision arithmetic, which is uncommon among
lightweight adaptive cubature packages.

# Complexity estimates

For a fixed dimension $d$, let $H$ be the number of nodes in the high-order rule and let
$N$ be the number of active subdomains. Evaluating one local embedded cubature requires $H$
function evaluations, while selecting the next subdomain to refine costs $O(\log N)$ due to
the max-heap. Since one refinement creates $2^d$ children for the generic simplex and
orthotope strategies used here, the computational cost per refinement is dominated by the
evaluation of $2^d$ new local cubatures, plus heap maintenance. In practice, function
evaluation usually dominates the bookkeeping overhead.

## Uniform regularity

For sufficiently smooth integrands, a standard heuristic is that the local quadrature error
on a cell of diameter $h$ behaves like $C h^p$, where $p$ reflects the effective order of
the embedded estimator. Under uniform isotropic refinement, $h$ is halved at each level and
the number of cells grows like $N \sim 2^{d\ell}$ after $\ell$ refinement levels. This leads
to the algebraic estimate
$$
E(N) \approx C N^{-p/d},
$$
so reaching a tolerance $\varepsilon$ requires $N = O(\varepsilon^{-d/p})$ cells in the
uniform regime. The adaptive algorithm does not improve the asymptotic order for globally
smooth problems, but it can reduce the constant by refining only the subdomains identified
as most relevant by the error estimator.

## Isotropic (nearly-)singularities

When the integrand is smooth on most of the domain but has a localized singularity or sharp
feature, uniform refinement spends too many function evaluations on regions that are already
well resolved. The heap-based adaptive strategy instead concentrates subdivisions near the
small subset of cells that dominate the error. Heuristically, the complexity then approaches
that of resolving the singular neighborhood plus a coarse background mesh, rather than the
cost of uniformly refining the entire domain. This is one of the main practical reasons to
prefer adaptive cubature for nearly singular kernels and corner singularities.

# Benchmarks

Representative comparisons against `HCubature.jl` on orthotopes confirm the intended scope
of the package. On a smooth integral over the unit square,
$f(x) = \exp(x_1 + x_2)$ with `rtol = 1e-8`, `HAdaptiveIntegration` required 125 function
evaluations versus 493 for `HCubature.jl`, and the observed runtime on a local Linux run was
about 1.1 microseconds versus 3.8 microseconds. On a more difficult corner-singular square
test,
$f(x) = (x_1^2 + x_2^2 + 10^{-8})^{-1/2}$ with `rtol = 1e-6`, the counts were 2125 versus
2941 evaluations, again in favor of `HAdaptiveIntegration`.

The same trend appears in three dimensions for low-dimensional orthotopes with specialized
rules: on the unit cuboid with $f(x) = \exp(x_1 + x_2 + x_3)$ and `rtol = 1e-6`,
`HAdaptiveIntegration` used 65 evaluations versus 429 for `HCubature.jl`. However, the
crossover expected from the state-of-the-field discussion is also visible: on a smooth
four-dimensional orthotope with `rtol = 1e-4`, `HCubature.jl` used fewer evaluations (171
versus 969) and was faster on the same machine. These measurements support the intended
positioning of `HAdaptiveIntegration`: low-dimensional simplices and orthotopes, especially
when geometry-specific tabulated rules are available, while `HCubature.jl` remains a strong
choice for higher-dimensional boxes.

# Extended precision

Low-dimensional tabulated cubatures are stored with finite decimal precision. For demanding
accuracy targets, this tabulation error can become the limiting factor even if the integrand
itself is evaluated in `BigFloat` arithmetic. `HAdaptiveIntegration` addresses this with the
optional `IncreasePrecisionExt` extension, which reconstructs higher-precision tabulated
rules by solving the polynomial exactness conditions with Newton iterations and automatic
differentiation through `ForwardDiff`.

This design keeps the default package lightweight while enabling high-precision workflows
when needed. In the test suite, deliberately reduced-precision tabulated rules fail to reach
relative accuracies around $10^{-14}$ after conversion to `Float64`, whereas their
reconstructed higher-precision counterparts recover the requested accuracy. The same
mechanism can then be used together with `BigFloat` domains and functions for validation
studies or sensitive numerical experiments.

# AI usage disclosure

Generative AI was used to help draft and edit parts of this manuscript. All text produced
with AI assistance was reviewed, revised, and verified by the authors against the source
code, tests, and benchmark outputs before inclusion in the paper.

# Acknowledgements

The authors thank the maintainers of `QuadGK.jl`, `HCubature.jl`, and `Cuba.jl` for
developing the Julia numerical integration ecosystem in which this package is situated.

# References
