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
A description of the high-level functionality and purpose of the software for a diverse, non-specialist audience.
-->

`HAdaptiveIntegration.jl` is a `Julia` [@Julia] package for automatic adaptive numerical integration on multidimensional simplices and axis-aligned orthotopes.
It approximates integrals of the form
$$
  I = \int_{\Omega} f(\boldsymbol{x}) \, \operatorname{d}\!\boldsymbol{x}
$$
where:

- $f \colon \mathbb{R}^d \to \mathbb{T}$ is any `Julia` function mapping $d$-dimensional vectors to elements of type $\mathbb{T}$ supporting multiplication by a real scalar, addition and a norm (*i.e.* a normed real vector space).
  Classic example are $\mathbb{T} = \mathbb{R},\ \mathbb{C},\ \mathbb{R}^n,\ \mathbb{R}^{n \times m}$.
- $\Omega \subset \mathbb{R}^d$ is the integration domain: simplices (triangle, tetrahedron, etc.) and axis-aligned orthotopes (rectangle, cuboid, etc.) are supported.

The package targets adaptive cubature in low to moderate dimension, with a particular emphasis on domains that arise naturally in mesh-based scientific computing.
`HAdaptiveIntegration` combines uniform subdivision rules with embedded cubature, so that it returns both an integral estimate and an a posteriori error estimate.

Its main features are:

- adaptive integration over simplices and orthotopes of arbitrary dimension;
- efficient tabulated cubature rules for low-dimensional simplices and orthotopes;
- support for user-defined embedded cubature rules and subdivision strategies;
- arbitrary-precision arithmetic.

## Usage

The package is centered around a single function, `integrate(f, domain)`, together with constructors for common domains.

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

The return value is a pair `(I, E)`, where `I` is the integral estimate and `E` is an a posteriori error estimate.

# Statement of need

<!--
A section that clearly illustrates the research purpose of the software and places it in the context of related work.
This should clearly state what problems the software is designed to solve, who the target audience is, and its relation to other work.
-->

Adaptive numerical integration is a routine building block in scientific computing, including finite-element and boundary-element methods, validation of semi-analytical formula, and parameter studies where repeated integrals must be computed over reference or physical cells.
In `Julia`, existing packages already cover important parts of this space: `QuadGK.jl` is the natural tool for one-dimensional adaptive integration [@QuadGK], `HCubature.jl` provides adaptive cubature on axis-aligned orthotopes [@HCubature], and `Cuba.jl` exposes a family of multidimensional integration methods, including stochastic ones, for rectangular domains [@Cuba].

The main contribution is that, to our knowledge, there were no `Julia` packages capable of doing adaptive integration on simplices (triangles, tetrahedra, etc.).
The second contribution is the use of specialized tabulated rules, which allow having higher order scheme than `HCubature.jl`.
We design the package in such way that there is a simple common interface for integrating over simplices and orthotopes.
And, it is easy to change the rule.
Therefore, extending the package with new rule is straightforward.
The package targets users who need deterministic adaptive integration with explicit error estimates, access to custom cubature and subdivision strategies, and the option to move beyond `Float64` arithmetic when validation or sensitive applications require it.

# State of the field

<!--
A description of how this software compares to other commonly-used packages in the research area.
If related tools exist, provide a clear “build vs. contribute” justification explaining your unique scholarly contribution and why existing alternatives are insufficient.
-->

`HAdaptiveIntegration` is best viewed as complementary to the existing `Julia` ecosystem rather than as a wholesale replacement.
`QuadGK.jl` remains the appropriate choice in one dimension because it is specialized for segments and offers very mature Gauss-Kronrod based integration [@QuadGK].
`HCubature.jl` is the closest comparison for multidimensional box domains and is often the better option for medium- and high-dimensional orthotopes [@HCubature].
`Cuba.jl` remains attractive for high-dimensional rectangular problems where stochastic or importance-sampling methods are more suitable than deterministic refinement [@Cuba].

The distinct contribution of `HAdaptiveIntegration` is its support for simplices of arbitrary dimension together with axis-aligned orthotopes under one API, and its use of efficient tabulated rules in low-dimensional cases.
That combination is not provided by the packages above.
A standalone package is therefore justified: domain geometry influences the reference mapping, admissible default cubature rules, subdivision schemes, and the practical user workflow.
Exposing these choices in a single package makes it easier to write research code that works uniformly on triangles, tetrahedra, rectangles, cuboids, and their higher-dimensional counterparts.

# Software design

<!--
An explanation of the trade-offs you weighed, the design/architecture you chose, and why it matters for your research application.
This should demonstrate meaningful design thinking beyond a superficial code structure description.
-->

The `HAdaptiveIntegration` package is organized with extending it in mind.
It has two submodules: the `Domain` module which define the possible integration domain and how to subdivide them; and the `Rule` module which define the integration rules.
There is different type of rules, the explicit rules which are computed using explicit formula, and tabulated rules which are stored in a decimal format.
The package has basically one entry point the `integrate` function.
The automatic adaptive algorithm it uses has two component: an embedded cubature and a subdivision strategy described in the following sections [Embedded cubature] and [The adaptive algorithm].

## Embedded cubature

At the core of adaptive numerical integration method, there is a cubature pair $(\mathcal{H}, \mathcal{L})$, define on the reference domain $\widehat{\Omega}$ ($\{\boldsymbol{x} \in \mathbb{R}_+^d \mid x_1 + \cdots + x_d \leq 1\}$ for simplices or $[0, 1]^d$ for orthotopes) by
$$
  \mathcal{H}(f) = \sum_{1 \leq i \leq \mathsf{H}} h_i \, f(\boldsymbol{x}_i)
  \quad \text{and} \quad
  \mathcal{L}(f) = \sum_{1 \leq i \leq \mathsf{L}} \ell_i \, f(\boldsymbol{x}_i),
  \qquad \forall f \in \mathscr{C}^0(\widehat{\Omega}).
$$
Where $\boldsymbol{x}_1, \ldots, \boldsymbol{x}_{\mathsf{H}} \in \widehat{\Omega}$ are the cubature point on the reference domain, $h_1, \ldots, h_{\mathsf{H}} \in \mathbb{R}$ are the $\mathcal{H}$ cubature weights, $\ell_1, \ldots, \ell_{\mathsf{L}} \in \mathbb{R}$ are the $\mathcal{L}$ cubature weights, and $\mathsf{H} > \mathsf{L}$ so the $\mathcal{H}$ rule has more points than the $\mathcal{L}$ rule.
The pair $(\mathcal{H}, \mathcal{L})$ is called an embedded cubature because the $\mathcal{L}$ rule use a subset of the points of the $\mathcal{H}$ rule.
In accordance with the number of point evaluation of the two cubature rule, $\mathcal{H}$ has order $d_{\mathcal{H}}$ and $\mathcal{L}$ has order $d_{\mathcal{L}}$ with $d_{\mathcal{H}} \geq d_{\mathcal{L}}$.
The order of a cubature rule is the highest integer $k \in \mathbb{N}$ such that the cubature is exact on the space of polynomials with total degree less or equal than $k$.

For a domain $\Omega$, we define the map $\phi \colon \widehat{\Omega} \to \Omega$ from the reference domain to the physical domain.
Using the embedded cubature, we define the estimated integral value $I_\Omega$ and the estimated error $E_\Omega$ by 
$$
  I_\Omega = \lvert\det \operatorname{J}_\phi\rvert \ \mathcal{H}(f \circ \phi)
  \quad \text{and} \quad
  E_\Omega = \lvert\det \operatorname{J}_\phi\rvert \ \lVert \mathcal{H}(f \circ \phi) - \mathcal{L}(f \circ \phi) \rVert,
$$
where $\operatorname{J}_\phi$ is the Jacobian of $\phi$ which is constant for simplices and axis-aligned orthotopes.

Each domain type has a default embedded cubature resume in Table 1.

|                                           |                                   |
| :---------------------------------------- | :-------------------------------- |
| 1d. `Segment` [@Laurie1997]               |                                   |
| 2d. `Triangle` [@Laurie1982]              | `Rectangle` [@CoolsHaegemans1989] |
| 3d. `Tetrahedron` [@GrundmannMoeller1978] | `Cuboid` [@BerntsenEspelid1988]   |
| $n$d. `Simplex` [@GrundmannMoeller1978]   | `Orthotope` [@GenzMalik1980]      |

Table 1: For each dimension and domain give the default rules used.

## The adaptive algorithm

Now, we have a way of computing an estimated integral value and error on a domain.
Given a function $f$ and an initial domain $\Omega$, the adaptive algorithm constructs a sequence of nested partitions.
It starts with $\mathcal{P}_0 = \{(\Omega, I_\Omega, E_\Omega)\}$ then
$$
  \mathcal{P}_{n+1} = \left[\mathcal{P}_n \setminus \left\{(\omega_0, I_{\omega_0}, E_{\omega_0})\right \}\right] \cup \left\{(\omega_1, I_{\omega_1}, E_{\omega_1}), \ldots, (\omega_{2^d}, I_{\omega_{2^d}}, E_{\omega_{2^d}})\right \},
  \qquad \forall n \in \mathbb{N},
$$
where $\omega_0$ is chosen such that $E_{\omega_0} = \max_{(\omega, I_\omega, E_\omega) \in \mathcal{P}_n} E_\omega$, and $\omega_1, \ldots, \omega_{2^d}$ are subdomain from the subdivision of $\omega_0$.
The subdivision follows the geometry of the domain.
In dimension $d$, orthotopes are bisected along each axis, producing $2^d$ subdomains, while simplices are subdivided by midpoint edge simplex refinement, again producing $2^d$ subdomain, see [@SimplexSubdiv].

Associated to the ${(\mathcal{P}_n)}_{n \in \mathbb{N}}$ sequence, we define the global integral value $I_n$ and error $E_n$ estimators by
$$
  I_n = \sum_{(\omega, I_\omega, E_\omega) \in \mathcal{P}_n} I_\omega
  \quad \text{and} \quad
  E_n = \sum_{(\omega, I_\omega, E_\omega) \in \mathcal{P}_n} E_\omega.
$$
For this type of algorithm, the stopping condition is control by three parameters: the absolute tolerance $\mathtt{atol} \geq 0$, the relative tolerance $\mathtt{rtol} \geq 0$, and the number of subdivision maximum $n_{\max} \in \mathbb{N}$.
The subdivision process stops when
$$
  E_n \leq \mathtt{atol}
  \quad \text{or} \quad
  E_n \leq \mathtt{rtol}\, \lVert I_n \rVert
  \quad \text{or} \quad
  n = n_{\max}.
$$
At the end $(I_n, E_n)$ is return as the integral value and error estimate.

## Implementation

The implementation follows closely the algorithm describe in [The adaptive algorithm].
One important implementation detail is that it uses a mutable binary heap to store the partition $\mathcal{P}_n$.
As this data structure allows efficient retrieval of a maximum local error, as well as efficient `pop` and `push` operation.
In addition, in the case of doing multiple integral on the same domain type and function signature, we can pre-allocate the heap and pass it to the `integrate` function via the `buffer` keyword to reduce the number of allocation.

## Extended precision

As said in the [Summary] section, `integrate` support arbitrary precision.
However, only the rules from [@GrundmannMoeller1978; @GenzMalik1980] are generated using explicit formula, all the other are tabulated with a finite decimal precision which is incompatible with arbitrary precision.
All the tabulated rule are stored in quadruple precision.
But, if we want to use tabulated rule with higher precision, it can be done by first increasing the precision of the rule.
`HAdaptiveIntegration` addresses this with the optional `IncreasePrecisionExt` extension, which reconstructs higher-precision tabulated rules by solving the polynomial exactness conditions with Newton iterations and automatic differentiation through `ForwardDiff` [@ForwardDiff2016].

This design keeps the default package lightweight while enabling high-precision workflows
when needed. In the test suite, deliberately reduced-precision tabulated rules fail to reach
relative accuracies around $10^{-14}$ after conversion to `Float64`, whereas their
reconstructed higher-precision counterparts recover the requested accuracy. The same
mechanism can then be used together with `BigFloat` domains and functions for validation
studies or sensitive numerical experiments.

# Research impact statement

<!--
Evidence of realized impact (publications, external use, integrations) or credible near-term
significance (benchmarks, reproducible materials, community-readiness signals). The evidence
should be compelling and specific, not aspirational.
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

# AI usage disclosure

Generative AI was used to help draft and edit parts of this manuscript. All text produced
with AI assistance was reviewed, revised, and verified by the authors against the source
code, tests, and benchmark outputs before inclusion in the paper.

# Acknowledgements

The authors thank the maintainers of `QuadGK.jl`, `HCubature.jl`, and `Cuba.jl` for
developing the `Julia` numerical integration ecosystem in which this package is situated.

# References
