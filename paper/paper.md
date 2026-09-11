---
title: '`HAdaptiveIntegration.jl`: Adaptive numerical integration over simplices and orthotopes'
tags:
  - Julia
  - numerical integration
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
date: 29 May 2026
bibliography: paper.bib
---

# Summary

`HAdaptiveIntegration.jl` is a `Julia` [@Julia] package that computes the numerical value of an integral over a geometric region in any number of dimensions, automatically refining where the integrand is hard to resolve, and returning both the value and an a posteriori error estimate.
Such integrals arise in finite-element and boundary-element methods and in parameter studies over integrands with localized features.
More precisely, it performs adaptive numerical integration on simplices (triangles, tetrahedra, and their higher-dimensional analogs) and axis-aligned orthotopes (rectangles, cuboids, and their higher-dimensional analogs), approximating integrals of the form
$$
  I = \int_{\Omega} f(\boldsymbol{x}) \, \operatorname{d}\!\boldsymbol{x}
$$
where $f \colon \mathbb{R}^d \to \mathbb{T}$ takes values in a normed real vector space (*e.g.* $\mathbb{T} = \mathbb{R},\ \mathbb{C},\ \mathbb{R}^n$), and $\Omega \subset \mathbb{R}^d$ is a simplex or an axis-aligned orthotope.

Its main features are:

- Adaptive integration over simplices and orthotopes of arbitrary dimension;
- Efficient tabulated cubature rules for low-dimensional domains;
- Support for user-defined cubature rules and subdivision strategies;
- Arbitrary-precision arithmetic (formula-based rules natively; tabulated rules via an optional extension).

## Usage

The package is centered on a single function, `integrate(f, domain; kwargs...)`, together with constructors for common domains.
The workflow has two steps:

**Create a domain** using `Simplex` (defined by its vertices) or `Orthotope` (defined by its lower and upper corners), with aliases `Segment`, `Triangle`, `Tetrahedron`, `Rectangle`, `Cuboid`.

```julia
using HAdaptiveIntegration
triangle  = Triangle((0,0), (1,0), (0,1))
rectangle = Rectangle((0,0), (1,1))
```

**Integrate** using the same interface for both domain types:

```julia
f(x) = cis(sum(x)) / (sum(abs2, x) + 1e-2)
I, E = integrate(f, triangle)
I, E = integrate(f, rectangle)
```

The function returns a pair `(I, E)`, where `I` is the integral estimate and `E` is an a posteriori error estimate. The main stopping-condition keywords are `atol`, `rtol`, and `maxsubdiv`.

# Statement of need

Adaptive numerical integration is a fundamental building block in scientific computing, including finite-element and boundary-element methods and parameter studies over reference or physical cells.
`HAdaptiveIntegration` provides this for simplices and orthotopes, with a single unified interface, user-extensible cubature, and optional arbitrary-precision arithmetic.

# State of the field

`HAdaptiveIntegration` is complementary to the existing `Julia` ecosystem.
`QuadGK.jl` [@QuadGK] remains the right choice in one dimension.
`HCubature.jl` [@HCubature] is often preferable for medium- to high-dimensional orthotopes, where its estimator-driven axis-aligned bisection can outperform isotropic subdivision if the integrand varies primarily along one direction.
`Cuba.jl` [@Cuba] is better suited for high-dimensional problems where stochastic methods dominate.
The distinct contribution of `HAdaptiveIntegration` is support for simplices of arbitrary dimension and orthotopes under one API with efficient tabulated rules in low dimensions, a combination not provided by the packages above.

# Software design

The package has two submodules: `Domain` (domains and subdivision) and `Rule` (cubature rules, formula-based or tabulated), with `integrate` as the single entry point.

## Method

### Embedded cubature

At the core of an adaptive numerical integration method is a cubature pair $(\mathcal{H}, \mathcal{L})$, defined on the reference domain $\widehat{\omega}$ (the unit simplex or $[0,1]^d$) by
$$
  \mathcal{H}(f) = \sum_{1 \leq i \leq \mathsf{H}} h_i \, f(\boldsymbol{x}_i)
  \quad \text{and} \quad
  \mathcal{L}(f) = \sum_{1 \leq i \leq \mathsf{L}} \ell_i \, f(\boldsymbol{x}_i),
$$
where $\boldsymbol{x}_i$ are the nodes, $h_i$ and $\ell_i$ the respective weights, and $\mathsf{H} > \mathsf{L}$.
The $\mathcal{L}$ rule reuses the first $\mathsf{L}$ nodes of $\mathcal{H}$ (hence *embedded*) with generally different weights $\ell_i \neq h_i$, so the pair costs only $\mathsf{H}$ evaluations rather than $\mathsf{H}+\mathsf{L}$.
The two rules have polynomial exactness orders $k_h > k_\ell$ (the highest total degree integrated exactly).

For a domain $\omega$ with reference map $\phi \colon \widehat{\omega} \to \omega$, the *local* estimated integral value $I_\omega$ and error $E_\omega$ are
$$
  I_\omega = \lvert\det \operatorname{J}_\phi\rvert \ \mathcal{H}(f \circ \phi)
  \quad \text{and} \quad
  E_\omega = \lvert\det \operatorname{J}_\phi\rvert \ \lVert \mathcal{H}(f \circ \phi) - \mathcal{L}(f \circ \phi) \rVert,
$$
where $\operatorname{J}_\phi$ is the (constant) Jacobian and $\lVert \cdot \rVert$ is the norm on $\mathbb{T}$.

Default rules per domain are summarized in \autoref{tbl:default-rule} and visualized in \autoref{fig:embedded_cubature}.

| Dimension | Domain | $k_h$ | $k_\ell$ | Reference |
| :-- | :-- | --: | --: | :-- |
| 1d | `Segment` | 23 | 13 | [@Laurie1997] |
| 2d | `Triangle` | 8 | 5 | [@Laurie1982] |
| 2d | `Rectangle` | 7 | 5 | [@GenzMalik1980] |
| 3d | `Tetrahedron` | 7 | 5 | [@GrundmannMoeller1978] |
| 3d | `Cuboid` | 9 | 7 | [@BerntsenEspelid1988] |
| $n$d | `Simplex` | 7 | 5 | [@GrundmannMoeller1978] |
| $n$d | `Orthotope` | 7 | 5 | [@GenzMalik1980] |
: Default embedded cubature rules by domain, where $k_h$ and $k_\ell$ are the polynomial exactness orders of the high and low rules. \label{tbl:default-rule}

![Plots of the default embedded cubature rules used for the triangle and rectangle. Cross markers show the $\mathcal{H}$ nodes and circle markers the $\mathsf{L}$ shared nodes reused by $\mathcal{L}$; color encodes the weight. \label{fig:embedded_cubature}](image/embedded_cubature.pdf){width=\textwidth}

### The adaptive algorithm

Given a function $f$ and an initial domain $\Omega$, the adaptive algorithm constructs a sequence of nested partitions.
It starts with $\mathcal{P}_0 = \{\Omega\}$ and iterates by
$$
  \mathcal{P}_{n+1} = \left[ \mathcal{P}_n \setminus \left\{\omega^*\right\} \right] \cup \left\{\omega^*_1, \ldots, \omega^*_{2^d}\right \},
  \qquad \forall n \in \mathbb{N},
$$
where $\omega^*$ is chosen such that $E_{\omega^*} = \max \{E_\omega : \omega \in \mathcal{P}_n\}$, and $\omega^*_1, \ldots, \omega^*_{2^d}$ are subdomains given by a partition of $\omega^*$.
In dimension $d$, orthotopes are bisected along each axis and simplices by midpoint edge refinement [@SimplexSubdiv], both producing $2^d$ subdomains.

For the sequence ${(\mathcal{P}_n)}_{n \in \mathbb{N}}$, we define the global integral and error estimators $\mathcal{I}_n$ and $\mathcal{E}_n$ by
$$
  \mathcal{I}_n = \sum_{\omega \in \mathcal{P}_n} I_\omega
  \quad \text{and} \quad
  \mathcal{E}_n = \sum_{\omega \in \mathcal{P}_n} E_\omega.
$$
The process stops when $\mathcal{E}_n \leq \mathtt{atol}$ (absolute tolerance) or $\mathcal{E}_n \leq \mathtt{rtol}\,\lVert \mathcal{I}_n \rVert$ (relative tolerance), or $n = \mathtt{maxsubdiv}$ (maximum number of subdivisions), and returns $(\mathcal{I}_n, \mathcal{E}_n)$.

![One level of isotropic subdivision for a triangle and tetrahedron split by joining edge midpoints, and a rectangle and cuboid bisected along each axis. \label{fig:subdivision}](image/subdivision.pdf){width=\textwidth}

## Implementation

The only domain structs are `Simplex{D,T,N}` (with $N = D+1$ vertices) and `Orthotope{D,T}` (low and high corners); `Segment`, `Triangle`, `Tetrahedron`, `Rectangle`, and `Cuboid` are type aliases, so dimension-specific defaults are selected by ordinary dispatch rather than run-time branching.
Vertices and cubature nodes are `SVector`s, so the dimension $D$ is part of the type and the inner evaluation loop is stack-allocated and specialized per $(D, T)$.
Default rules are built by `@generated` functions, so the `EmbeddedCubature` for a given element type is constructed once at compile time.
An `EmbeddedCubature` stores the shared nodes and separate high- and low-order weight vectors, the first $\mathsf{L}$ nodes being those of the low-order rule.
Rules are tabulated on the reference domain and mapped by `map_from_reference`, which returns $(\phi, \lvert\det \operatorname{J}_\phi\rvert)$ with constant Jacobian.

The adaptive set lives in a max binary heap from [DataStructures.jl](https://github.com/JuliaCollections/DataStructures.jl), ordered by local error, and `integrate` maintains $(\mathcal{I}_n, \mathcal{E}_n)$ incrementally.
When computing many integrals of the same type, the heap can be pre-allocated and passed via the `buffer` keyword so that, provided the integrand does not allocate, the only allocation is the returned pair; see the documentation for a zero-allocation benchmark.
Defaults are `atol = 0`, `rtol = sqrt(eps(T))` for element type `T`, and $\mathtt{maxsubdiv} = 2^{13+D}$; if `maxsubdiv` is reached, a warning is issued and the current estimate is returned.
The same code path runs with `Float64`, `BigFloat`, and [`Unitful.jl`](https://github.com/JuliaPhysics/Unitful.jl) quantities.
Extensibility is provided through the `rule` and `norm` keywords on `integrate`, which let users plug in their own cubature and error norm.

## Extended precision

As noted above, `integrate` supports arbitrary precision.
Only the rules from @GrundmannMoeller1978 and @GenzMalik1980 are formula-based and already arbitrary-precision; the others are tabulated at quadruple precision and therefore incompatible with it.
`HAdaptiveIntegration` addresses this with the optional `IncreasePrecisionExt` extension, which refines a tabulated rule's nodes and weights by a least-squares Newton iteration on the polynomial-exactness (moment) conditions for both $\mathcal{H}$ and $\mathcal{L}$, with Jacobians from `ForwardDiff.jl` [@ForwardDiff2016], following [@XiaoGimbutas2010, section 2.3].
The full construction, including the definition of the refinement operator and its stopping criteria, is given in the documentation.

**Remark.** The monomial basis is convenient but poorly conditioned, so achieving precision $\varepsilon$ requires `BigFloat` arithmetic at a higher internal precision $\eta < \varepsilon$; an $\mathrm{L}^2$-orthogonal basis would improve conditioning.

**Remark.** Symmetry is not enforced in the Newton iteration; a rule refined from precision $\varepsilon$ to $\eta < \varepsilon$ therefore retains its symmetry only up to $\varepsilon$.

# Research impact statement

The repository includes API documentation, examples (buffer pre-allocation, `callback` mechanism, custom cubature, arbitrary-precision workflows), an automated test suite, and continuous integration covering tests, documentation, and linting.
`HAdaptiveIntegration` has been featured in the [Julia world newsletter](https://discourse.julialang.org/t/this-month-in-julia-world-2026-02/136110), interfaced via [Integrals.jl](https://github.com/SciML/Integrals.jl) (SciML ecosystem), used as a backend in @Inti, and adopted as a dependency of [`Meshes.jl`](https://github.com/JuliaGeometry/Meshes.jl/releases/tag/v0.57.0).

# Example gallery

We showcase the package on integrands with localized features constructed from the mollified integrand
$$
  \rho_\delta(\psi(\boldsymbol{x})) = \frac{1}{\delta^c} \rho \left(\frac{\psi(\boldsymbol{x})}{\delta}\right),
  \qquad \rho(r) = \frac{1}{\sqrt{2\pi}} e^{-r^2/2},
$$
where $\psi$ is a level-set function vanishing on $\Gamma = \psi^{-1}(0)$, $c$ is the codimension of $\Gamma$, so that the total mass stays $\mathcal{O}(1)$ as $\delta \to 0$; $\rho_\delta \circ \psi$ varies rapidly near $\Gamma$.
We consider three canonical geometries:

- **Point:** $\psi(\boldsymbol{x}) = \lVert\boldsymbol{x} - \boldsymbol{z}_0\rVert_2$, codimension $c = d$ (point refinement);
- **Hypersphere:** $\psi(\boldsymbol{x}) = \lVert\boldsymbol{x} - \boldsymbol{y}\rVert_2^2 - r^2$, codimension $c = 1$ (curve/surface refinement);
- **Hyperplane:** $\psi(\boldsymbol{x}) = x_1 - a$, codimension $c = 1$ (axis-aligned hyperplane refinement).

We set $\delta = 0.05$, $\boldsymbol{z}_0 = \tfrac{1}{\pi}(1, \ldots, 1)$, $\boldsymbol{z}_1 = (1, 0, \ldots, 0)$, $r = 2/3$, $a = 1/\pi$; the axis-aligned hyperplane benchmarks isotropic versus estimator-driven subdivision.
Convergence plots sweep $\mathtt{rtol} = 10^{-i}$ ($i = 1, \ldots, 10$; up to $8$ in 3D), recording the total number of evaluations $N$, the returned error estimate, and the actual error against a reference solution computed at $\mathtt{rtol} = 10^{-12}$.

When one evaluation of $f$ dominates the per-node overhead of the algorithm (a few clock cycles), the wall-clock time is roughly $N$ times the cost of a single evaluation, so $N$ is a hardware-independent proxy for runtime.

## Simplices

\autoref{fig:cvg_simplex} shows convergence on the unit triangle ($d=2$, Radon-Laurie rule [@Laurie1982]) and tetrahedron ($d=3$, Grundmann-Möller rule [@GrundmannMoeller1978], available in arbitrary dimensions).
Two observations hold across all features and both geometries: the estimated error reliably tracks the actual error, confirming a sound a posteriori indicator; and once the feature is resolved, errors follow $\mathcal{O}(N^{-(k+1)/d})$ with $k = k_h$ or $k_\ell$ and the exponent implied by $N \propto h^{-d}$.

![](image/cvg_triangle.pdf){width=\textwidth}

![Convergence of the actual and estimated errors for the point, hypersphere, and hyperplane features on the unit triangle (top) and unit tetrahedron (bottom) versus the number of evaluations ($N$). Insets show the sub-domains and the integrand as max-intensity projections onto the coordinate planes. \label{fig:cvg_simplex}](image/cvg_tetrahedron.pdf){width=\textwidth}

## Orthotopes and comparison with `HCubature.jl`

\autoref{fig:cvg_orthotope} compares `HAdaptiveIntegration` against `HCubature.jl` on the unit square and cube.
In 2D, both use the Genz-Malik rule [@GenzMalik1980], isolating the subdivision strategy; in 3D, the cubature rules also differ (Berntsen-Espelid [@BerntsenEspelid1988] for `HAdaptiveIntegration`).
For the point and hypersphere features the solvers are comparable in both dimensions.
For the hyperplane, `HCubature.jl` has a clear advantage because its estimator refines exclusively along $x_1$ rather than bisecting uniformly in all $d$ directions, and this gap grows with $d$, motivating anisotropic splitting as future work.

![](image/cvg_rectangle.pdf){width=\textwidth}

![Convergence of the actual and estimated errors for the point, hypersphere, and hyperplane features on the unit square (top) and the unit cube (bottom), comparing `HAdaptiveIntegration` (HAI) with `HCubature.jl`. Insets show the sub-domains and the integrand as max-intensity projections onto the coordinate planes. \label{fig:cvg_orthotope}](image/cvg_cuboid.pdf){width=\textwidth}

# AI usage disclosure

Generative AI was used to help develop the code and to draft and edit parts of this manuscript.
All text produced with AI assistance was reviewed, revised, and verified by the authors before inclusion in the code or paper.

# References
