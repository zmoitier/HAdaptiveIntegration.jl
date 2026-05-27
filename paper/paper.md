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

The package is centered around a single function, `integrate(f, domain; kwargs...)`, together with constructors for common domains. The usual workflow is comprised of the following two steps:

**1. Create a domain.**
There are two main domain constructors:

- `Simplex`: defined by its vertices;
- `Orthotope`: defined by its low and high corner.

For convenience, the following aliases are provided for the most common cases: `Segment`,
`Triangle`, `Tetrahedron`, `Rectangle`, `Cuboid`.

For example:

```julia
using HAdaptiveIntegration
segment     = Segment(0, 1)
triangle    = Triangle((0,0), (1,0), (0,1))
rectangle   = Rectangle((0,0), (1,1))
simplex4    = Simplex((0,0,0,0), (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1))
```

**Integrate.**
Once the domain is defined, the same interface is used for simplices and orthotopes:

```julia
f(x) = cis(sum(x)) / (sum(abs2, x) + 1e-2)
I, E = integrate(f, triangle)
I, E = integrate(f, rectangle)
```

The return value is a pair `(I, E)`, where `I` is the integral estimate and `E` is an a posteriori error estimate.

Additionally, the `integrate` function accepts keyword arguments to control the stopping
condition of the adaptive algorithm and the choice of embedded cubature. The most important ones, controlling the stopping condition, are:

- `atol`: absolute tolerance for the error estimate;
- `rtol`: relative tolerance for the error estimate;
- `maxsubdiv`: maximum number of subdivisions.

# Statement of need

<!--
A section that clearly illustrates the research purpose of the software and places it in the context of related work.
This should clearly state what problems the software is designed to solve, who the target audience is, and its relation to other work.
-->

Adaptive numerical integration is a routine building block in scientific computing, including finite-element and boundary-element methods, validation of semi-analytical formula, and parameter studies where repeated integrals must be computed over reference or physical cells.
In `Julia`, existing packages already cover important parts of this space: `QuadGK.jl` [@QuadGK] is the package to use for one-dimensional adaptive integration, `HCubature.jl` [@HCubature] provides adaptive cubature on axis-aligned orthotopes, and `Cuba.jl` [@Cuba] exposes a family of multidimensional integration methods, including stochastic ones, for rectangular domains.

The main contribution is that, to our knowledge, there were no `Julia` packages capable of doing adaptive integration on simplices (triangles, tetrahedra, etc.).
The second contribution is the use of specialized tabulated rules, which allow higher-order schemes than `HCubature.jl`.
The package exposes a single common interface for both simplices and orthotopes, and adding a new cubature rule requires only implementing a small interface, making it straightforward to extend.
The package targets users who need deterministic adaptive integration with explicit error estimates, access to custom cubature, and the option to move beyond `Float64` arithmetic when validation or sensitive applications require it.

# State of the field

<!--
A description of how this software compares to other commonly-used packages in the research area.
If related tools exist, provide a clear “build vs. contribute” justification explaining your unique scholarly contribution and why existing alternatives are insufficient.
-->

`HAdaptiveIntegration` is best viewed as complementary to the existing `Julia` ecosystem rather than as a wholesale replacement.
`QuadGK.jl` [@QuadGK] remains the appropriate choice in one dimension because it is specialized for segments and offers very mature Gauss-Kronrod integration.
`HCubature.jl` [@HCubature] is the closest comparison for multidimensional box domains and is often the better option for medium- and high-dimensional orthotopes due to a different refinement strategy (uniform for `HAdaptiveIntegration` versus estimator based for `HCubature`).
`Cuba.jl` [@Cuba] remains attractive for high-dimensional rectangular problems where stochastic or importance-sampling methods are more suitable than deterministic refinement.

The distinct contribution of `HAdaptiveIntegration` is its support for simplices of arbitrary dimension together with axis-aligned orthotopes under one API, and its use of efficient tabulated rules in low-dimensional cases.
That combination is not provided by the packages above.
A standalone package is therefore justified: domain geometry influences the reference mapping, admissible default cubature rules, subdivision schemes, and the practical user workflow.
Exposing these choices in a single package makes it easier to write research code that works uniformly on triangles, tetrahedra, rectangles, cuboids, and their higher-dimensional counterparts.

# Software design

<!--
An explanation of the trade-offs you weighed, the design/architecture you chose, and why it matters for your research application.
This should demonstrate meaningful design thinking beyond a superficial code structure description.
-->

The `HAdaptiveIntegration` package is organized so as to make it easy to use and extend.
It has two submodules: the `Domain` module which defines the possible integration domains and how to subdivide them; and the `Rule` module which defines the integration rules.
There are different types of rules: explicit rules which are computed using explicit formulas, and tabulated rules which are stored in a decimal format.
The package has essentially one entry point: the `integrate` function.
The automatic adaptive algorithm it uses has two components: an embedded cubature and a subdivision strategy described in the following sections [Embedded cubature] and [The adaptive algorithm].

## Embedded cubature

At the core of adaptive numerical integration method, there is a cubature pair $(\mathcal{H}, \mathcal{L})$, define on the reference domain $\widehat{\omega}$ (equal to $\{\boldsymbol{x} \in \mathbb{R}_+^d \mid x_1 + \cdots + x_d \leq 1\}$ for simplices or $[0, 1]^d$ for orthotopes) by
$$
  \mathcal{H}(f) = \sum_{1 \leq i \leq \mathsf{H}} h_i \, f(\boldsymbol{x}_i)
  \quad \text{and} \quad
  \mathcal{L}(f) = \sum_{1 \leq i \leq \mathsf{L}} \ell_i \, f(\boldsymbol{x}_i),
  \qquad \forall f \in \mathscr{C}^0(\widehat{\omega}).
$$
where $\boldsymbol{x}_1, \ldots, \boldsymbol{x}_{\mathsf{H}} \in \widehat{\omega}$ are the cubature points on the reference domain, $h_1, \ldots, h_{\mathsf{H}} \in \mathbb{R}$ are the $\mathcal{H}$ cubature weights, $\ell_1, \ldots, \ell_{\mathsf{L}} \in \mathbb{R}$ are the $\mathcal{L}$ cubature weights, and $\mathsf{H} > \mathsf{L}$ so the $\mathcal{H}$ rule has more points than the $\mathcal{L}$ rule.
The pair $(\mathcal{H}, \mathcal{L})$ is called an embedded cubature because the $\mathcal{L}$ rule uses a subset of the points of the $\mathcal{H}$ rule.
In accordance with the number of point evaluations of the two cubature rules, $\mathcal{H}$ has order $k_{\mathcal{H}}$ and $\mathcal{L}$ has order $k_{\mathcal{L}}$ with $k_{\mathcal{H}} > k_{\mathcal{L}}$.
The order of a cubature rule is the highest integer $k \in \mathbb{N}$ such that the cubature is exact on the space of polynomials with total degree less or equal than $k$.

For a domain $\omega$, we define the map $\phi \colon \widehat{\omega} \to \omega$ from the reference domain to the physical domain.
Using the embedded cubature, we define the *local* estimated integral value $I_\omega$ and the *local* estimated error $E_\omega$ by
$$
  I_\omega = \lvert\det \operatorname{J}_\phi\rvert \ \mathcal{H}(f \circ \phi)
  \quad \text{and} \quad
  E_\omega = \lvert\det \operatorname{J}_\phi\rvert \ \lVert \mathcal{H}(f \circ \phi) - \mathcal{L}(f \circ \phi) \rVert,
$$
where $\operatorname{J}_\phi$ is the Jacobian of $\phi$ (which is constant for simplices and axis-aligned orthotopes) and $\lVert \cdot \rVert$ is the chosen norm on $\mathbb{T}$.

Each domain type has a default embedded cubature summarized in \autoref{tbl:default-rule}.

|                                           |                                 |
| :---------------------------------------- | :------------------------------ |
| 1d. `Segment` [@Laurie1997]               |                                 |
| 2d. `Triangle` [@Laurie1982]              | `Rectangle` [@GenzMalik1980]    |
| 3d. `Tetrahedron` [@GrundmannMoeller1978] | `Cuboid` [@BerntsenEspelid1988] |
| $n$d. `Simplex` [@GrundmannMoeller1978]   | `Orthotope` [@GenzMalik1980]    |
: For each dimension and domain give the default rules used.\label{tbl:default-rule}

## The adaptive algorithm

Now, we have a way of computing an estimated integral value and error on a domain.
Given a function $f$ and an initial domain $\Omega$, the adaptive algorithm constructs a sequence of nested partitions.
It starts with $\mathcal{P}_0 = \{\Omega\}$ then
$$
  \mathcal{P}_{n+1} = \left[ \mathcal{P}_n \setminus \left\{\omega_0\right\} \right] \cup \left\{\omega_1, \ldots, \omega_{2^d}\right \},
  \qquad \forall n \in \mathbb{N},
$$
where $\omega_0$ is chosen such that $E_{\omega_0} = \max \{E_\omega : \omega \in \mathcal{P}_n\}$, and $\omega_1, \ldots, \omega_{2^d}$ are subdomains given by a subdivision of $\omega_0$.
The subdivision follows the geometry of the domain.
In dimension $d$, orthotopes are bisected along each axis, producing $2^d$ subdomains, while simplices are subdivided by midpoint edge simplex refinement, again producing $2^d$ subdomains, see [@SimplexSubdiv].

Associated to the ${(\mathcal{P}_n)}_{n \in \mathbb{N}}$ sequence, we define the global integral value $I_n$ and error $E_n$ estimators by
$$
  I_n = \sum_{\omega \in \mathcal{P}_n} I_\omega
  \quad \text{and} \quad
  E_n = \sum_{\omega \in \mathcal{P}_n} E_\omega.
$$
For this type of algorithm, the stopping condition is controlled by three parameters: the absolute tolerance $\mathtt{atol} \geq 0$, the relative tolerance $\mathtt{rtol} \geq 0$, and the number of maximum subdivisions $n_{\max} \in \mathbb{N}$.
The subdivision process stops when
$$
  E_n \leq \mathtt{atol}
  \quad \text{or} \quad
  E_n \leq \mathtt{rtol}\, \lVert I_n \rVert
  \quad \text{or} \quad
  n = n_{\max}.
$$
At the end $(I_n, E_n)$ is returned as the integral value and error estimate.

## Implementation

The implementation follows closely the algorithm described in section [The adaptive algorithm].
One important implementation detail is the use of a binary heap to store the partition $\mathcal{P}_n$.
More precisely, a max binary heap storing $\{(\omega, I_\omega, E_\omega) : \omega \in \mathcal{P}_n \}$ is used, and the comparison only involves $E_\omega$.
This data structure allows efficient retrieval of a maximum local error, as well as efficient `pop` and `push` operations.
In addition, in the case of computing multiple integrals on the same domain type and function signature, the heap can be pre-allocated and passed to the `integrate` function via the `buffer` keyword to reduce allocations.

## Extended precision

As noted in the [Summary] section, `integrate` supports arbitrary precision.
However, only the rules from [@GrundmannMoeller1978; @GenzMalik1980] are generated using explicit formulas; all others are tabulated with finite decimal precision (quadruple precision), which is incompatible with arbitrary precision.
`HAdaptiveIntegration` addresses this with the optional `IncreasePrecisionExt` extension, which computes higher-precision rules from a tabulated rule by solving the polynomial exactness conditions with Newton iterations and automatic differentiation through `ForwardDiff` [@ForwardDiff2016].

To be more precise, let $\widehat{\omega}$ be a reference domain and $(\mathcal{H}, \mathcal{L})$ be an embedded cubature on $\widehat{\omega}$ with orders $k_{\mathcal{H}} > k_{\mathcal{L}}$.
Let's call $\mathbb{P}_k$ the space of polynomials with total degree less or equal than $k$, and denote $b_1, \ldots, b_{K_{\mathcal{H}}}$ a basis of $\mathbb{P}_{k_{\mathcal{H}}}$ such that $b_1, \ldots, b_{K_{\mathcal{L}}}$ is a basis of $\mathbb{P}_{k_{\mathcal{L}}}$.
By definition, we have $K_{\mathcal{H}} = \dim\mathbb{P}_{k_{\mathcal{H}}}$ and $K_{\mathcal{L}} = \dim\mathbb{P}_{k_{\mathcal{L}}}$.
Let's define $\boldsymbol{u}^{\mathcal{H}, \mathcal{L}} = (\boldsymbol{x}_1, \ldots, \boldsymbol{x}_{\mathsf{H}}, h_1, \ldots, h_{\mathsf{H}}, \ell_1, \ldots, \ell_{\mathsf{L}})$ the embedded cubature data.
We define the function $F \colon \mathbb{R}^{(d+1) \mathsf{H} + \mathsf{L}} \to \mathbb{R}^{K_{\mathcal{H}} + K_{\mathcal{L}}}$ by
$$
  F\left( \boldsymbol{u}^{\mathcal{H}, \mathcal{L}} \right) =
  \begin{pmatrix}
    \mathcal{H}(b_1) - \int_{\widehat{\Omega}} b_1(\boldsymbol{x}) \operatorname{d}\!\boldsymbol{x}
    \\[1ex]
    \vdots
    \\[1ex]
    \mathcal{H}(b_{K_{\mathcal{H}}}) - \int_{\widehat{\Omega}} b_{K_{\mathcal{H}}}(\boldsymbol{x}) \operatorname{d}\!\boldsymbol{x}
    \\[2ex]
    \mathcal{L}(b_1) - \int_{\widehat{\Omega}} b_1(\boldsymbol{x}) \operatorname{d}\!\boldsymbol{x}
    \\[1ex]
    \vdots
    \\[1ex]
    \mathcal{L}(b_{K_{\mathcal{L}}}) - \int_{\widehat{\Omega}} b_{K_{\mathcal{L}}}(\boldsymbol{x}) \operatorname{d}\!\boldsymbol{x}
  \end{pmatrix}.
$$
The goal is to find a root of the function $F$ with higher precision than the stored precision. For a tabulated rule $\boldsymbol{u}^{\mathcal{H}, \mathcal{L}}$, we already have $\lVert F(\boldsymbol{u}^{\mathcal{H}, \mathcal{L}}) \rVert_2 = \varepsilon \ll 1$, and we want to find a $\tilde{\boldsymbol{u}}$ such that $\lVert F(\tilde{\boldsymbol{u}}) \rVert_2 = \eta < \varepsilon$. To do that we use a least-square Newton method, [@XiaoGimbutas2010{section 2.3}].
In detail, set $\boldsymbol{u}_0 = \boldsymbol{u}^{\mathcal{H}, \mathcal{L}}$ and define the iteration
$$
  \boldsymbol{u}_{p+1} = \boldsymbol{u}_p - \boldsymbol{\delta}_p
  \qquad \text{where} \
  \boldsymbol{\delta}_p = \arg\min \left\{ \lVert \boldsymbol{w} \rVert_2 \mid \operatorname{J}_F(\boldsymbol{u}_p) \boldsymbol{w} = F(\boldsymbol{u}_p) \right\},\tag{N}
$$
and $\operatorname{J}_F(\boldsymbol{u}_p)$ is the Jacobian matrix of $F$ at the point $\boldsymbol{u}_p$.
The iteration stops when $\lVert \boldsymbol{u}_{p+1} - \boldsymbol{u}_p \rVert_2 \leq \mathtt{x\_atol}$ (absolute tolerance of successive iterates) or $\lVert F(\boldsymbol{u}_p) \rVert_2 \leq \mathtt{f\_atol}$ (absolute tolerance of function value) or $p = p_{\max}$ (maximum number of iterations).

**Remark.** In (N), $\operatorname{J}_F(\boldsymbol{u}_p) \boldsymbol{w} = F(\boldsymbol{u}_p)$ is solved via the least square method (just "`\`" in `Julia`).
It is written that way because, we usually have $K_{\mathcal{H}} + K_{\mathcal{L}} < (d+1) \mathsf{H} + \mathsf{L}$, so $\operatorname{J}_F(\boldsymbol{u}_p)$ is a rectangular matrix with more columns than rows.
Therefore, we choose the vector in the null space of $\operatorname{J}_F(\boldsymbol{u}_p)$ with the smallest magnitude.

**Remark.** `IncreasePrecisionExt` uses the monomial basis because of the ease of uses and the exact values of the integrals is known explicitly.
However, this basis is not well conditioned, meaning that in practice, to obtain a rule at precision $\varepsilon = 10^{-m}$ one must work in `BigFloat` with a higher internal precision $\varepsilon^\alpha$ with $\alpha > 1$. Using a better basis, such as the $L^2$-orthogonal basis, should improve the conditioning number of the problem.

**Remark.** The tabulated rules have the same symmetries as the reference domain.
However, for simplicity, we did not impose the symmetry in this method, so it does not preserve symmetry.
If a rule is tabulated with precision $\varepsilon$ and then its precision is increased to $\eta < \varepsilon$ using `IncreasePrecisionExt`,
The rule will respect the symmetry up to precision $\varepsilon$ and not $\eta$.

# Research impact statement

<!--
Evidence of realized impact (publications, external use, integrations) or credible near-term significance (benchmarks, reproducible materials, community-readiness signals).
The evidence should be compelling and specific, not aspirational.
-->

`HAdaptiveIntegration` is ready for direct use.
The repository includes API documentation, advanced examples covering the pre-allocation of the buffer, the `callback` mechanism, the custom cubature, an extension for arbitrary-precision workflows, and an automated test suite spanning all supported domains.
Continuous integration checks tests, documentation, and linting.

`HAdaptiveIntegration` was featured in the [This month in Julia world - 2026-02](https://discourse.julialang.org/t/this-month-in-julia-world-2026-02/136110) Newsletter. It has been interfaced in [@Integrals], which is part of the SciML ecosystem. It is also used as the backend for a nearly-singular integration method in [@Inti].

# Example gallery

We now showcase the package on integrands with localized features, which is where adaptive integration is most useful.
We begin with simplices, which are the most common domain in mesh-based scientific computing and the most novel feature of the package, then move to orthotopes, where we provide a brief comparison with `HCubature.jl`.

All examples are constructed from the standard Gaussian kernel

$$
  \eta(r) = \frac{1}{\sqrt{2\pi}} e^{-r^2/2},
$$

which satisfies $\int_{-\infty}^{\infty} \eta(r) \operatorname{d}\!r = 1$ and is smooth, non-negative, and peaks at $r = 0$.
Given a level-set function $\phi \colon \mathbb{R}^d \to \mathbb{R}$ that vanishes on the desired feature location $\Gamma = \phi^{-1}(0)$, we define the mollified integrand

$$
  \delta_\epsilon(\phi(\boldsymbol{x})) = \frac{1}{\varepsilon^c} \eta \left(\frac{\phi(\boldsymbol{x})}{\varepsilon}\right),
$$

where $\varepsilon > 0$ controls the width and $c$ is the **codimension** of $\Gamma$ in $\mathbb{R}^d$ (chosen so that the total mass of the integrand remains $\mathcal{O}(1)$ as $\varepsilon \to 0$).
The qualitative behavior of $\delta_\epsilon(\phi(\boldsymbol{x}))$ is that of a smooth, yet sharply peaked function, concentrated around $\Gamma$.

In the examples below we use three canonical feature geometries:

- **Point:** $\phi(\boldsymbol{x}) = \lVert\boldsymbol{x} - \boldsymbol{x}_0\rVert_2$, codimension $c = d$, producing a sharp isotropic peak at $\boldsymbol{x}_0$.
- **Hypersphere:** $\phi(\boldsymbol{x}) = \lVert\boldsymbol{x}\rVert_2^2 - r^2$, codimension $c = 1$, producing a ridge along the hypersphere of radius $r$.
- **Hyperplane:** $\phi(\boldsymbol{x}) = x_1 - c_0$, codimension $c = 1$, producing a ridge along the axis-aligned hyperplane $x_1 = c_0$.

The hyperplane case is geometrically special: because the feature is aligned with a coordinate axis, axis-aligned subdivision strategies can exploit this structure, which makes it a useful benchmark for comparing `HAdaptiveIntegration` (uniform subdivision) against `HCubature.jl` (estimator-driven directional subdivision).
We pick $\boldsymbol{x}_0 = \tfrac{1}{\pi}(1, \ldots, 1)$, $r = 1/2$, and $c_0 = 1/\pi$ so that the features are well contained in the unit simplex and unit cube.

Each convergence plot is produced by sweeping the requested relative tolerance $\mathtt{rtol} = 10^{-i}$ for $i = 1, \ldots, 10$ (and $i = 1, \ldots, 8$ for the three-dimensional cases, where each refinement step is more expensive), and recording $N$, the total number of integrand evaluations, alongside the returned integral estimate and error estimate.
The *actual* error at each tolerance is measured against a reference solution obtained by running `HAdaptiveIntegration` with the same embedded rule at $\mathtt{rtol} = 10^{-12}$.
The *estimated* error is the a posteriori quantity returned directly by the solver.

## Simplices

\autoref{fig:cvg_triangle} shows convergence results for the three features on the unit triangle using the default `RadonLaurie` embedded rule [@Laurie1982].
The inset in each panel shows the adaptive sub-triangles at default tolerance: refinement concentrates tightly around the feature, leaving the rest of the domain sparsely covered.
Two observations are consistent across all three features.
First, the estimated error returned by the solver reliably tracks the actual error throughout the refinement, confirming that the embedded rule provides a sound a posteriori indicator.
Second, once the mesh is fine enough to resolve the feature, the solver recovers the asymptotic convergence rate predicted by the cubature order: the actual error follows the high-order slope $\mathcal{O}(N^{-(k_{\mathcal{H}}+1)/2})$ and the estimated error follows the low-order slope $\mathcal{O}(N^{-(k_{\mathcal{L}}+1)/2})$, where $N$ is the number of function evaluations.
The exponent $1/2$ comes from the relation $h \sim N^{-1/d}$ with $d = 2$: in two dimensions each cell has area $h^2$, so $N \propto h^{-2}$, and a rule of order $k$ with error $\mathcal{O}(h^{k+1})$ gives $\mathcal{O}(N^{-(k+1)/2})$ overall.
All three features show qualitatively similar convergence behavior.

![Convergence of the actual and estimated errors for the point, hypersphere, and hyperplane features on the unit triangle as a function of the number of evaluations ($N$).
Insets show the adaptive sub-triangulation at default tolerance.\label{fig:cvg_triangle}](cvg_triangle.png)

The same behavior extends to higher-dimensional simplices.
\autoref{fig:cvg_tetrahedron} repeats the study on the unit tetrahedron using the `GrundmannMoeller` embedded rule [@GrundmannMoeller1978], which is available for simplices of arbitrary dimension.
The convergence rates now follow $\mathcal{O}(N^{-(k+1)/3})$ by the same argument with $d = 3$: each cell has volume $h^3$, so $N \propto h^{-3}$ and the error scales as $\mathcal{O}(N^{-(k+1)/3})$.
The estimated error again tracks the actual error reliably across all three features.

![Convergence of the actual and estimated errors for the point, hypersphere, and hyperplane features on the unit tetrahedron as a function of the number of evaluations ($N$).\label{fig:cvg_tetrahedron}](cvg_tetrahedron.png)

## Orthotopes and comparison with `HCubature.jl`

The same `integrate` interface applies to orthotopes without any change.
The convergence plots in \autoref{fig:cvg_rectangle} compare `HAdaptiveIntegration` against `HCubature.jl` on the unit square for three representative integrands.
Both solvers use the same Genz–Malik rule [@GenzMalik1980], so the comparison isolates the effect of the subdivision strategy: uniform bisection in all dimensions for `HAdaptiveIntegration` versus estimator-driven axis-aligned bisection for `HCubature.jl`.
For the point and hypersphere features, where the localized region is not aligned with any coordinate axis, the two strategies perform similarly.
For the axis-aligned line feature at $x_1 = 1/\pi$, `HCubature.jl` has a clear advantage: its estimator detects that the error is concentrated in the $x_1$ direction and refines exclusively along that axis, while `HAdaptiveIntegration` bisects uniformly in both directions, multiplying the number of evaluations by $2^d$ at each step.
This case illustrates a known limitation of uniform subdivision, and anisotropic splitting remains a direction for future development.

![Convergence of the actual and estimated errors for the point, hypersphere, and hyperplane features on the unit square, comparing `HAdaptiveIntegration` (HAI) with `HCubature.jl`.\label{fig:cvg_rectangle}](cvg_rectangle.png)

The same comparison is repeated for the unit cube in \autoref{fig:cvg_cube}.
Here `HAdaptiveIntegration` uses the default `BerntsenEspelid` rule [@BerntsenEspelid1988] while `HCubature.jl` uses its own built-in rule, so unlike the 2D case the comparison does not isolate subdivision strategy alone — the cubature rules also differ.
Nevertheless, the qualitative picture is the same: for the point and hypersphere features the two solvers are comparable, while for the hyperplane feature `HCubature.jl` retains its advantage from estimator-driven refinement.
The cost of uniform bisection grows as $2^d$ per step, which makes the gap more pronounced in 3D and reinforces anisotropic splitting as a priority for future work.

![Convergence of the actual and estimated errors for the point, hypersphere, and hyperplane features on the unit cube, comparing `HAdaptiveIntegration` (HAI) with `HCubature.jl`.\label{fig:cvg_cube}](cvg_cube.png)

# AI usage disclosure

Generative AI was used for developing the code and to help draft and edit parts of this manuscript.
All text produced with AI assistance was reviewed, revised, and verified by the authors before inclusion in the code or paper.

<!-- # Acknowledgements -->

# References
