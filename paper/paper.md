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
In `Julia`, existing packages already cover important parts of this space: `QuadGK.jl` [@QuadGK] is the natural tool for one-dimensional adaptive integration, `HCubature.jl` [@HCubature] provides adaptive cubature on axis-aligned orthotopes, and `Cuba.jl` [@Cuba] exposes a family of multidimensional integration methods, including stochastic ones, for rectangular domains.

The main contribution is that, to our knowledge, there were no `Julia` packages capable of doing adaptive integration on simplices (triangles, tetrahedra, etc.).
The second contribution is the use of specialized tabulated rules, which allow having higher order scheme than `HCubature.jl`.
We design the package in such way that there is a simple common interface for integrating over simplices and orthotopes.
And, it is easy to change the rule.
Therefore, extending the package with new rule is straightforward.
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

The `HAdaptiveIntegration` package is organized with extending it in mind.
It has two submodules: the `Domain` module which define the possible integration domain and how to subdivide them; and the `Rule` module which define the integration rules.
There is different type of rules, the explicit rules which are computed using explicit formula, and tabulated rules which are stored in a decimal format.
The package has basically one entry point the `integrate` function.
The automatic adaptive algorithm it uses has two component: an embedded cubature and a subdivision strategy described in the following sections [Embedded cubature] and [The adaptive algorithm].

## Embedded cubature

At the core of adaptive numerical integration method, there is a cubature pair $(\mathcal{H}, \mathcal{L})$, define on the reference domain $\widehat{\omega}$ (equal to $\{\boldsymbol{x} \in \mathbb{R}_+^d \mid x_1 + \cdots + x_d \leq 1\}$ for simplices or $[0, 1]^d$ for orthotopes) by
$$
  \mathcal{H}(f) = \sum_{1 \leq i \leq \mathsf{H}} h_i \, f(\boldsymbol{x}_i)
  \quad \text{and} \quad
  \mathcal{L}(f) = \sum_{1 \leq i \leq \mathsf{L}} \ell_i \, f(\boldsymbol{x}_i),
  \qquad \forall f \in \mathscr{C}^0(\widehat{\omega}).
$$
Where $\boldsymbol{x}_1, \ldots, \boldsymbol{x}_{\mathsf{H}} \in \widehat{\Omega}$ are the cubature point on the reference domain, $h_1, \ldots, h_{\mathsf{H}} \in \mathbb{R}$ are the $\mathcal{H}$ cubature weights, $\ell_1, \ldots, \ell_{\mathsf{L}} \in \mathbb{R}$ are the $\mathcal{L}$ cubature weights, and $\mathsf{H} > \mathsf{L}$ so the $\mathcal{H}$ rule has more points than the $\mathcal{L}$ rule.
The pair $(\mathcal{H}, \mathcal{L})$ is called an embedded cubature because the $\mathcal{L}$ rule use a subset of the points of the $\mathcal{H}$ rule.
In accordance with the number of point evaluation of the two cubature rule, $\mathcal{H}$ has order $k_{\mathcal{H}}$ and $\mathcal{L}$ has order $k_{\mathcal{L}}$ with $k_{\mathcal{H}} > k_{\mathcal{L}}$.
The order of a cubature rule is the highest integer $k \in \mathbb{N}$ such that the cubature is exact on the space of polynomials with total degree less or equal than $k$.

For a domain $\omega$, we define the map $\phi \colon \widehat{\omega} \to \omega$ from the reference domain to the physical domain.
Using the embedded cubature, we define the *local* estimated integral value $I_\omega$ and the *local* estimated error $E_\omega$ by 
$$
  I_\omega = \lvert\det \operatorname{J}_\phi\rvert \ \mathcal{H}(f \circ \phi)
  \quad \text{and} \quad
  E_\omega = \lvert\det \operatorname{J}_\phi\rvert \ \lVert \mathcal{H}(f \circ \phi) - \mathcal{L}(f \circ \phi) \rVert,
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
It starts with $\mathcal{P}_0 = \{\Omega\}$ then
$$
  \mathcal{P}_{n+1} = \left[ \mathcal{P}_n \setminus \left\{\omega_0\right\} \right] \cup \left\{\omega_1, \ldots, \omega_{2^d}\right \},
  \qquad \forall n \in \mathbb{N},
$$
where $\omega_0$ is chosen such that $E_{\omega_0} = \max \{E_\omega : \omega \in \mathcal{P}_n\}$, and $\omega_1, \ldots, \omega_{2^d}$ are subdomain from the subdivision of $\omega_0$.
The subdivision follows the geometry of the domain.
In dimension $d$, orthotopes are bisected along each axis, producing $2^d$ subdomains, while simplices are subdivided by midpoint edge simplex refinement, again producing $2^d$ subdomain, see [@SimplexSubdiv].

Associated to the ${(\mathcal{P}_n)}_{n \in \mathbb{N}}$ sequence, we define the global integral value $I_n$ and error $E_n$ estimators by
$$
  I_n = \sum_{\omega \in \mathcal{P}_n} I_\omega
  \quad \text{and} \quad
  E_n = \sum_{\omega \in \mathcal{P}_n} E_\omega.
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

The implementation follows closely the algorithm describe in section [The adaptive algorithm].
One important implementation detail is the uses a binary heap to store the partition $\mathcal{P}_n$.
More precisely, a max binary heap storing $\{(\omega, I_\omega, E_\omega) : \omega \in \mathcal{P}_n \}$ is used, and the comparison only involve $E_\omega$.
As this data structure allows efficient retrieval of a maximum local error, as well as efficient `pop` and `push` operation.
In addition, in the case of doing multiple integral on the same domain type and function signature, we can pre-allocate the heap and pass it to the `integrate` function via the `buffer` keyword to reduce the number of allocation.

## Extended precision

As said in the [Summary] section, `integrate` support arbitrary precision.
However, only the rules from [@GrundmannMoeller1978; @GenzMalik1980] are generated using explicit formula, all the other are tabulated with a finite decimal precision (quadruple precision) which is incompatible with arbitrary precision.
`HAdaptiveIntegration` addresses this with the optional `IncreasePrecisionExt` extension, which compute higher-precision rules from a tabulated rule by solving the polynomial exactness conditions with Newton iterations and automatic differentiation through `ForwardDiff` [@ForwardDiff2016].

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
The goal is to find a root of the function $F$ with higher precision than the stored precision, for a tabulated rule $\boldsymbol{u}^{\mathcal{H}, \mathcal{L}}$, we already have $\lVert F(\boldsymbol{u}^{\mathcal{H}, \mathcal{L}}) \rVert_2 = \varepsilon \ll 0$.
We want to find a $\tilde{\boldsymbol{u}}$ such that $\lVert F(\tilde{\boldsymbol{u}}) \rVert_2 = \eta < \varepsilon$, to do that we use a least-square Newton method, [@XiaoGimbutas2010{section 2.3}].
Set $\boldsymbol{u}_0 = \boldsymbol{u}^{\mathcal{H}, \mathcal{L}}$ and define the iteration
$$
  \boldsymbol{u}_{p+1} = \boldsymbol{u}_p - \boldsymbol{\delta}_p
  \qquad \text{where} \
  \boldsymbol{\delta}_p = \arg\min \left\{ \lVert \boldsymbol{w} \rVert_2 \mid \operatorname{J}_F(\boldsymbol{u}_p) \boldsymbol{w} = F(\boldsymbol{u}_p) \right\},\tag{N}
$$
and $\operatorname{J}_F(\boldsymbol{u}_p)$ is the Jacobian matrix of $F$ at the point $\boldsymbol{u}_p$.
The iteration stop when $\lVert \boldsymbol{u}_{p+1} - \boldsymbol{u}_p \rVert_2 \leq \mathtt{x\_atol}$ (absolute tolerance of successive iterate) or $\lVert F(\boldsymbol{u}_p) \rVert_2 \leq \mathtt{f\_atol}$ (absolute tolerance of function value) or $p = p_{\max}$ (maximum number of iterations).

**Remark.** In (N), $\operatorname{J}_F(\boldsymbol{u}_p) \boldsymbol{w} = F(\boldsymbol{u}_p)$ is solved via the least square method (just "`\`" in `Julia`).
It is written that way because, we usually have $K_{\mathcal{H}} + K_{\mathcal{L}} < (d+1) \mathsf{H} + \mathsf{L}$, so $\operatorname{J}_F(\boldsymbol{u}_p)$ is a rectangular matrix with more column than rows.
Therefore, we choose the vector in the null space of $\operatorname{J}_F(\boldsymbol{u}_p)$ with the smallest magnitude.

**Remark.** `IncreasePrecisionExt` uses the monomial basis because of the ease of uses and the exact values of the integrals is known explicitly.
However, this basis is not well condition, that means that in practice in order to have a rule at precision $\varepsilon = 10^{-m}$ you need to increase the precision using `BigFloat` with a higher precision $\varepsilon^\alpha$ with $\alpha > 1$.

**Remark.** The tabulated rules have the same symmetries as the reference domain.
However, for simplicity, we did not impose the symmetry in this method, so it does not preserve symmetry.
If a rule is tabulated with precision $\varepsilon$ and then increase the precision to $\eta < \varepsilon$ using `IncreasePrecisionExt`.
The rule will respect the symmetry up to precision $\varepsilon$ and not $\eta$.

# Research impact statement

<!--
Evidence of realized impact (publications, external use, integrations) or credible near-term significance (benchmarks, reproducible materials, community-readiness signals).
The evidence should be compelling and specific, not aspirational.
-->

`HAdaptiveIntegration` is ready for direct use.
The repository includes API documentation, advanced examples covering the pre-allocation of the buffer, the `callback` mechanism, the custom cubature, an extension for arbitrary-precision workflows, and an automated test suite spanning all supported domain.
Continuous integration checks tests, documentation, and linting.

`HAdaptiveIntegration` has been interface in [@Integrals] which is part of the SciML ecosystem.
And, it was featured in the [This month in Julia world - 2026-02](https://discourse.julialang.org/t/this-month-in-julia-world-2026-02/136110) Newsletter.

<!--
# Complexity estimates

## Uniform regularity

## Isotropic (nearly-)singularities

## Submanifold (nearly-)singularities
-->

# Benchmarks (comparison with `HCubature.jl`)

## Uniform regularity

## Isotropic (nearly-)singularities

## Sub-manifold (nearly-)singularities

Representative comparisons against `HCubature.jl` on orthotopes confirm the intended scope of the package. On a smooth integral over the unit square, $f(x) = \exp(x_1 + x_2)$ with `rtol = 1e-8`, `HAdaptiveIntegration` required 125 function evaluations versus 493 for `HCubature.jl`, and the observed runtime on a local Linux run was about 1.1 microseconds versus 3.8 microseconds. On a more difficult corner-singular square test, $f(x) = (x_1^2 + x_2^2 + 10^{-8})^{-1/2}$ with `rtol = 1e-6`, the counts were 2125 versus 2941 evaluations, again in favor of `HAdaptiveIntegration`.

The same trend appears in three dimensions for low-dimensional orthotopes with specialized rules: on the unit cuboid with $f(x) = \exp(x_1 + x_2 + x_3)$ and `rtol = 1e-6`, `HAdaptiveIntegration` used 65 evaluations versus 429 for `HCubature.jl`. However, the crossover expected from the state-of-the-field discussion is also visible: on a smooth four-dimensional orthotope with `rtol = 1e-4`, `HCubature.jl` used fewer evaluations (171 versus 969) and was faster on the same machine. These measurements support the intended positioning of `HAdaptiveIntegration`: low-dimensional simplices and orthotopes, especially when geometry-specific tabulated rules are available, while `HCubature.jl` remains a strong choice for higher-dimensional boxes.

# AI usage disclosure

Generative AI was used for developing the code and to help draft and edit parts of this manuscript.
All text produced with AI assistance was reviewed, revised, and verified by the authors before inclusion in the code or paper.

<!-- # Acknowledgements -->

# References
