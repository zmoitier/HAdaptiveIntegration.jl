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

The algorithm adaptively subdivides the domain and uses embedded cubature rules to estimate
errors, targeting high accuracy with fewer function evaluations.

**Key Features:**

- Adaptive integration over simplices and orthotopes of arbitrary dimension.
- Efficient tabulated cubature rules for low-dimensional simplices and orthotopes.
- Support for custom embedded cubature rules.
- Arbitrary-precision arithmetic.

## Usage

...

### Create a domain

...

### Integrate

...

# Statement of need

<!--
A section that clearly illustrates the research purpose of the software and places it in the context of related work.
This should clearly state what problems the software is designed to solve, who the target audience is, and its relation to other work.
-->

`QuadGK.jl` [@QuadGK]

`HCubature` [@HCubature]

`Cuba.jl` [@Cuba]

# State of the field

<!--
A description of how this software compares to other commonly-used packages in the research area.
If related tools exist, provide a clear “build vs. contribute” justification explaining your unique scholarly contribution and why existing alternatives are insufficient.
-->

[@QuadGK; @HCubature]

# Software design <!-- Mathematical background -->

<!--
An explanation of the trade-offs you weighed, the design/architecture you chose, and why it matters for your research application.
This should demonstrate meaningful design thinking beyond a superficial code structure description.
-->

## Embedded cubatures

- Two quadratures
- A posteriori error estimate

## Adaptive algorithm

- Subdivision of domains
- Selection of domains to refine
- Stopping criteria

# Research impact statement

<!--
Evidence of realized impact (publications, external use, integrations) or credible near-term significance (benchmarks, reproducible materials, community-readiness signals).
The evidence should be compelling and specific, not aspirational.
-->

# Complexity estimates

...

## Uniform regularity

...

## Isotropic (nearly-)singularities

...

# Benchmarks

...

# Extended precision

...

# AI usage disclosure

<!--
Transparent disclosure of any use of generative AI in the software creation, documentation, or paper authoring.
If no AI tools were used, state this explicitly.
If AI tools were used, describe how they were used and how the quality and correctness of AI-generated content was verified.
-->

# Acknowledgements

# References
