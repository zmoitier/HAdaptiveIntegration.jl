# HAdaptiveIntegration

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://zmoitier.github.io/HAdaptiveIntegration.jl/stable)
[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://zmoitier.github.io/HAdaptiveIntegration.jl/dev)
[![Test workflow status](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Lint workflow Status](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/zmoitier/HAdaptiveIntegration.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/zmoitier/HAdaptiveIntegration.jl)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

* list of schemes:
  * <https://github.com/sigma-py/quadpy>
  * <https://github.com/zfergus/legacy-quadpy>

* Possibly relevant Julia packages:
  * <https://github.com/JuliaMath/Cubature.jl>
  * <https://github.com/JuliaMath/HCubature.jl>
  * <https://github.com/JuliaMath/QuadGK.jl>

  * <https://github.com/SciML/Integrals.jl>

  * <https://github.com/stevengj/cubature>
  * <https://github.com/eschnett/GrundmannMoeller.jl>
  * <https://github.com/eschnett/SimplexQuad.jl>

## To-Do

1. Make it work!
2. Make it right!
3. Make it fast!

* 2-simplex (a.k.a. triangle):
  * try log singularity and nearly singular

* 3-simplex (a.k.a. tetrahedron):
  * try to do something

* n-simplex:
  * hope for the best

## Remark

* If you show the measure of the simplex in the integrate function, they seem to be way smaller than the `rtol`. Might need to add a test on the size which could speed up the computation.

## Related packages
