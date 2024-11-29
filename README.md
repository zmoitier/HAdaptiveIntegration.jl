# Adaptive quadrature on simplices

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://zmoitier.github.io/AdaptiveSimplexQuadrature.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://zmoitier.github.io/AdaptiveSimplexQuadrature.jl/dev/)
[![Build Status](https://github.com/zmoitier/AdaptiveSimplexQuadrature.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/zmoitier/AdaptiveSimplexQuadrature.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/zmoitier/AdaptiveSimplexQuadrature.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/zmoitier/AdaptiveSimplexQuadrature.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

* list of schemes:
  - https://github.com/sigma-py/quadpy
  - https://github.com/zfergus/legacy-quadpy

* Possibly relevant Julia packages:
  - https://github.com/JuliaMath/Cubature.jl
  - https://github.com/JuliaMath/HCubature.jl
  - https://github.com/JuliaMath/QuadGK.jl

  - https://github.com/SciML/Integrals.jl

  - https://github.com/stevengj/cubature
  - https://github.com/eschnett/GrundmannMoeller.jl
  - https://github.com/eschnett/SimplexQuad.jl

# To-Do

1. Make it work!
2. Make it right!
3. Make it fast!

* 2-simplex (a.k.a. triangle):
  - try log singularity and nearly singular

* 3-simplex (a.k.a. tetrahedron):
  - try to do something

* n-simplex:
  - hope for the best

# Remark

* If you show the measure of the simplex in the integrate function, they seem to be way smaller than the `rtol`. Might need to add a test on the size which could speed up the computation.
  
## Related packages