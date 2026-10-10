# [Citing and contributing](@id citing-contributing)

## [Citing](@id citing)

If you use this package in your research, please cite our article in the Journal of Open Source Software ([doi:10.21105/joss.10833](https://doi.org/10.21105/joss.10833)):

```bibtex
@article{FariaMoitier2026,
  doi = {10.21105/joss.10833},
  url = {https://doi.org/10.21105/joss.10833},
  year = {2026},
  publisher = {The Open Journal},
  volume = {11},
  number = {126},
  pages = {10833},
  author = {Faria, Luiz and Moitier, Zoïs},
  title = {HAdaptiveIntegration.jl: Adaptive numerical integration over simplices and orthotopes},
  journal = {Journal of Open Source Software}
}
```

## [Contributing](@id contributing)

Contributions of all kinds are welcome: bug reports, documentation fixes, new cubature rules, performance improvements, and new features.
See [`CONTRIBUTING.md`](https://github.com/zmoitier/HAdaptiveIntegration.jl/blob/main/CONTRIBUTING.md) for the full guidelines.

- **Report a bug:** open an issue on the [issue tracker](https://github.com/zmoitier/HAdaptiveIntegration.jl/issues) with a minimal, self-contained example and the output of `versioninfo()`.
  For numerical issues, also give the integrand, the domain, the requested tolerances, and the value you expect.
- **Suggest a feature:** open an issue describing the use case before writing code, so the API and the location of the feature can be agreed on.
- **Ask a question:** unclear documentation is a documentation bug, so please open an issue.
  General questions are also welcome on the [Julia Discourse](https://discourse.julialang.org) forum.
- **Submit a pull request:** fork the repository, run the test suite with `julia --project=. -e 'using Pkg; Pkg.test()'`, and open a pull request.
