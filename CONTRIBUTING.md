# Contributing to HAdaptiveIntegration.jl

Thanks for contributing to `HAdaptiveIntegration.jl`. Contributions of all kinds are
welcome: bug reports, documentation fixes, new cubature rules, performance
improvements, and new features. By participating, you agree to interact respectfully
and constructively with other contributors.

## Ways to contribute

### Reporting a bug

Please open an issue on the
[issue tracker](https://github.com/zmoitier/HAdaptiveIntegration.jl/issues) with:

- what you expected to happen and what happened instead;
- a minimal, self-contained example reproducing the problem;
- the output of `versioninfo()` and your `HAdaptiveIntegration` version
  (`] status HAdaptiveIntegration`).

For numerical issues, also state the integrand, the domain, the requested tolerances,
and the value you believe is correct (and how you obtained it).

### Suggesting a feature

Open an issue describing the use case before writing code, so we can agree on the API
and where the feature belongs. Good candidates: new embedded cubature rules, support
for additional domain types, and interoperability with other packages.

### Asking a question

If the documentation is unclear, that is a documentation bug — please open an issue.
General questions about numerical integration in Julia are also welcome on the
[Julia Discourse](https://discourse.julialang.org) forum.

## Development setup

1. Fork the repository on GitHub and clone your fork:

   ```sh
   git clone https://github.com/<your-username>/HAdaptiveIntegration.jl
   cd HAdaptiveIntegration.jl
   ```

2. Instantiate the package environment:

   ```sh
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   ```

3. Run the test suite to check that everything works before you change anything:

   ```sh
   julia --project=. -e 'using Pkg; Pkg.test()'
   ```

The test suite lives in `test/` and is driven by `test/runtests.jl`. It includes code
quality checks via [Aqua.jl](https://github.com/JuliaTesting/Aqua.jl), numerical tests,
and extension tests under `test/test_ext_incr_prec.jl` for the `ForwardDiff`-based
`IncreasePrecisionExt`. To run a single file: `julia --project=. -e
'include("test/test_rule.jl")'`.

### Building the documentation

The documentation is built with
[Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) from the sources in
`docs/src`:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The rendered HTML is written to `docs/build`; `make.jl` only builds locally and does
not deploy.

## Making a change

1. Create a branch: `git switch -c my-descriptive-branch-name`.

2. Make your changes, adding tests that cover them and updating docstrings and
   documentation as needed.

3. Make sure the full test suite passes locally.

4. Push your branch and open a pull request against `main`. Describe what the change
   does and why; link the related issue. For non-trivial changes, open an issue first
   to discuss it. Small, focused pull requests are easier to review and merge.

### Code style and formatting

Julia code is formatted with [Runic.jl](https://github.com/fredrikekre/Runic.jl), and
other files are checked by [pre-commit](https://pre-commit.com) hooks (see
`.pre-commit-config.yaml` for the full list: TOML/YAML/JSON validation, markdown
linting, trailing whitespace, line endings). The same checks run in CI, so install the
hooks locally:

```sh
pip install pre-commit   # or: pipx install pre-commit
pre-commit install
```

Formatting is then applied automatically on each commit. To run all hooks on the whole
repository:

```sh
pre-commit run --all-files
```

Beyond formatting, match the surrounding code: exported functions and types carry
docstrings, and new numerical routines should be generic in the element type where
reasonable.

### Adding a cubature rule

New embedded cubature rules are especially welcome. When contributing one, also
provide:

- the source of the nodes and weights (paper, technical report, or the script that
  generated them), under `script/` if applicable;
- the polynomial degree of both the high- and low-order rule;
- tests verifying that the rule integrates monomials up to its stated degree exactly
  (to round-off), following the existing tests in `test/test_rule.jl`.

### Performance improvements

Performance changes are welcome but should come with evidence. Compare before and
after with [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl) and report
timings and allocations (`@benchmark` or `BenchmarkTools.@btime`). State the workload,
the Julia version, and whether any change trades accuracy for speed.

## Continuous integration

Every pull request runs the test suite on Linux (`ubuntu-latest`, x64) for the current
stable Julia release, plus linting, formatting, and a documentation build. All must
pass before merging. Maintainers may ask for changes; this is a normal part of review.

## License

`HAdaptiveIntegration.jl` is distributed under the terms of the license in
[LICENSE](LICENSE). By contributing, you agree that your contributions will be licensed
under the same terms.
