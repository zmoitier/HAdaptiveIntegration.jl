# Contributing to HAdaptiveIntegration.jl

Thanks for your interest in contributing! `HAdaptiveIntegration.jl` is an open-source
Julia package, and contributions of all kinds are welcome: bug reports, documentation
fixes, new cubature rules, performance improvements, and new features.

By participating in this project, you agree to interact respectfully and
constructively with other contributors.

## Ways to contribute

### Reporting a bug

Please open an issue on the
[issue tracker](https://github.com/zmoitier/HAdaptiveIntegration.jl/issues) and include:

- what you expected to happen and what happened instead;
- a minimal, self-contained code example that reproduces the problem;
- the output of `versioninfo()` and the version of `HAdaptiveIntegration` you are using
  (`] status HAdaptiveIntegration`).

For numerical issues, it helps a lot to state the integrand, the domain, the requested
tolerances, and the value you believe is correct (and how you obtained it).

### Suggesting a feature

Open an issue describing the use case before writing code. This avoids duplicated work
and lets us agree on the API and on where the feature belongs. Good candidates include
new embedded cubature rules, support for additional domain types, and interoperability
with other packages.

### Asking a question

If something in the documentation is unclear, that is a documentation bug — please open
an issue. General questions about numerical integration in Julia are also welcome on
the [Julia Discourse](https://discourse.julialang.org) forum.

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
quality checks via [Aqua.jl](https://github.com/JuliaTesting/Aqua.jl) in addition to the
numerical tests.

### Building the documentation

The documentation is built with
[Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) from the sources in
`docs/src`:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The rendered HTML is written to `docs/build`.

## Making a change

1. Create a branch for your work:

   ```sh
   git switch -c my-descriptive-branch-name
   ```

2. Make your changes, adding tests that cover them and updating the docstrings and
   documentation as needed.

3. Make sure the full test suite passes locally.

4. Push your branch and open a pull request against `main`. Describe what the change
   does and why; link any related issue.

Small, focused pull requests are easier to review and get merged faster. If you are
planning a large change, please open an issue first to discuss it.

### Code style and formatting

Julia code in this repository is formatted with
[Runic.jl](https://github.com/fredrikekre/Runic.jl), and other files are checked by a
set of [pre-commit](https://pre-commit.com) hooks (TOML/YAML validation, markdown
linting, trailing whitespace, line endings). The same checks run in CI, so it is
convenient to install the hooks locally once:

```sh
pip install pre-commit   # or: pipx install pre-commit
pre-commit install
```

Then formatting is applied automatically on each commit. You can also run all hooks on
the whole repository at any time:

```sh
pre-commit run --all-files
```

Beyond formatting, please try to match the surrounding code: exported functions and
types carry docstrings, and new numerical routines should be generic in the element type
where reasonable.

### Adding a cubature rule

New embedded cubature rules are especially welcome. When contributing one, please also
provide:

- the source of the nodes and weights (paper, technical report, or the script used to
  generate them), added under `script/` if applicable;
- the polynomial degree of both the high- and low-order rule;
- tests verifying that the rule integrates monomials up to its stated degree exactly (to
  within round-off), following the existing tests in `test/test_rule.jl`.

## Continuous integration

Every pull request runs the test suite on the current stable Julia release, the linting
and formatting checks, and a documentation build. These are all expected to pass before
a pull request is merged. Maintainers may ask for changes; this is a normal part of the
review process and not a judgement on your work.

## Licence

`HAdaptiveIntegration.jl` is distributed under the terms of the licence in
[LICENSE](LICENSE). By contributing, you agree that your contributions will be licensed
under the same terms.
