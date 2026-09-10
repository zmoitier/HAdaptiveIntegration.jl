```@meta
CurrentModule = HAdaptiveIntegration
```

# [Extensions](@id extensions)

## `IncreasePrecisionExt`

This extension allows increasing the precision of tabulated embedded cubature. To use this extension you must add the package [`ForwardDiff`](https://juliadiff.org/ForwardDiff.jl/stable/).

```@example increase_precision
using HAdaptiveIntegration
using ForwardDiff
```

For example if we want to increase the precision of the rule [`HAdaptiveIntegration.SQUARE_CH21`](@ref), we do:

```@example increase_precision
tec0 = HAdaptiveIntegration.SQUARE_CH21
tec0.nodes
```

```@example increase_precision
tec1 = HAdaptiveIntegration.increase_precision(
    tec0, BigFloat; x_atol=big"1e-64", f_atol=big"1e-64"
)
display(tec1) # hide
```

```@example increase_precision
tec1.nodes
```

## Complete workflow: Arbitrary precision integration

Now let's use the increased precision rule in an actual integration. First define a function to integrate over a domain:

```@example increase_precision
domain = Rectangle((big"0", big"0"), (big"1", big"1"))
f = x -> x[1]^2 + x[2]^2
nothing # hide
```

Create an embedded cubature with the high-precision rule and a lower-order pair:

```@example increase_precision
using HAdaptiveIntegration.Rule: embedded_cubature
ec_big = embedded_cubature(tec1, BigFloat)
nothing # hide
```

Now integrate with arbitrary precision:

```@example increase_precision
I, E = HAdaptiveIntegration.integrate(f, domain; rule = ec_big, rtol = big"1e-64")
println("I = $I") # hide
println("E = $E") # hide
```

The result is now computed with arbitrary precision (BigFloat), enabling high-precision numerical integration when needed for sensitive applications or validation studies.

## Mathematical details

Let $\widehat{\omega}$ be a reference domain of dimension $d$ and $(\mathcal{H}, \mathcal{L})$ be an embedded cubature on $\widehat{\omega}$ with orders $k_h > k_\ell$. Denoting by $\mathsf{H}$ and $\mathsf{L}$ the number of nodes of the high and low-order rules, the cubature operators are

```math
\mathcal{H}(f) = \sum_{i=1}^{\mathsf{H}} h_i \, f(\boldsymbol{x}_i),
\qquad
\mathcal{L}(f) = \sum_{i=1}^{\mathsf{L}} \ell_i \, f(\boldsymbol{x}_i),
```

where $\boldsymbol{x}_1, \ldots, \boldsymbol{x}_{\mathsf{H}}$ are the shared nodes, $h_1, \ldots, h_{\mathsf{H}}$ the high-order weights, and $\ell_1, \ldots, \ell_{\mathsf{L}}$ the low-order weights; the low-order rule reuses the first $\mathsf{L}$ nodes. Let $K_h = \dim \mathbb{P}_{k_h}$ and $K_\ell = \dim \mathbb{P}_{k_\ell}$, and let $b_1, \ldots, b_{K_h}$ be a basis of $\mathbb{P}_{k_h}$ (polynomials of total degree $\leq k_h$) such that $b_1, \ldots, b_{K_\ell}$ is a basis of $\mathbb{P}_{k_\ell}$. Define

```math
\boldsymbol{u}^{\mathcal{H}, \mathcal{L}}
= (\boldsymbol{x}_1, \ldots, \boldsymbol{x}_{\mathsf{H}}, h_1, \ldots, h_{\mathsf{H}}, \ell_1, \ldots, \ell_{\mathsf{L}})
```

as the embedded cubature data, and the function $F \colon \mathbb{R}^{(d+1) \mathsf{H} + \mathsf{L}} \to \mathbb{R}^{K_h + K_\ell}$ by

```math
F\left( \boldsymbol{u}^{\mathcal{H}, \mathcal{L}} \right) =
\begin{pmatrix}
  \mathcal{H}(b_1) - \int_{\widehat{\omega}} b_1(\boldsymbol{x}) \operatorname{d}\!\boldsymbol{x}
  \\[1ex]
  \vdots
  \\[1ex]
  \mathcal{H}(b_{K_h}) - \int_{\widehat{\omega}} b_{K_h}(\boldsymbol{x}) \operatorname{d}\!\boldsymbol{x}
  \\[2ex]
  \mathcal{L}(b_1) - \int_{\widehat{\omega}} b_1(\boldsymbol{x}) \operatorname{d}\!\boldsymbol{x}
  \\[1ex]
  \vdots
  \\[1ex]
  \mathcal{L}(b_{K_\ell}) - \int_{\widehat{\omega}} b_{K_\ell}(\boldsymbol{x}) \operatorname{d}\!\boldsymbol{x}
\end{pmatrix}.
```

Starting from a tabulated rule with $\lVert F(\boldsymbol{u}^{\mathcal{H}, \mathcal{L}}) \rVert_2 = \varepsilon \ll 1$, we seek $\tilde{\boldsymbol{u}}$ with $\lVert F(\tilde{\boldsymbol{u}}) \rVert_2 = \eta < \varepsilon$ via a least-squares Newton method ([Xiao and Gimbutas, 2010, section 2.3](https://doi.org/10.1016/j.camwa.2009.10.027)). Setting $\boldsymbol{u}_0 = \boldsymbol{u}^{\mathcal{H}, \mathcal{L}}$, the iteration is

```math
\boldsymbol{u}_{p+1} = \boldsymbol{u}_p - \operatorname{J}_F(\boldsymbol{u}_p)^\dagger F(\boldsymbol{u}_p),
```

where $\operatorname{J}_F(\boldsymbol{u}_p)$ is the Jacobian matrix of $F$ at $\boldsymbol{u}_p$ and $\operatorname{J}_F(\boldsymbol{u}_p)^\dagger$ is the pseudo-inverse. In Julia, this is implemented with the `\` operator. The iteration stops when $\lVert \boldsymbol{u}_{p+1} - \boldsymbol{u}_p \rVert_2 \leq \mathtt{x\_atol}$ (absolute tolerance of successive iterates), or $\lVert F(\boldsymbol{u}_p) \rVert_2 \leq \mathtt{f\_atol}$ (absolute tolerance of function value), or $p = \mathtt{maxiter}$ (maximum number of iterations).
