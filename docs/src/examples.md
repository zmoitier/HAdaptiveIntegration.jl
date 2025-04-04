```@meta
CurrentModule = HAdaptiveIntegration
```

# Examples and benchmarks

Let's start with a simple example using `HCubature`:

```@example hcubature-square
using HCubature, LinearAlgebra
a, b = (0, 0), (1, 1)
const counter = Ref(0)
f = x -> (counter[]+=1; 1 / (norm(x) + 1))
I, E = hcubature(f, a, b)
println("I = $I, E = $E, counter = $(counter[])")
```

Now, let's do the same with `HAdaptiveIntegration`:

```@example hcubature-square
using HAdaptiveIntegration
domain = rectangle(a, b)
counter[] = 0
I, E = integrate(f, domain)
println("I = $I, E = $E, counter = $(counter[])")
```

Let's look at performance now:

```@example hcubature-square
using BenchmarkTools
counter[] = 0
b1 = @benchmark hcubature($f, $a, $b)
```

```@example hcubature-square
ec = HAdaptiveIntegration.embedded_cubature(Float64, HAdaptiveIntegration.SQUARE_CH25)
counter[] = 0
b2 = @benchmark integrate($f, $domain; embedded_cubature = $ec)
```

Let's do the same comparison for the 3d-cube.

```@example hcubature-cube
using HCubature, LinearAlgebra
a, b = (0, 0, 0), (1, 1, 1)
const counter = Ref(0)
f = x -> (counter[]+=1; cos(5*prod(x)))
I, E = hcubature(f, a, b)
println("I = $I, E = $E, counter = $(counter[])")
```

```@example hcubature-cube
using HAdaptiveIntegration
domain = cuboid(a, b)
counter[] = 0
I, E = integrate(f, domain)
println("I = $I, E = $E, counter = $(counter[])")
```

```@example hcubature-cube
using BenchmarkTools
counter[] = 0
b1 = @benchmark hcubature($f, $a, $b)
```

```@example hcubature-cube
ec = HAdaptiveIntegration.embedded_cubature(Float64, HAdaptiveIntegration.CUBE_BE65)
counter[] = 0
b2 = @benchmark integrate($f, $domain; embedded_cubature = $ec)
```
