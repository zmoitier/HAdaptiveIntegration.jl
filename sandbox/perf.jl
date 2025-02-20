using HCubature, LinearAlgebra
import HAdaptiveIntegration as HAI

a, b = (0.0, 0.0), (1.0, 1.0)
f = x -> cos(200 * prod(x))
buffer = hcubature_buffer(f, a, b)
I, E = hcubature(f, a, b; buffer=buffer)
@btime hcubature($f, $a, $b; buffer=$buffer)

domain = HAI.Square(a, b)
heap = HAI.allocate_buffer(f, domain)
I, E = HAI.integrate(f, domain; heap)
@btime HAI.integrate($f, $domain; heap=$heap)
