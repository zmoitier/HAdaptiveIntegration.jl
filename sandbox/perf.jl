using HCubature, LinearAlgebra
import HAdaptiveIntegration as ASQ

a, b = (0.0, 0.0), (1.0, 1.0)
f = x -> cos(200 * prod(x))
buffer = hcubature_buffer(f, a, b)
I, E = hcubature(f, a, b; buffer=buffer)
@btime hcubature($f, $a, $b; buffer=$buffer)

domain = ASQ.Square(a, b)
heap = ASQ.allocate_buffer(f, domain)
I, E = ASQ.integrate(f, domain; heap)
@btime ASQ.integrate($f, $domain; heap=$heap)
