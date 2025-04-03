using HAdaptiveIntegration: GenzMalik, embedded_cubature
using Printf: Format, format

ec = embedded_cubature(BigFloat, GenzMalik{3}())

fmt_node = Format("[\"%.36e\", \"%.36e\", \"%.36e\"],")
for x in ec.nodes
    println(format(fmt_node, x[1], x[2], x[3]))
end
println()

fmt_weight = Format("\"%.36e\",")
for w in ec.weights_high
    println(format(fmt_weight, w))
end
println()

for w in ec.weights_low
    println(format(fmt_weight, w))
end
println()
