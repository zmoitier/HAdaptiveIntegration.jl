using StaticArrays
import Printf: Format, format

import HAdaptiveIntegration as hai

setprecision(BigFloat, 40; base=10)

ec = hai.embedded_cubature(BigFloat, hai.GrundmannMoeller(3, 7))
n = 33

fmt_node = Format("[\"%.$(n)e\", \"%.$(n)e\", \"%.$(n)e\"],")
for x in ec.nodes
    println(
        format(
            fmt_node,
            round(x[1]; sigdigits=n + 1, base=10),
            round(x[2]; sigdigits=n + 1, base=10),
            round(x[3]; sigdigits=n + 1, base=10),
        ),
    )
end
println()

fmt_weight = Format("\"%.$(n)e\",")
for w in ec.weights_high
    println(format(fmt_weight, round(w; sigdigits=n + 1, base=10)))
end
println()

for w in ec.weights_low
    println(format(fmt_weight, round(w; sigdigits=n + 1, base=10)))
end
println()
