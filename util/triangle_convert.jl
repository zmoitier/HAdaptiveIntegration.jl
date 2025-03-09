using StaticArrays
import Printf: Format, format

import HAdaptiveIntegration as hai

setprecision(30; base=10)

r7 = (
    nodes=[
        [big"0.3333333333333333", big"0.3333333333333333"],
        [big"0.7974269853530872", big"0.1012865073234563"],
        [big"0.1012865073234563", big"0.7974269853530872"],
        [big"0.1012865073234563", big"0.1012865073234563"],
        [big"0.0597158717897698", big"0.4701420641051151"],
        [big"0.4701420641051151", big"0.0597158717897698"],
        [big"0.4701420641051151", big"0.4701420641051151"],
    ],
    weights=[
        big"0.22500000000000000",
        big"0.12593918054482717",
        big"0.12593918054482717",
        big"0.12593918054482717",
        big"0.13239415278850616",
        big"0.13239415278850616",
        big"0.13239415278850616",
    ],
)

l19 = (
    nodes=[
        [big"0.3333333333333333", big"0.3333333333333333"],
        [big"0.7974269853530872", big"0.1012865073234563"],
        [big"0.1012865073234563", big"0.7974269853530872"],
        [big"0.1012865073234563", big"0.1012865073234563"],
        [big"0.0597158717897698", big"0.4701420641051151"],
        [big"0.4701420641051151", big"0.0597158717897698"],
        [big"0.4701420641051151", big"0.4701420641051151"],
        [big"0.5357953464498992", big"0.2321023267750504"],
        [big"0.2321023267750504", big"0.5357953464498992"],
        [big"0.2321023267750504", big"0.2321023267750504"],
        [big"0.9410382782311209", big"0.0294808608844396"],
        [big"0.0294808608844396", big"0.9410382782311209"],
        [big"0.0294808608844396", big"0.0294808608844396"],
        [big"0.7384168123405100", big"0.2321023267750504"],
        [big"0.7384168123405100", big"0.0294808608844396"],
        [big"0.2321023267750504", big"0.7384168123405100"],
        [big"0.2321023267750504", big"0.0294808608844396"],
        [big"0.0294808608844396", big"0.7384168123405100"],
        [big"0.0294808608844396", big"0.2321023267750504"],
    ],
    weights=[
        big"0.0378610912003147",
        big"0.0376204254131829",
        big"0.0376204254131829",
        big"0.0376204254131829",
        big"0.0783573522441174",
        big"0.0783573522441174",
        big"0.0783573522441174",
        big"0.1162714796569659",
        big"0.1162714796569659",
        big"0.1162714796569659",
        big"0.0134442673751655",
        big"0.0134442673751655",
        big"0.0134442673751655",
        big"0.0375097224552317",
        big"0.0375097224552317",
        big"0.0375097224552317",
        big"0.0375097224552317",
        big"0.0375097224552317",
        big"0.0375097224552317",
    ],
)

for (name, cbt, n) in [("R7", r7, 15), ("L19", l19, 15)]
    println(">> $name <<")

    local j = big"0.5"

    fmt_node = Format("[\"%.$(n)e\", \"%.$(n)e\"],")
    for x in cbt[:nodes]
        local v = SVector{2}(x)
        println(
            format(
                fmt_node,
                round(v[1]; sigdigits=n + 1, base=10),
                round(v[2]; sigdigits=n + 1, base=10),
            ),
        )
    end
    println()

    fmt_weight = Format("\"%.$(n)e\",")
    for w in cbt[:weights]
        println(format(fmt_weight, round(j * w; sigdigits=n + 1, base=10)))
    end
    println()
end
