# TODO: extended precision version
"""
    const TRIANGLE_RADON_O5_N7

Reference:
https://www.math.unipd.it/~alvise/SETS_CUBATURE_TRIANGLE/radon/set_radon_standard.m
"""
const TRIANGLE_RADON_O5_N7 = Quadrature(
    [
        SVector(0.3333333333333333, 0.3333333333333333),
        SVector(0.7974269853530872, 0.1012865073234563),
        SVector(0.1012865073234563, 0.7974269853530872),
        SVector(0.1012865073234563, 0.1012865073234563),
        SVector(0.0597158717897698, 0.4701420641051151),
        SVector(0.4701420641051151, 0.0597158717897698),
        SVector(0.4701420641051151, 0.4701420641051151),
    ],
    [
        0.22500000000000000 / 2,
        0.12593918054482717 / 2,
        0.12593918054482717 / 2,
        0.12593918054482717 / 2,
        0.13239415278850616 / 2,
        0.13239415278850616 / 2,
        0.13239415278850616 / 2,
    ],
    "TRIANGLE_RADON_O5_N7",
)

# TODO: extended precision version
"""
    const TRIANGLE_LAURIE_O8_N19

Reference:
https://www.math.unipd.it/~alvise/SETS_CUBATURE_TRIANGLE/laurie/set_laurie_standard.m
"""
const TRIANGLE_LAURIE_O8_N19 = Quadrature(
    [
        SVector(0.3333333333333333, 0.3333333333333333),
        SVector(0.7974269853530872, 0.1012865073234563),
        SVector(0.1012865073234563, 0.7974269853530872),
        SVector(0.1012865073234563, 0.1012865073234563),
        SVector(0.0597158717897698, 0.4701420641051151),
        SVector(0.4701420641051151, 0.0597158717897698),
        SVector(0.4701420641051151, 0.4701420641051151),
        #
        SVector(0.5357953464498992, 0.2321023267750504),
        SVector(0.2321023267750504, 0.5357953464498992),
        SVector(0.2321023267750504, 0.2321023267750504),
        SVector(0.9410382782311209, 0.0294808608844396),
        SVector(0.0294808608844396, 0.9410382782311209),
        SVector(0.0294808608844396, 0.0294808608844396),
        SVector(0.7384168123405100, 0.2321023267750504),
        SVector(0.7384168123405100, 0.0294808608844396),
        SVector(0.2321023267750504, 0.7384168123405100),
        SVector(0.2321023267750504, 0.0294808608844396),
        SVector(0.0294808608844396, 0.7384168123405100),
        SVector(0.0294808608844396, 0.2321023267750504),
    ],
    [
        0.0378610912003147 / 2,
        0.0376204254131829 / 2,
        0.0376204254131829 / 2,
        0.0376204254131829 / 2,
        0.0783573522441174 / 2,
        0.0783573522441174 / 2,
        0.0783573522441174 / 2,
        0.1162714796569659 / 2,
        0.1162714796569659 / 2,
        0.1162714796569659 / 2,
        0.0134442673751655 / 2,
        0.0134442673751655 / 2,
        0.0134442673751655 / 2,
        0.0375097224552317 / 2,
        0.0375097224552317 / 2,
        0.0375097224552317 / 2,
        0.0375097224552317 / 2,
        0.0375097224552317 / 2,
        0.0375097224552317 / 2,
    ],
    "TRIANGLE_RADON_O5_N7",
)
