# Radon' 7 point rule of order 5 for triangle from
# https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
const TRIANGLE_R5N7 = (
    nodes = [
        SVector(0.33333333333333333, 0.33333333333333333),
        SVector(0.79742698535308720, 0.10128650732345633),
        SVector(0.10128650732345633, 0.79742698535308720),
        SVector(0.10128650732345633, 0.10128650732345633),
        SVector(0.05971587178976981, 0.47014206410511505),
        SVector(0.47014206410511505, 0.05971587178976981),
        SVector(0.47014206410511505, 0.47014206410511505),
    ],
    weights = [
        0.22500000000000000 / 2,
        0.12593918054482717 / 2,
        0.12593918054482717 / 2,
        0.12593918054482717 / 2,
        0.13239415278850616 / 2,
        0.13239415278850616 / 2,
        0.13239415278850616 / 2,
    ],
)

# Laurie's 19 point rule of order 8 from
# https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
const TRIANGLE_L8N19 = (
    SVector{19}(
        SVector(0.3333333333333333, 0.3333333333333333),
        SVector(0.7974269853530872, 0.1012865073234563),
        SVector(0.1012865073234563, 0.7974269853530872),
        SVector(0.1012865073234563, 0.1012865073234563),
        SVector(0.0597158717897698, 0.4701420641051151),
        SVector(0.4701420641051151, 0.0597158717897698),
        SVector(0.4701420641051151, 0.4701420641051151),
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
    ),
    SVector{19}(
        0.0378610912003147,
        0.0376204254131829,
        0.0376204254131829,
        0.0376204254131829,
        0.0783573522441174,
        0.0783573522441174,
        0.0783573522441174,
        0.1162714796569659,
        0.1162714796569659,
        0.1162714796569659,
        0.0134442673751655,
        0.0134442673751655,
        0.0134442673751655,
        0.0375097224552317,
        0.0375097224552317,
        0.0375097224552317,
        0.0375097224552317,
        0.0375097224552317,
        0.0375097224552317,
    ) / 2,
)
