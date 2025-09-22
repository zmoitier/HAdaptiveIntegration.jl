using ForwardDiff
using HAdaptiveIntegration.Domain
using HAdaptiveIntegration.Rule
using HAdaptiveIntegration: increase_precision, integrate

cbt = SQUARE_CH25
ec_f32 = embedded_cubature(cbt, Float32)
tec0 = TabulatedEmbeddedCubature{Rectangle}(;
    description="",
    reference="",
    order_high=cbt.order_high,
    order_low=cbt.order_low,
    precision=7,
    nodes=[Vector([string(v) for v in p]) for p in ec_f32.nodes],
    weights_high=string.(ec_f32.weights_high),
    weights_low=string.(ec_f32.weights_low),
)

setprecision(BigFloat, 128; base=10)
tec1 = increase_precision(tec0, BigFloat)
nothing
