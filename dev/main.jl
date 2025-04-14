using LinearAlgebra
import HAdaptiveIntegration as hai

rule = hai.embedded_cubature(Float64, hai.TRIANGLE_GM20)

const n = length(rule.weights_low)
const N = 3n

U = vcat([x[1] for x in rule.nodes[1:n]], [x[2] for x in rule.nodes[1:n]], rule.weights_low)
