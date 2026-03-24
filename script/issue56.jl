using HAdaptiveIntegration
using HAdaptiveIntegration.Domain: Rectangle, Triangle, map_from_reference
using HAdaptiveIntegration.Rule:
    EmbeddedCubature,
    GenzMalik,
    GrundmannMoeller,
    RadonLaurie,
    SQUARE_CH21,
    embedded_cubature
using StaticArrays
using Test
using Unitful

function bug1()
    I, E = @inferred(integrate(p -> 1u"g", Triangle((0, 0), (1, 0), (1, 1))))
    @show I
    @show E

    return nothing
end

function bug2()
    T = typeof(0.0u"m")

    # dom = Triangle((0u"L", 0u"L"), (1u"L", 0u"L"), (1u"L", 1u"L"))
    dom = Rectangle((0u"m", 0u"m"), (1u"m", 1u"m"))

    ec = embedded_cubature(SQUARE_CH21, typeof(one(T)))
    I, E = @inferred integrate(p -> 1, dom; embedded_cubature=ec)
    @show I
    @show E

    return nothing
end

function bug3()
    dom = Triangle((0u"m", 0u"m"), (1u"m", 0u"m"), (1u"m", 1u"m"))
    I, E = @inferred (integrate(p -> 1u"g", dom))
    @show I
    @show E

    return nothing
end
