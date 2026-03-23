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

function bug1() # work
    I, E = @inferred(integrate(p -> sum(p) * u"A", Triangle((0, 0), (1, 0), (1, 1))))
    return I, E
end

function bug2() # does not work
    T = typeof(0.0u"A")
    dom = Triangle((0u"A", 0u"A"), (1u"A", 0u"A"), (1u"A", 1u"A"))

    # unit not working generated rule (RadonLaurie, GenzMalik, GrundmannMoeller)
    # ec = embedded_cubature(RadonLaurie(), T)
    # ec = embedded_cubature(GenzMalik{2}(), T)
    # ec = embedded_cubature(GrundmannMoeller{2}(3, 1), T)

    # uint not working in TabulatedEmbeddedCubature rule because of parse
    # ec = embedded_cubature(SQUARE_CH21, T)

    ec = embedded_cubature(RadonLaurie(), T.types[1])
    I, E = @inferred(integrate(p -> sum(p), dom; embedded_cubature=ec))
    @show I
    @show E

    return nothing
end

function bug3() # does not work
    T = typeof(0.0u"A")
    dom = Triangle((0u"A", 0u"A"), (1u"A", 0u"A"), (1u"A", 1u"A"))

    ec = embedded_cubature(RadonLaurie(), T.types[1])
    I, E = @inferred (integrate(p -> sum(p) * u"B", dom; embedded_cubature=ec))
    @show I
    @show E
end
