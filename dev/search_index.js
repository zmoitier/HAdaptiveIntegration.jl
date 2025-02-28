var documenterSearchIndex = {"docs":
[{"location":"95-reference/#reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"95-reference/#Contents","page":"Reference","title":"Contents","text":"","category":"section"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Pages = [\"95-reference.md\"]","category":"page"},{"location":"95-reference/#Index","page":"Reference","title":"Index","text":"","category":"section"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Pages = [\"95-reference.md\"]","category":"page"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Modules = [HAdaptiveIntegration]","category":"page"},{"location":"95-reference/#HAdaptiveIntegration.Cuboid","page":"Reference","title":"HAdaptiveIntegration.Cuboid","text":"A axes-aligned cuboid given by two 3d-points low_corner and high_corner.\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#HAdaptiveIntegration.Domain","page":"Reference","title":"HAdaptiveIntegration.Domain","text":"abstract type Domain{D,T<:Real}\n\nAbstract type for integration domains in D dimensions with type T<:Real.\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#HAdaptiveIntegration.EmbeddedCubature","page":"Reference","title":"HAdaptiveIntegration.EmbeddedCubature","text":"struct EmbeddedCubature{H,L,D,T}\n\nAn embedded cubature rule consisting of a high order cubature rule with H nodes and a low order cubature rule with L nodes. The cubature nodes and weights are assume to be for the reference simplex or orthotope. Note that the low order cubature uses nodes[1:L] as its nodes.\n\nFields:\n\nnodes::SVector{H,SVector{D,T}}: the cubature nodes;\nweights_high::SVector{H,T}: the cubature weights for the high order cubature;\nweights_low::SVector{L,T}: the cubature weights for the low order cubature.\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#HAdaptiveIntegration.EmbeddedCubatureData","page":"Reference","title":"HAdaptiveIntegration.EmbeddedCubatureData","text":"@kwdef struct EmbeddedCubatureData\n\nAn embedded cubature rule consisting of a high order cubature rule nodes and a low order cubature rule. The cubature nodes and weights are assume to be for the reference simplex or orthotope. Note that the low order cubature uses nodes[1:L] as its nodes.\n\nFields:\n\nname::String: name of the embedded cubature;\nreference::String: where the values are found;\nnb_significant_digits::Int: number of significant digits on the node and weight values, 10^-nb_significant_digits give the relative precision of the values;\nnodes::Vector{Vector{String}}: the cubature nodes;\nweights_high::Vector{String}: the cubature weights for the high order cubature;\norder_high::Int: order of the high order cubature;\nweights_low::Vector{String}: the cubature weights for the low order cubature;\norder_low::Int: order of the low order cubature.\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#HAdaptiveIntegration.Orthotope","page":"Reference","title":"HAdaptiveIntegration.Orthotope","text":"struct Orthotope{D,T}\n\nAxes-aligned Orthotope in D dimensions given by two points low_corner and high_corner.\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#HAdaptiveIntegration.Rectangle","page":"Reference","title":"HAdaptiveIntegration.Rectangle","text":"An axes-aligned rectangle given by two 2d-points low_corner and high_corner.\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#HAdaptiveIntegration.Segment","page":"Reference","title":"HAdaptiveIntegration.Segment","text":"A segment in 1 dimensions given by two value xmin and xmax.\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#HAdaptiveIntegration.Simplex","page":"Reference","title":"HAdaptiveIntegration.Simplex","text":"struct Simplex{D,T,N}\n\nA simplex in D dimensions with N=D+1 points of type T.\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#HAdaptiveIntegration.Tetrahedron","page":"Reference","title":"HAdaptiveIntegration.Tetrahedron","text":"A tetrahedron in 3 dimensions with 4 points of type T.\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#HAdaptiveIntegration.Triangle","page":"Reference","title":"HAdaptiveIntegration.Triangle","text":"A triangle in 2 dimensions with 3 vertices of type T.\n\n\n\n\n\n","category":"type"},{"location":"95-reference/#HAdaptiveIntegration.abs_det_jacobian-Union{Tuple{HAdaptiveIntegration.Domain{D, T}}, Tuple{T}, Tuple{D}} where {D, T<:Real}","page":"Reference","title":"HAdaptiveIntegration.abs_det_jacobian","text":"abs_det_jacobian(d::Domain)\n\nThe absolute value of the Jacobian's determinant of the map from the reference domain to the physical domain d.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#HAdaptiveIntegration.check_subdivision-Union{Tuple{T}, Tuple{D}, Tuple{HAdaptiveIntegration.Domain{D, T}, Any}} where {D, T<:Real}","page":"Reference","title":"HAdaptiveIntegration.check_subdivision","text":"check_subdivision(\n    domain::Domain{D,T},\n    subdiv_algo;\n    atol::T=zero(T),\n    rtol::T=(atol > zero(T)) ? zero(T) : 10 * eps(T),\n) where {D,T<:Real}\n\nReturn nothing if the sum of the volume of the subdomain by the subdiv_algo is equal to the volume of the domain else throw an error.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#HAdaptiveIntegration.cuboid-Tuple{Any, Any}","page":"Reference","title":"HAdaptiveIntegration.cuboid","text":"cuboid(low_corner, high_corner)\n\nA axes-aligned cuboid given by two 3d-points low_corner and high_corner.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#HAdaptiveIntegration.embedded_cubature","page":"Reference","title":"HAdaptiveIntegration.embedded_cubature","text":"embedded_cubature(ecr::EmbeddedCubatureData, T=Float64)\n\nReturn the EmbeddedCubature with type T from an EmbeddedCubatureData.\n\n\n\n\n\n","category":"function"},{"location":"95-reference/#HAdaptiveIntegration.embedded_cubature-Union{Tuple{T}, Tuple{Array{Vector{T}, 1}, Vector{T}, Vector{T}}} where T<:Real","page":"Reference","title":"HAdaptiveIntegration.embedded_cubature","text":"embedded_cubature(nodes, weights_high, weights_low)\n\nReturn an embedded cubature form a vector of nodes and two vector of weights for the high order and low order cubature.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#HAdaptiveIntegration.map_from_reference-Union{Tuple{HAdaptiveIntegration.Domain{D, T}}, Tuple{T}, Tuple{D}} where {D, T<:Real}","page":"Reference","title":"HAdaptiveIntegration.map_from_reference","text":"map_from_reference(d::Domain)::Function\n\nReturn an anonymous function that maps the reference domain to the physical domain d.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#HAdaptiveIntegration.map_to_reference-Union{Tuple{HAdaptiveIntegration.Domain{D, T}}, Tuple{T}, Tuple{D}} where {D, T<:Real}","page":"Reference","title":"HAdaptiveIntegration.map_to_reference","text":"map_to_reference(d::Domain)::Function\n\nReturn an anonymous function that maps the physical domain d to the reference domain.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#HAdaptiveIntegration.orthotope-Tuple{Any, Any}","page":"Reference","title":"HAdaptiveIntegration.orthotope","text":"orthotope(low_corner, high_corner)\n\nAn axes-aligned orthotope in D dimensions given by two vectors low_corner and high_corner.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#HAdaptiveIntegration.rectangle-Tuple{Any, Any}","page":"Reference","title":"HAdaptiveIntegration.rectangle","text":"rectangle(low_corner, high_corner)\n\nAn axes-aligned rectangle given by two 2d-points low_corner and high_corner.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#HAdaptiveIntegration.reference_domain-Tuple{DataType}","page":"Reference","title":"HAdaptiveIntegration.reference_domain","text":"reference_domain(domain_type::DataType)\n\nReturn the reference domain for domain_type.\n\nReturn the reference D-simplex given by the vertices (0,...,0), (1,0,...,0), (0,1,0,...,0), (0,...,0,1).\n\nReturn the reference orthotope in D dimensions, representing [0, 1]ᴰ.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#HAdaptiveIntegration.segment-Union{Tuple{T}, Tuple{T, T}} where T<:Real","page":"Reference","title":"HAdaptiveIntegration.segment","text":"segment(xmin, xmax)\n\nA segment in 1 dimensions representing [xmin, xmax].\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#HAdaptiveIntegration.simplex-Tuple","page":"Reference","title":"HAdaptiveIntegration.simplex","text":"simplex(points...)\n\nA simplex in D dimensions with N=D+1 points of type T.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#HAdaptiveIntegration.subdivide_cuboid8-Union{Tuple{Cuboid{T}}, Tuple{T}} where T","page":"Reference","title":"HAdaptiveIntegration.subdivide_cuboid8","text":"subdivide_cuboid8(c::Cuboid)\n\nDivide the cuboid c into 8 cuboid by connecting the center of the cuboid to the midpoints of the edges.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#HAdaptiveIntegration.subdivide_rectangle4-Tuple{Rectangle}","page":"Reference","title":"HAdaptiveIntegration.subdivide_rectangle4","text":"subdivide_rectangle4(r::Rectangle)\n\nDivide the rectangle r into four squares by connecting the center of the square to the midpoints of the edges.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#HAdaptiveIntegration.subdivide_segment2-Tuple{Segment}","page":"Reference","title":"HAdaptiveIntegration.subdivide_segment2","text":"subdivide_segment2(s::Segment)\n\nDivide the segment s into two segments of equal length.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#HAdaptiveIntegration.subdivide_segment3-Tuple{Segment}","page":"Reference","title":"HAdaptiveIntegration.subdivide_segment3","text":"subdivide_segment3(s::Segment)\n\nDivide the segment s into three segments of equal length.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#HAdaptiveIntegration.subdivide_tetrahedron8-Tuple{Tetrahedron}","page":"Reference","title":"HAdaptiveIntegration.subdivide_tetrahedron8","text":"subdivide_tetrahedron8(t::Tetrahedron)\n\nDivide the tetrahedron t into eight tetrahedra by connecting the midpoints of the edges.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#HAdaptiveIntegration.subdivide_triangle2-Tuple{Triangle}","page":"Reference","title":"HAdaptiveIntegration.subdivide_triangle2","text":"subdivide_triangle2(s::Triangle)\n\nDivide the triangle t into two triangles by connecting the first point of t to the midpoints of the two other points.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#HAdaptiveIntegration.subdivide_triangle4-Tuple{Triangle}","page":"Reference","title":"HAdaptiveIntegration.subdivide_triangle4","text":"subdivide_triangle4(t::Triangle)\n\nDivide the triangle t into four triangles by connecting the midpoints of the edges.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#HAdaptiveIntegration.tetrahedron-NTuple{4, Any}","page":"Reference","title":"HAdaptiveIntegration.tetrahedron","text":"tetrahedron(a, b, c, d)\n\nA tetrahedron in 3 dimensions given by the 3d-points a, b, c, and d.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#HAdaptiveIntegration.triangle-Tuple{Any, Any, Any}","page":"Reference","title":"HAdaptiveIntegration.triangle","text":"triangle(a, b, c)\n\nA triangle in 2 dimensions given by the 2d-points a, b, and c.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = HAdaptiveIntegration","category":"page"},{"location":"#HAdaptiveIntegration","page":"Home","title":"HAdaptiveIntegration","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for HAdaptiveIntegration.","category":"page"},{"location":"#Basic-usage","page":"Home","title":"Basic usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The main function of this package is integrate, which computes the integral of a Domain","category":"page"},{"location":"","page":"Home","title":"Home","text":"using InteractiveUtils # hide\nusing HAdaptiveIntegration\nsubtypes(HAdaptiveIntegration.Domain)","category":"page"},{"location":"#Comparisson-with-HCubature.jl","page":"Home","title":"Comparisson with HCubature.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package shares many similarities with HCubature.jl; there are, however, a few important differences:","category":"page"},{"location":"","page":"Home","title":"Home","text":"HAdaptiveIntegration.jl uses tabulated embedded cubature rules such as the ones found in Cubature.jl, whereas HCubature.jl implemented the Genz-Malik algorithm valid axis-aligned rectangles in any dimension.\nHAdaptiveIntegration.jl supports simplicies in low dimensions, whereas HCubature.jl supports axis-aligned rectangles in any dimension.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Let's start with a simple example using HCubature:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using HCubature, LinearAlgebra\na, b = (0.0, 0.0), (1.0,1.0)\nconst counter = Ref(0)\n# f = x -> (counter[]+=1; 1 / (norm(x) + 1e-0))\nf = x -> (counter[]+=1; cos(20*prod(x)))\nI, E = hcubature(f, a, b)\nprintln(\"I = $I, E = $E, counter = $(counter[])\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"Now, let's do the same with HAdaptiveIntegration:","category":"page"},{"location":"","page":"Home","title":"Home","text":"import HAdaptiveIntegration as HAI\ndomain = HAI.rectangle(a, b)\ncounter[] = 0\nI, E = HAI.integrate(f, domain)\nprintln(\"I = $I, E = $E, counter = $(counter[])\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"Lets look at performance now:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using BenchmarkTools\ncounter[] = 0\nb1 = @benchmark hcubature($f, $a, $b)","category":"page"},{"location":"","page":"Home","title":"Home","text":"counter[] = 0\nb2 = @benchmark HAI.integrate($f, $domain)","category":"page"},{"location":"","page":"Home","title":"Home","text":"<!– function Base.parse(T::Type{MultiFloat{Float64,N}}, str::String) where {N}     return T(str) end –>","category":"page"}]
}
