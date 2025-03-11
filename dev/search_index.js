var documenterSearchIndex = {"docs":
[{"location":"examples/","page":"Examples and benchmarks","title":"Examples and benchmarks","text":"CurrentModule = HAdaptiveIntegration","category":"page"},{"location":"examples/#Examples-and-benchmarks","page":"Examples and benchmarks","title":"Examples and benchmarks","text":"","category":"section"},{"location":"examples/","page":"Examples and benchmarks","title":"Examples and benchmarks","text":"Let's start with a simple example using HCubature:","category":"page"},{"location":"examples/","page":"Examples and benchmarks","title":"Examples and benchmarks","text":"using HCubature, LinearAlgebra\na, b = (0.0, 0.0), (1.0,1.0)\nconst counter = Ref(0)\nf = x -> (counter[]+=1; 1 / (norm(x) + 1e-0))\n# f = x -> (counter[]+=1; cos(20*prod(x)))\nI, E = hcubature(f, a, b)\nprintln(\"I = $I, E = $E, counter = $(counter[])\")","category":"page"},{"location":"examples/","page":"Examples and benchmarks","title":"Examples and benchmarks","text":"Now, let's do the same with HAdaptiveIntegration:","category":"page"},{"location":"examples/","page":"Examples and benchmarks","title":"Examples and benchmarks","text":"using HAdaptiveIntegration\ndomain = rectangle(a, b)\ncounter[] = 0\nI, E = integrate(f, domain)\nprintln(\"I = $I, E = $E, counter = $(counter[])\")","category":"page"},{"location":"examples/","page":"Examples and benchmarks","title":"Examples and benchmarks","text":"Lets look at performance now:","category":"page"},{"location":"examples/","page":"Examples and benchmarks","title":"Examples and benchmarks","text":"using BenchmarkTools\ncounter[] = 0\nb1 = @benchmark hcubature($f, $a, $b)","category":"page"},{"location":"examples/","page":"Examples and benchmarks","title":"Examples and benchmarks","text":"ec = HAdaptiveIntegration.embedded_cubature(HAdaptiveIntegration.SQUARE_CHG25, Float64)\ncounter[] = 0\nb2 = @benchmark integrate($f, $domain; embedded_cubature = $ec)","category":"page"},{"location":"examples/","page":"Examples and benchmarks","title":"Examples and benchmarks","text":"Let's do the same comparison for the 3d-cube.","category":"page"},{"location":"examples/","page":"Examples and benchmarks","title":"Examples and benchmarks","text":"using HCubature, LinearAlgebra\na, b = (0.0, 0.0, 0.0), (1.0, 1.0, 1.0)\nconst counter = Ref(0)\nf = x -> (counter[]+=1; cos(5*prod(x)))\nI, E = hcubature(f, a, b)\nprintln(\"I = $I, E = $E, counter = $(counter[])\")","category":"page"},{"location":"examples/","page":"Examples and benchmarks","title":"Examples and benchmarks","text":"using HAdaptiveIntegration\ndomain = cuboid(a, b)\ncounter[] = 0\nI, E = integrate(f, domain)\nprintln(\"I = $I, E = $E, counter = $(counter[])\")","category":"page"},{"location":"examples/","page":"Examples and benchmarks","title":"Examples and benchmarks","text":"using BenchmarkTools\ncounter[] = 0\nb1 = @benchmark hcubature($f, $a, $b)","category":"page"},{"location":"examples/","page":"Examples and benchmarks","title":"Examples and benchmarks","text":"ec = HAdaptiveIntegration.embedded_cubature(HAdaptiveIntegration.CUBE_BE65, Float64)\ncounter[] = 0\nb2 = @benchmark integrate($f, $domain; embedded_cubature = $ec)","category":"page"},{"location":"advanced/","page":"Advanced usage","title":"Advanced usage","text":"CurrentModule = HAdaptiveIntegration","category":"page"},{"location":"advanced/#advanced-usage","page":"Advanced usage","title":"Advanced usage","text":"","category":"section"},{"location":"advanced/","page":"Advanced usage","title":"Advanced usage","text":"We now cover the options available for the integrate function.","category":"page"},{"location":"advanced/#Buffering","page":"Advanced usage","title":"Buffering","text":"","category":"section"},{"location":"advanced/","page":"Advanced usage","title":"Advanced usage","text":"When calling integrate(f, domain), the package allocates memory for storing the various subregions that are generated during the adaptive integration process. Here is what it looks like in practice:","category":"page"},{"location":"advanced/","page":"Advanced usage","title":"Advanced usage","text":"using HAdaptiveIntegration\nusing BenchmarkTools\nt = triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))\nf = x -> 1 / (x[1]^2 + x[2]^2 + 1e-2)\n@benchmark integrate($f, $t)","category":"page"},{"location":"advanced/","page":"Advanced usage","title":"Advanced usage","text":"While the overhead associated to these (small) allocations are usually negligible, there are circumstances where one may want to avoid allocations altogether. This can achieved by passing a buffer to the integrate using allocate_buffer:","category":"page"},{"location":"advanced/","page":"Advanced usage","title":"Advanced usage","text":"using HAdaptiveIntegration: allocate_buffer, integrate, triangle\nbuffer = allocate_buffer(f, t)\nintegrate(f,t; buffer)\n@benchmark integrate($f, $t; buffer = $buffer)","category":"page"},{"location":"advanced/#Embedded-cubature-formulas","page":"Advanced usage","title":"Embedded cubature formulas","text":"","category":"section"},{"location":"advanced/","page":"Advanced usage","title":"Advanced usage","text":"By default, when calling integrate(f, domain), the package uses a default embedded cubature formula for the given domain by calling default_embedded_cubature. Although these are generally good choices, you can also specify a custom embedded cubature formula by passing it as the third argument to integrate. For example, in the case of a triangle, the package defaults to a Radon-Laurie embedded cubature formula of orders 5 and 8 (see the rules_simplex.jl file for more details). If you want e.g. to use an embedded cubature based on the GrundmannMoeller rule of order 13, you can do","category":"page"},{"location":"advanced/","page":"Advanced usage","title":"Advanced usage","text":"using HAdaptiveIntegration: GrundmannMoeller, embedded_cubature, integrate, triangle\nt = triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))\nf = x -> 1 / (x[1]^2 + x[2]^2 + 1e-2)\nec = embedded_cubature(GrundmannMoeller(2,13), Float64)\nI,E = integrate(f, t; embedded_cubature = ec)","category":"page"},{"location":"advanced/","page":"Advanced usage","title":"Advanced usage","text":"tip: Available embedded cubature formulas\nYou can find a list of available embedded cubature formulas in the rule_simplex.jl file and rule_orthotope.jl files.","category":"page"},{"location":"advanced/","page":"Advanced usage","title":"Advanced usage","text":"To add a custom embedded quadrature for a given domain, you must write a constructor e.g. my_custom_cubature(args...) that returns a valid EmbeddedCubature object (see embedded_cubature.jl for some examples on how this is done). PRs with new schemes are more than welcome!","category":"page"},{"location":"advanced/#Subdivision-strategies","page":"Advanced usage","title":"Subdivision strategies","text":"","category":"section"},{"location":"advanced/","page":"Advanced usage","title":"Advanced usage","text":"The package uses a default subdivision strategy for the given domain by calling default_subdivision. For example, by default triangles are subidivided into 4 smaller triangles by connecting the midpoints of the edges:","category":"page"},{"location":"advanced/","page":"Advanced usage","title":"Advanced usage","text":"using HAdaptiveIntegration\nt = triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))\nsubdiv_algo = HAdaptiveIntegration.default_subdivision(t)","category":"page"},{"location":"advanced/","page":"Advanced usage","title":"Advanced usage","text":"Here are the subdivided triangles:","category":"page"},{"location":"advanced/","page":"Advanced usage","title":"Advanced usage","text":"subdiv_algo(t)","category":"page"},{"location":"advanced/","page":"Advanced usage","title":"Advanced usage","text":"But is is also possible (and maybe desirable) to split the triangle into 2 smaller triangles instead. The following function accomplishes this:","category":"page"},{"location":"advanced/","page":"Advanced usage","title":"Advanced usage","text":"using HAdaptiveIntegration: subdivide_triangle2\nsubdivide_triangle2(t)","category":"page"},{"location":"advanced/","page":"Advanced usage","title":"Advanced usage","text":"Passing subdivide_triangle2 as the subdiv_algo to integrate will use this instead of the default:","category":"page"},{"location":"advanced/","page":"Advanced usage","title":"Advanced usage","text":"f = x -> 1 / (x[1]^2 + x[2]^2 + 1e-2)\nI, E = integrate(f, t; subdiv_algo = subdivide_triangle2)","category":"page"},{"location":"advanced/","page":"Advanced usage","title":"Advanced usage","text":"Which subdivision strategy is best depends on the function being integrated; for the example presented above, it turns out the default strategy is more efficient!","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = HAdaptiveIntegration","category":"page"},{"location":"#HAdaptiveIntegration","page":"Home","title":"HAdaptiveIntegration","text":"","category":"section"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"HAdaptiveIntegration.jl is a Julia package for for approximating integrals of functions over various predefined Domains. It uses embedded quadrature rules to build error estimates, and refines the integration domain by splitting its mesh elements until a certain tolerance is reached. Features include:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Adaptive integration over simplices of any dimension\nUse of efficient (tabulated) cubatures for low-dimensional cuboids and simplices\nSupport for custom cubature rules\nArbitrary precision arithmetic","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package can be installed with the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add HAdaptiveIntegration","category":"page"},{"location":"","page":"Home","title":"Home","text":"Or, equivalently, via the Julia Pkg API","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> import Pkg; Pkg.add(\"HAdaptiveIntegration\")","category":"page"},{"location":"#Basic-usage","page":"Home","title":"Basic usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The main function exported by this package is integrate(f,Ω), which is used to approximate","category":"page"},{"location":"","page":"Home","title":"Home","text":"I = int_Omega f(x)  dx","category":"page"},{"location":"","page":"Home","title":"Home","text":"where Omega subset mathbbR^d is a Domain object, and f  mathbbR^d to mathbbF is a function/functor. Here is a simple example:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using HAdaptiveIntegration\ndomain = triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))\nf      = x -> 1 / (x[1]^2 + x[2]^2 + 1e-2)\nI, E   = integrate(f, domain)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The result I is the integral of f over a triangle with vertices (0,0), (1,0), and (0,1), and E is an error estimate.","category":"page"},{"location":"","page":"Home","title":"Home","text":"warning: Function signature\nThe function f must accept a single argument x which is an abstract vector of length d, the dimension of the domain (concretely, f is called through f(::SVector)). The return type T of f can be any type that supports the operations +(T,T), norm(T), and multiplication by a scalar (e.g. vectors or matrices).","category":"page"},{"location":"","page":"Home","title":"Home","text":"The keyword arguments atol and rtol can be used to control the desired absolute and relative error tolerances, respectively:","category":"page"},{"location":"","page":"Home","title":"Home","text":"I, E   = integrate(f, domain; rtol = 1e-12)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Finally, to integrate the same function over a three-dimensional axis-aligned cuboid we can use","category":"page"},{"location":"","page":"Home","title":"Home","text":"domain = cuboid((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))\nI, E   = integrate(f, domain)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Domains are constructed using the following functions (see their respective docstrings for more details):","category":"page"},{"location":"","page":"Home","title":"Home","text":"triangle: a triangle in 2D\ntetrahedron: a tetrahedron in 3D\nsimplex: a simplex in arbitrary dimension\nrectangle: a rectangle in 2D\ncuboid: a cuboid in 3D","category":"page"},{"location":"","page":"Home","title":"Home","text":"tip: N-cuboids and `HCubature.jl`\nIf you are looking for a package that supports adaptive integration over artbitrarily high-dimensional axis-aligned cuboids, you may want to check out HCubature.jl.","category":"page"},{"location":"#Going-further","page":"Home","title":"Going further","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"In the previous examples we covered the basic usage of the integrate function. There are, however, other options that can be passed to integrate in order to customize various aspects of the underlying algorithm (e.g. passing a different cubature rule, using a buffer to avoid memory allocations, etc.). For more details, see the docstring of the integrate function, as well as the next section on advanced usage.","category":"page"},{"location":"docstrings/#docstrings","page":"Docstrings","title":"Docstrings","text":"","category":"section"},{"location":"docstrings/","page":"Docstrings","title":"Docstrings","text":"Modules = [HAdaptiveIntegration]","category":"page"},{"location":"docstrings/#HAdaptiveIntegration.CUBE_BE65","page":"Docstrings","title":"HAdaptiveIntegration.CUBE_BE65","text":"Berntsen and Espelid with 65 nodes.\n\n\n\n\n\n","category":"constant"},{"location":"docstrings/#HAdaptiveIntegration.SEGMENT_GK15","page":"Docstrings","title":"HAdaptiveIntegration.SEGMENT_GK15","text":"Gauss-Kronrod with 15 nodes.\n\n\n\n\n\n","category":"constant"},{"location":"docstrings/#HAdaptiveIntegration.SEGMENT_GK31","page":"Docstrings","title":"HAdaptiveIntegration.SEGMENT_GK31","text":"Gauss-Kronrod with 31 nodes.\n\n\n\n\n\n","category":"constant"},{"location":"docstrings/#HAdaptiveIntegration.SQUARE_CHG25","page":"Docstrings","title":"HAdaptiveIntegration.SQUARE_CHG25","text":"Cools and Haegemans with 25 nodes.\n\n\n\n\n\n","category":"constant"},{"location":"docstrings/#HAdaptiveIntegration.TETRAHEDRON_GM35","page":"Docstrings","title":"HAdaptiveIntegration.TETRAHEDRON_GM35","text":"Grundmann and Möller with 35 nodes.\n\n\n\n\n\n","category":"constant"},{"location":"docstrings/#HAdaptiveIntegration.TRIANGLE_GM20","page":"Docstrings","title":"HAdaptiveIntegration.TRIANGLE_GM20","text":"Grundmann and Möller with 20 nodes.\n\n\n\n\n\n","category":"constant"},{"location":"docstrings/#HAdaptiveIntegration.TRIANGLE_RL19","page":"Docstrings","title":"HAdaptiveIntegration.TRIANGLE_RL19","text":"Laurie Radon with 19 nodes.\n\n\n\n\n\n","category":"constant"},{"location":"docstrings/#HAdaptiveIntegration.Cuboid","page":"Docstrings","title":"HAdaptiveIntegration.Cuboid","text":"Cuboid{T} = Orthotope{3,T}\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#HAdaptiveIntegration.Domain","page":"Docstrings","title":"HAdaptiveIntegration.Domain","text":"abstract type Domain{D,T<:Real}\n\nAbstract type for integration on domains in D dimensions with value type T.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#HAdaptiveIntegration.EmbeddedCubature","page":"Docstrings","title":"HAdaptiveIntegration.EmbeddedCubature","text":"struct EmbeddedCubature{H,L,D,T}\n\nAn embedded cubature rule consisting of a high order cubature rule with H nodes and a low order cubature rule with L nodes. Note that the low order cubature uses nodes[1:L] as its nodes. The cubature nodes and weights are assume to be for the reference domain.\n\nFields:\n\nnodes::SVector{H,SVector{D,T}}: the cubature nodes.\nweights_high::SVector{H,T}: the cubature weights for the high order cubature.\nweights_low::SVector{L,T}: the cubature weights for the low order cubature.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#HAdaptiveIntegration.EmbeddedCubature-Union{Tuple{T}, Tuple{D}, Tuple{L}, Tuple{H}, Tuple{Any, HAdaptiveIntegration.Domain{D, T}}, Tuple{Any, HAdaptiveIntegration.Domain{D, T}, Any}} where {H, L, D, T}","page":"Docstrings","title":"HAdaptiveIntegration.EmbeddedCubature","text":"(ec::EmbeddedCubature{H,L,D,T})(\n    fct, domain::Domain{D,T}, norm=x -> LinearAlgebra.norm(x, Inf)\n) where {H,L,D,T}\n\nReturn I_high and norm(I_high - I_low) where I_high and I_low are the result of the high order cubature and the low order cubature on domain. The function fct must take a SVector{D,T} to a return type K, and K must support the multiplication by a scalar and the addition. Note that there is no check, beyond compatibility of dimension and type, that the embedded cubature is for the right domain.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.GrundmannMoeller","page":"Docstrings","title":"HAdaptiveIntegration.GrundmannMoeller","text":"struct GrundmannMoeller\n\nCubature rule for a dim-simplex of degree deg.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#HAdaptiveIntegration.Orthotope","page":"Docstrings","title":"HAdaptiveIntegration.Orthotope","text":"struct Orthotope{D,T} <: Domain{D,T}\n\nAxes-aligned Orthotope in D dimensions given by two points low_corner and high_corner. Note that, we must have low_corner .≤ high_corner.\n\nFields:\n\nlow_corner::SVector{D,T}: the low corner.\nhigh_corner::SVector{D,T}: the high corner.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#HAdaptiveIntegration.Rectangle","page":"Docstrings","title":"HAdaptiveIntegration.Rectangle","text":"Rectangle{T} = Orthotope{2,T}\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#HAdaptiveIntegration.Segment","page":"Docstrings","title":"HAdaptiveIntegration.Segment","text":"Segment{T} = Orthotope{1,T}\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#HAdaptiveIntegration.Simplex","page":"Docstrings","title":"HAdaptiveIntegration.Simplex","text":"struct Simplex{D,T,N} <: Domain{D,T}\n\nA simplex in D dimensions with N=D+1 vertices of value type T.\n\nFields:\n\nvertices::SVector{N,SVector{D,T}}: vertices of the simplex.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#HAdaptiveIntegration.TabulatedEmbeddedCubature","page":"Docstrings","title":"HAdaptiveIntegration.TabulatedEmbeddedCubature","text":"struct TabulatedEmbeddedCubature\n\nAn embedded cubature rule consisting of a high order cubature rule and a low order cubature rule. Note that the low order cubature uses nodes[1:L] as its nodes where L is the length of the weights_low. The cubature nodes and weights are assume to be for the reference domain.\n\nFields:\n\ndescription::String: description of the embedded cubature.\ndomain::String: domain of the cubature, the possible values are \"reference segment\", \"reference rectangle\", \"reference cuboid\", \"reference triangle\", and \"reference tetrahedron\".\nreference::String: where the values are found.\nnb_significant_digits::Int: number of significant digits on the node and weight values, 10^-nb_significant_digits give the relative precision of the values.\nnodes::Vector{Vector{String}}: the cubature nodes.\nweights_high::Vector{String}: the cubature weights for the high order cubature.\norder_high::Int: order of the high order cubature.\nweights_low::Vector{String}: the cubature weights for the low order cubature.\norder_low::Int: order of the low order cubature.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#HAdaptiveIntegration.Tetrahedron","page":"Docstrings","title":"HAdaptiveIntegration.Tetrahedron","text":"Tetrahedron{T} = Simplex{3,T,4}\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#HAdaptiveIntegration.Triangle","page":"Docstrings","title":"HAdaptiveIntegration.Triangle","text":"Triangle{T} = Simplex{2,T,3}\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#HAdaptiveIntegration.abs_det_jac-Union{Tuple{HAdaptiveIntegration.Orthotope{D, T}}, Tuple{T}, Tuple{D}} where {D, T}","page":"Docstrings","title":"HAdaptiveIntegration.abs_det_jac","text":"abs_det_jac(domain::Domain)\n\nReturn the absolute value of the Jacobian's determinant of the map from the reference domain to the physical domain domain.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.allocate_buffer-Union{Tuple{DOM}, Tuple{Any, DOM}, Tuple{Any, DOM, HAdaptiveIntegration.EmbeddedCubature}} where DOM<:HAdaptiveIntegration.Domain","page":"Docstrings","title":"HAdaptiveIntegration.allocate_buffer","text":"allocate_buffer(\n    fct, domain::DOM, ec::EmbeddedCubature=default_embedded_cubature(domain)\n) where {DOM<:Domain}\n\nTBW\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.check_order-Union{Tuple{T}, Tuple{D}, Tuple{HAdaptiveIntegration.TabulatedEmbeddedCubature, HAdaptiveIntegration.Domain{D, T}}} where {D, T}","page":"Docstrings","title":"HAdaptiveIntegration.check_order","text":"check_order(\n    tec::TabulatedEmbeddedCubature,\n    domain::Domain{D,T};\n    atol=zero(T),\n    rtol=(atol > zero(T)) ? zero(T) : 10 * eps(T),\n) where {D,T}\n\nReturn 0 if the tabulated embedded cubature on the reference domain of domain integrate exactly (within tolerance) the monomials up to degree order_high for the high oder cubature and order_low for the low order cubature. Else return 1.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.check_order-Union{Tuple{T}, Tuple{D}, Tuple{L}, Tuple{H}, Tuple{HAdaptiveIntegration.EmbeddedCubature{H, L, D, T}, HAdaptiveIntegration.Domain{D, T}, Int64, Int64}} where {H, L, D, T}","page":"Docstrings","title":"HAdaptiveIntegration.check_order","text":"check_order(\n    ec::EmbeddedCubature{H,L,D,T},\n    domain::Domain{D,T},\n    order_high::Int,\n    order_low::Int;\n    atol=zero(T),\n    rtol=(atol > zero(T)) ? zero(T) : 10 * eps(T),\n) where {H,L,D,T}\n\nReturn 0 if the embedded cubature on the reference domain of domain integrate exactly (within tolerance) the monomials up to degree order_high for the high oder cubature and order_low for the low order cubature. Else return 1.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.check_subdivision-Union{Tuple{T}, Tuple{D}, Tuple{Any, HAdaptiveIntegration.Domain{D, T}}} where {D, T}","page":"Docstrings","title":"HAdaptiveIntegration.check_subdivision","text":"check_subdivision(\n    subdiv_algo,\n    domain::Domain{D,T};\n    atol=zero(T),\n    rtol=(atol > zero(T)) ? zero(T) : 10 * eps(T),\n) where {D,T}\n\nReturn 0 if the sum of the volume of the subdomain by the subdiv_algo is equal to the volume of the domain else return 1.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.combinations-Tuple{Int64, Int64}","page":"Docstrings","title":"HAdaptiveIntegration.combinations","text":"combinations(n::Int, k::Int)\n\nHelper function to generate all combinations of k elements from 1:n, similar to calling combinations(1:n, k) from Combinatorics.jl.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.cuboid-Tuple{Any, Any}","page":"Docstrings","title":"HAdaptiveIntegration.cuboid","text":"cuboid(low_corner, high_corner)\n\nReturn an axes-aligned cuboid given by two 3d-points low_corner and high_corner. Note that, we must have low_corner .≤ high_corner.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.default_embedded_cubature-Union{Tuple{HAdaptiveIntegration.Orthotope{1, T}}, Tuple{T}} where T","page":"Docstrings","title":"HAdaptiveIntegration.default_embedded_cubature","text":"default_embedded_cubature(domain::Domain)\n\nReturn a default embedded cubature for the domains:\n\ndimension 1:\nsegment: SEGMENT_GK15\ndimension 2:\nrectangle: SQUARE_CHG25\ntriangle: TRIANGLE_RL19\ndimension 3:\ncuboid: CUBE_BE65\ntetrahedron: TETRAHEDRON_GM35\ndimension d:\nsimplex: [GrundmannMoeller](@ref)(d, 7)\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.default_subdivision-Tuple{HAdaptiveIntegration.Orthotope{1}}","page":"Docstrings","title":"HAdaptiveIntegration.default_subdivision","text":"default_subdivision(domain::Domain)\n\nReturn the default algorithm to subdivide domain.\n\ndimension 1:\nsegment: subdivide_segment2\ndimension 2:\nrectangle: subdivide_rectangle4\ntriangle: subdivide_triangle4\ndimension 3:\ncuboid: subdivide_cuboid8\ntetrahedron: subdivide_tetrahedron8\ndimension d:\nsimplex: subdivide_simplex\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.embedded_cubature","page":"Docstrings","title":"HAdaptiveIntegration.embedded_cubature","text":"embedded_cubature(gm::GrundmannMoeller, T::DataType=Float64)\n\nReturn the embedded cubature of Grundmann and Möller.\n\n\n\n\n\n","category":"function"},{"location":"docstrings/#HAdaptiveIntegration.embedded_cubature-2","page":"Docstrings","title":"HAdaptiveIntegration.embedded_cubature","text":"embedded_cubature(tec::TabulatedEmbeddedCubature, T::DataType=Float64)\n\nReturn the embedded cubature with value type T from a TabulatedEmbeddedCubature.\n\n\n\n\n\n","category":"function"},{"location":"docstrings/#HAdaptiveIntegration.embedded_cubature-Union{Tuple{VT}, Tuple{Vector{VT}, VT, VT}} where VT<:(AbstractVector{<:Real})","page":"Docstrings","title":"HAdaptiveIntegration.embedded_cubature","text":"embedded_cubature(\n    nodes::Vector{VT}, weights_high::VT, weights_low::VT\n) where {VT<:AbstractVector{<:Real}}\n\nReturn an embedded cubature form a vector of nodes and two vector of weights for the high order and low order cubature.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.integral_monomial-Union{Tuple{T}, Tuple{D}, Tuple{HAdaptiveIntegration.Orthotope{D, T}, Int64}} where {D, T}","page":"Docstrings","title":"HAdaptiveIntegration.integral_monomial","text":"integral_monomial(domain::Domain{D,T}, tot_deg_max::Int) where {D,T}\n\nReturn the values of monomial's integral over the reference domain. It return a Vector{ Vector{ Pair{ Ntuple{dim,Int}, Rational{Int} } } }. The outer vector is index by the total degree + 1, for the total degree form 0 to tot_deg_max. The inner vector contain Pair{ Ntuple{dim,Int}, Rational{Int} } where the Ntuple{D,Int} is the multi-index of the monomial and Rational{Int} is the value of the integral.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.integrate-Union{Tuple{T}, Tuple{D}, Tuple{L}, Tuple{H}, Tuple{Any, HAdaptiveIntegration.Domain{D, T}}} where {H, L, D, T}","page":"Docstrings","title":"HAdaptiveIntegration.integrate","text":"integrate(\n    fct,\n    domain::Domain{D,T};\n    embedded_cubature::EmbeddedCubature{H,L,D,T}=default_embedded_cubature(domain),\n    subdiv_algo=default_subdivision(domain),\n    buffer=nothing,\n    norm=x -> LinearAlgebra.norm(x, Inf),\n    atol=zero(T),\n    rtol=(atol > zero(T)) ? zero(T) : sqrt(eps(T)),\n    maxsubdiv=512 * 2^D,\n) where {H,L,D,T}\n\nReturn I and E where I is the integral of the function fct over domain and E is an error estimate.\n\nArguments\n\nfct: a function that must take a SVector{D,T} to a return type K, with K must  support the multiplication by a scalar of type T and the addition.\ndomain::Domain{D,T}: the integration domain. Currently, we support segment, triangle, rectangle, tetrahedron, cuboid, and d-dimensional simplex.\n\nOptional arguments\n\nec::EmbeddedCubature{H,L,D,T}=default_embedded_cubature(domain): the embedded cubature,  each supported domain has a default_embedded_cubature.\nsubdiv_algo=default_subdivision(domain): the subdivision algorithm, each domain has a default_subdivision.\nbuffer=nothing: heap use to do the adaptive algorithm, can be allocated using  allocate_buffer, which might result in performance gain if multiple call to  integrate is perform.\nnorm=x -> LinearAlgebra.norm(x, Inf): norm used to estimate the error.\natol=zero(T): absolute tolerance.\nrtol=(atol > zero(T)) ? zero(T) : sqrt(eps(T)): relative tolerance.\nmaxsubdiv=512 * 2^D: maximum number of subdivision.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.map_from_reference-Union{Tuple{HAdaptiveIntegration.Orthotope{D, T}}, Tuple{T}, Tuple{D}} where {D, T}","page":"Docstrings","title":"HAdaptiveIntegration.map_from_reference","text":"map_from_reference(domain::Domain)\n\nReturn an anonymous function that maps the reference domain to the physical domain domain.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.map_to_reference-Union{Tuple{HAdaptiveIntegration.Orthotope{D, T}}, Tuple{T}, Tuple{D}} where {D, T}","page":"Docstrings","title":"HAdaptiveIntegration.map_to_reference","text":"map_to_reference(domain::Domain)\n\nReturn an anonymous function that maps the physical domain domain to the reference domain.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.orthotope-Tuple{Any, Any}","page":"Docstrings","title":"HAdaptiveIntegration.orthotope","text":"orthotope(low_corner, high_corner)\n\nReturn an axes-aligned orthotope in D dimensions given by two points low_corner and high_corner. Note that, we must have low_corner .≤ high_corner.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.rectangle-Tuple{Any, Any}","page":"Docstrings","title":"HAdaptiveIntegration.rectangle","text":"rectangle(low_corner, high_corner)\n\nReturn an axes-aligned rectangle given by two 2d-points low_corner and high_corner. Note that, we must have low_corner .≤ high_corner.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.reference_orthotope","page":"Docstrings","title":"HAdaptiveIntegration.reference_orthotope","text":"reference_orthotope(D::Int, T::DataType=Float64)\n\nReturn the reference D-dimensional orthotope [0, 1]ᴰ with value type T.\n\n\n\n\n\n","category":"function"},{"location":"docstrings/#HAdaptiveIntegration.reference_simplex","page":"Docstrings","title":"HAdaptiveIntegration.reference_simplex","text":"reference_simplex(D::Int, T::DataType=Float64)\n\nReturn the reference D-dimensional simplex with value type T, which is the convex hull of the N=D+1 points (0,...,0), (1,0,...,0), (0,1,0,...,0), ..., (0,...,0,1).\n\n\n\n\n\n","category":"function"},{"location":"docstrings/#HAdaptiveIntegration.segment-Union{Tuple{T}, Tuple{T, T}} where T<:Real","page":"Docstrings","title":"HAdaptiveIntegration.segment","text":"segment(xmin, xmax)\n\nReturn a segment in 1 dimensions representing the interval [xmin, xmax] with xmin ≤ xmax.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.simplex-Tuple","page":"Docstrings","title":"HAdaptiveIntegration.simplex","text":"simplex(vertices...)\n\nReturn a D-simplex from a collection of vertices. Note that all points must have the same length D and there must be N=D+1 points.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.subdivide_cuboid8-Union{Tuple{HAdaptiveIntegration.Orthotope{3, T}}, Tuple{T}} where T","page":"Docstrings","title":"HAdaptiveIntegration.subdivide_cuboid8","text":"subdivide_cuboid8(c::Cuboid)\n\nDivide the cuboid c into 8 cuboid by connecting the center of the cuboid to the midpoints of the edges.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.subdivide_rectangle4-Union{Tuple{HAdaptiveIntegration.Orthotope{2, T}}, Tuple{T}} where T","page":"Docstrings","title":"HAdaptiveIntegration.subdivide_rectangle4","text":"subdivide_rectangle4(r::Rectangle)\n\nDivide the rectangle r into four squares by connecting the center of the square to the midpoints of the edges.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.subdivide_reference_simplex-Union{Tuple{Val{D}}, Tuple{T}, Tuple{D}, Tuple{Val{D}, Type{T}}} where {D, T}","page":"Docstrings","title":"HAdaptiveIntegration.subdivide_reference_simplex","text":"subdivide_reference_simplex(::Val{D}, ::Type{T}=Float64) where {D,T}\n\nLike subdivide_simplex, but operates on the reference simplex. Since the output depends only on the dimension D, and the type T used to represent coordinates, this function is generated for each combination of D and T.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.subdivide_segment2-Union{Tuple{HAdaptiveIntegration.Orthotope{1, T}}, Tuple{T}} where T","page":"Docstrings","title":"HAdaptiveIntegration.subdivide_segment2","text":"subdivide_segment2(s::Segment)\n\nDivide the segment s into two segments of equal length.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.subdivide_simplex-Union{Tuple{HAdaptiveIntegration.Simplex{D, T, N}}, Tuple{N}, Tuple{T}, Tuple{D}} where {D, T, N}","page":"Docstrings","title":"HAdaptiveIntegration.subdivide_simplex","text":"subdivide_simplex(s::Simplex)\n\nSubdivide a D-simplex into 2ᴰ simplices by using the Freudenthal triangulation.\n\nImplements the RedRefinementND algorithm in Simplicial grid refinement: on Freudenthal's algorithm and the optimal number of congruence classes.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.subdivide_tetrahedron8-Union{Tuple{HAdaptiveIntegration.Simplex{3, T, 4}}, Tuple{T}} where T","page":"Docstrings","title":"HAdaptiveIntegration.subdivide_tetrahedron8","text":"subdivide_tetrahedron8(t::Tetrahedron)\n\nDivide the tetrahedron t into eight tetrahedra by connecting the midpoints of the edges.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.subdivide_triangle2-Union{Tuple{HAdaptiveIntegration.Simplex{2, T, 3}}, Tuple{T}} where T","page":"Docstrings","title":"HAdaptiveIntegration.subdivide_triangle2","text":"subdivide_triangle2(s::Triangle)\n\nDivide the triangle t into two triangles by connecting the first point of t to the midpoints of the two other points.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.subdivide_triangle4-Union{Tuple{HAdaptiveIntegration.Simplex{2, T, 3}}, Tuple{T}} where T","page":"Docstrings","title":"HAdaptiveIntegration.subdivide_triangle4","text":"subdivide_triangle4(t::Triangle)\n\nDivide the triangle t into four triangles by connecting the midpoints of the edges.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.tetrahedron-NTuple{4, Any}","page":"Docstrings","title":"HAdaptiveIntegration.tetrahedron","text":"tetrahedron(a, b, c, d)\n\nReturn a tetrahedron in 3 dimensions given by four 3d-points a, b, c, and d.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#HAdaptiveIntegration.triangle-Tuple{Any, Any, Any}","page":"Docstrings","title":"HAdaptiveIntegration.triangle","text":"triangle(a, b, c)\n\nReturn a triangle in 2 dimensions given by three 2d-points a, b, and c.\n\n\n\n\n\n","category":"method"}]
}
