"""
Grundmann-Möller with 19 nodes.
"""
const TRIANGLE_GM19 = TabulatedEmbeddedCubature{Triangle}(;
    description="Grundmann-Möller with 19 nodes (TRIANGLE_GM19)",
    reference="https://epubs.siam.org/doi/10.1137/0715019",
    precision=35,
    nodes=[
        [
            "3.333333333333333333333333333333333333e-01",
            "3.333333333333333333333333333333333333e-01",
        ],
        [
            "2.000000000000000000000000000000000000e-01",
            "6.000000000000000000000000000000000000e-01",
        ],
        [
            "6.000000000000000000000000000000000000e-01",
            "2.000000000000000000000000000000000000e-01",
        ],
        [
            "2.000000000000000000000000000000000000e-01",
            "2.000000000000000000000000000000000000e-01",
        ],
        [
            "1.428571428571428571428571428571428571e-01",
            "7.142857142857142857142857142857142857e-01",
        ],
        [
            "4.285714285714285714285714285714285714e-01",
            "4.285714285714285714285714285714285714e-01",
        ],
        [
            "1.428571428571428571428571428571428571e-01",
            "4.285714285714285714285714285714285714e-01",
        ],
        [
            "7.142857142857142857142857142857142857e-01",
            "1.428571428571428571428571428571428571e-01",
        ],
        [
            "4.285714285714285714285714285714285714e-01",
            "1.428571428571428571428571428571428571e-01",
        ],
        [
            "1.428571428571428571428571428571428571e-01",
            "1.428571428571428571428571428571428571e-01",
        ],
        [
            "1.111111111111111111111111111111111111e-01",
            "7.777777777777777777777777777777777778e-01",
        ],
        [
            "3.333333333333333333333333333333333333e-01",
            "5.555555555555555555555555555555555556e-01",
        ],
        [
            "1.111111111111111111111111111111111111e-01",
            "5.555555555555555555555555555555555556e-01",
        ],
        [
            "5.555555555555555555555555555555555556e-01",
            "3.333333333333333333333333333333333333e-01",
        ],
        [
            "1.111111111111111111111111111111111111e-01",
            "3.333333333333333333333333333333333333e-01",
        ],
        [
            "7.777777777777777777777777777777777778e-01",
            "1.111111111111111111111111111111111111e-01",
        ],
        [
            "5.555555555555555555555555555555555556e-01",
            "1.111111111111111111111111111111111111e-01",
        ],
        [
            "3.333333333333333333333333333333333333e-01",
            "1.111111111111111111111111111111111111e-01",
        ],
        [
            "1.111111111111111111111111111111111111e-01",
            "1.111111111111111111111111111111111111e-01",
        ],
    ],
    weights_high=[
        "1.980364118303571428571428571428571429e-01",
        "1.211015004960317460317460317460317460e-01",
        "1.211015004960317460317460317460317460e-01",
        "1.211015004960317460317460317460317460e-01",
        "-3.191433376736111111111111111111111111e-01",
        "-3.191433376736111111111111111111111111e-01",
        "-3.191433376736111111111111111111111111e-01",
        "-3.191433376736111111111111111111111111e-01",
        "-3.191433376736111111111111111111111111e-01",
        "-3.191433376736111111111111111111111111e-01",
        "2.059465680803571428571428571428571429e-01",
        "2.059465680803571428571428571428571429e-01",
        "2.059465680803571428571428571428571429e-01",
        "2.059465680803571428571428571428571429e-01",
        "2.059465680803571428571428571428571429e-01",
        "2.059465680803571428571428571428571429e-01",
        "2.059465680803571428571428571428571429e-01",
        "2.059465680803571428571428571428571429e-01",
        "2.059465680803571428571428571428571429e-01",
    ],
    order_high=7,
    weights_low=[
        "6.328125000000000000000000000000000000e-02",
        "-2.712673611111111111111111111111111111e-01",
        "-2.712673611111111111111111111111111111e-01",
        "-2.712673611111111111111111111111111111e-01",
        "2.084201388888888888888888888888888889e-01",
        "2.084201388888888888888888888888888889e-01",
        "2.084201388888888888888888888888888889e-01",
        "2.084201388888888888888888888888888889e-01",
        "2.084201388888888888888888888888888889e-01",
        "2.084201388888888888888888888888888889e-01",
    ],
    order_low=5,
)

"""
Radon-Laurie with 19 nodes.
"""
const TRIANGLE_RL19 = TabulatedEmbeddedCubature{Triangle}(;
    description="Radon-Laurie with 19 nodes (TRIANGLE_RL19)",
    reference="""https://doi.org/10.1145/355993.356001""",
    precision=35,
    nodes=[
        [
            "3.333333333333333333333333333333333333e-01",
            "3.333333333333333333333333333333333333e-01",
        ],
        [
            "7.974269853530873223980252761697523439e-01",
            "1.012865073234563388009873619151238281e-01",
        ],
        [
            "1.012865073234563388009873619151238281e-01",
            "7.974269853530873223980252761697523439e-01",
        ],
        [
            "1.012865073234563388009873619151238281e-01",
            "1.012865073234563388009873619151238281e-01",
        ],
        [
            "5.971587178976982045911758097310479897e-02",
            "4.701420641051150897704412095134476005e-01",
        ],
        [
            "4.701420641051150897704412095134476005e-01",
            "5.971587178976982045911758097310479897e-02",
        ],
        [
            "4.701420641051150897704412095134476005e-01",
            "4.701420641051150897704412095134476005e-01",
        ],
        [
            "9.410382782311208665596303797019388487e-01",
            "2.948086088443956672018481014903057564e-02",
        ],
        [
            "2.948086088443956672018481014903057564e-02",
            "9.410382782311208665596303797019388487e-01",
        ],
        [
            "2.948086088443956672018481014903057564e-02",
            "2.948086088443956672018481014903057564e-02",
        ],
        [
            "5.357953464498992646629508988845634681e-01",
            "2.321023267750503676685245505577182659e-01",
        ],
        [
            "2.321023267750503676685245505577182659e-01",
            "5.357953464498992646629508988845634681e-01",
        ],
        [
            "2.321023267750503676685245505577182659e-01",
            "2.321023267750503676685245505577182659e-01",
        ],
        [
            "2.948086088443956672018481014903057564e-02",
            "2.321023267750503676685245505577182659e-01",
        ],
        [
            "2.948086088443956672018481014903057564e-02",
            "7.384168123405100656112906392932511584e-01",
        ],
        [
            "2.321023267750503676685245505577182659e-01",
            "2.948086088443956672018481014903057564e-02",
        ],
        [
            "2.321023267750503676685245505577182659e-01",
            "7.384168123405100656112906392932511584e-01",
        ],
        [
            "7.384168123405100656112906392932511584e-01",
            "2.948086088443956672018481014903057564e-02",
        ],
        [
            "7.384168123405100656112906392932511584e-01",
            "2.321023267750503676685245505577182659e-01",
        ],
    ],
    weights_high=[
        "1.893054560015734165415411067800806149e-02",
        "1.881021270659148607215700260658383995e-02",
        "1.881021270659148607215700260658383995e-02",
        "1.881021270659148607215700260658383995e-02",
        "3.917867612205866877772300044712309003e-02",
        "3.917867612205866877772300044712309003e-02",
        "3.917867612205866877772300044712309003e-02",
        "6.722133687582700949055532526590556751e-03",
        "6.722133687582700949055532526590556751e-03",
        "6.722133687582700949055532526590556751e-03",
        "5.813573982848294819737435282749222628e-02",
        "5.813573982848294819737435282749222628e-02",
        "5.813573982848294819737435282749222628e-02",
        "1.875486122761587439281937068310379991e-02",
        "1.875486122761587439281937068310379991e-02",
        "1.875486122761587439281937068310379991e-02",
        "1.875486122761587439281937068310379991e-02",
        "1.875486122761587439281937068310379991e-02",
        "1.875486122761587439281937068310379991e-02",
    ],
    order_high=8,
    weights_low=[
        "1.125000000000000000000000000000000000e-01",
        "6.296959027241357629784197275009066683e-02",
        "6.296959027241357629784197275009066683e-02",
        "6.296959027241357629784197275009066683e-02",
        "6.619707639425309036882469391657599984e-02",
        "6.619707639425309036882469391657599984e-02",
        "6.619707639425309036882469391657599984e-02",
    ],
    order_low=5,
)

"""
Grundmann-Möller with 35 nodes.
"""
const TETRAHEDRON_GM35 = TabulatedEmbeddedCubature{Tetrahedron}(;
    description="Grundmann-Möller with 35 nodes (TETRAHEDRON_GM35)",
    reference="https://epubs.siam.org/doi/10.1137/0715019",
    precision=35,
    nodes=[
        [
            "2.500000000000000000000000000000000000e-01",
            "2.500000000000000000000000000000000000e-01",
            "2.500000000000000000000000000000000000e-01",
        ],
        [
            "1.666666666666666666666666666666666667e-01",
            "1.666666666666666666666666666666666667e-01",
            "5.000000000000000000000000000000000000e-01",
        ],
        [
            "1.666666666666666666666666666666666667e-01",
            "5.000000000000000000000000000000000000e-01",
            "1.666666666666666666666666666666666667e-01",
        ],
        [
            "5.000000000000000000000000000000000000e-01",
            "1.666666666666666666666666666666666667e-01",
            "1.666666666666666666666666666666666667e-01",
        ],
        [
            "1.666666666666666666666666666666666667e-01",
            "1.666666666666666666666666666666666667e-01",
            "1.666666666666666666666666666666666667e-01",
        ],
        [
            "1.250000000000000000000000000000000000e-01",
            "1.250000000000000000000000000000000000e-01",
            "6.250000000000000000000000000000000000e-01",
        ],
        [
            "1.250000000000000000000000000000000000e-01",
            "3.750000000000000000000000000000000000e-01",
            "3.750000000000000000000000000000000000e-01",
        ],
        [
            "3.750000000000000000000000000000000000e-01",
            "1.250000000000000000000000000000000000e-01",
            "3.750000000000000000000000000000000000e-01",
        ],
        [
            "1.250000000000000000000000000000000000e-01",
            "1.250000000000000000000000000000000000e-01",
            "3.750000000000000000000000000000000000e-01",
        ],
        [
            "1.250000000000000000000000000000000000e-01",
            "6.250000000000000000000000000000000000e-01",
            "1.250000000000000000000000000000000000e-01",
        ],
        [
            "3.750000000000000000000000000000000000e-01",
            "3.750000000000000000000000000000000000e-01",
            "1.250000000000000000000000000000000000e-01",
        ],
        [
            "1.250000000000000000000000000000000000e-01",
            "3.750000000000000000000000000000000000e-01",
            "1.250000000000000000000000000000000000e-01",
        ],
        [
            "6.250000000000000000000000000000000000e-01",
            "1.250000000000000000000000000000000000e-01",
            "1.250000000000000000000000000000000000e-01",
        ],
        [
            "3.750000000000000000000000000000000000e-01",
            "1.250000000000000000000000000000000000e-01",
            "1.250000000000000000000000000000000000e-01",
        ],
        [
            "1.250000000000000000000000000000000000e-01",
            "1.250000000000000000000000000000000000e-01",
            "1.250000000000000000000000000000000000e-01",
        ],
        [
            "1.000000000000000000000000000000000000e-01",
            "1.000000000000000000000000000000000000e-01",
            "7.000000000000000000000000000000000000e-01",
        ],
        [
            "1.000000000000000000000000000000000000e-01",
            "3.000000000000000000000000000000000000e-01",
            "5.000000000000000000000000000000000000e-01",
        ],
        [
            "3.000000000000000000000000000000000000e-01",
            "1.000000000000000000000000000000000000e-01",
            "5.000000000000000000000000000000000000e-01",
        ],
        [
            "1.000000000000000000000000000000000000e-01",
            "1.000000000000000000000000000000000000e-01",
            "5.000000000000000000000000000000000000e-01",
        ],
        [
            "1.000000000000000000000000000000000000e-01",
            "5.000000000000000000000000000000000000e-01",
            "3.000000000000000000000000000000000000e-01",
        ],
        [
            "3.000000000000000000000000000000000000e-01",
            "3.000000000000000000000000000000000000e-01",
            "3.000000000000000000000000000000000000e-01",
        ],
        [
            "1.000000000000000000000000000000000000e-01",
            "3.000000000000000000000000000000000000e-01",
            "3.000000000000000000000000000000000000e-01",
        ],
        [
            "5.000000000000000000000000000000000000e-01",
            "1.000000000000000000000000000000000000e-01",
            "3.000000000000000000000000000000000000e-01",
        ],
        [
            "3.000000000000000000000000000000000000e-01",
            "1.000000000000000000000000000000000000e-01",
            "3.000000000000000000000000000000000000e-01",
        ],
        [
            "1.000000000000000000000000000000000000e-01",
            "1.000000000000000000000000000000000000e-01",
            "3.000000000000000000000000000000000000e-01",
        ],
        [
            "1.000000000000000000000000000000000000e-01",
            "7.000000000000000000000000000000000000e-01",
            "1.000000000000000000000000000000000000e-01",
        ],
        [
            "3.000000000000000000000000000000000000e-01",
            "5.000000000000000000000000000000000000e-01",
            "1.000000000000000000000000000000000000e-01",
        ],
        [
            "1.000000000000000000000000000000000000e-01",
            "5.000000000000000000000000000000000000e-01",
            "1.000000000000000000000000000000000000e-01",
        ],
        [
            "5.000000000000000000000000000000000000e-01",
            "3.000000000000000000000000000000000000e-01",
            "1.000000000000000000000000000000000000e-01",
        ],
        [
            "3.000000000000000000000000000000000000e-01",
            "3.000000000000000000000000000000000000e-01",
            "1.000000000000000000000000000000000000e-01",
        ],
        [
            "1.000000000000000000000000000000000000e-01",
            "3.000000000000000000000000000000000000e-01",
            "1.000000000000000000000000000000000000e-01",
        ],
        [
            "7.000000000000000000000000000000000000e-01",
            "1.000000000000000000000000000000000000e-01",
            "1.000000000000000000000000000000000000e-01",
        ],
        [
            "5.000000000000000000000000000000000000e-01",
            "1.000000000000000000000000000000000000e-01",
            "1.000000000000000000000000000000000000e-01",
        ],
        [
            "3.000000000000000000000000000000000000e-01",
            "1.000000000000000000000000000000000000e-01",
            "1.000000000000000000000000000000000000e-01",
        ],
        [
            "1.000000000000000000000000000000000000e-01",
            "1.000000000000000000000000000000000000e-01",
            "1.000000000000000000000000000000000000e-01",
        ],
    ],
    weights_high=[
        "-8.465608465608465608465608465608465608e-03",
        "5.424107142857142857142857142857142857e-02",
        "5.424107142857142857142857142857142857e-02",
        "5.424107142857142857142857142857142857e-02",
        "5.424107142857142857142857142857142857e-02",
        "-9.029982363315696649029982363315696649e-02",
        "-9.029982363315696649029982363315696649e-02",
        "-9.029982363315696649029982363315696649e-02",
        "-9.029982363315696649029982363315696649e-02",
        "-9.029982363315696649029982363315696649e-02",
        "-9.029982363315696649029982363315696649e-02",
        "-9.029982363315696649029982363315696649e-02",
        "-9.029982363315696649029982363315696649e-02",
        "-9.029982363315696649029982363315696649e-02",
        "-9.029982363315696649029982363315696649e-02",
        "4.305831128747795414462081128747795414e-02",
        "4.305831128747795414462081128747795414e-02",
        "4.305831128747795414462081128747795414e-02",
        "4.305831128747795414462081128747795414e-02",
        "4.305831128747795414462081128747795414e-02",
        "4.305831128747795414462081128747795414e-02",
        "4.305831128747795414462081128747795414e-02",
        "4.305831128747795414462081128747795414e-02",
        "4.305831128747795414462081128747795414e-02",
        "4.305831128747795414462081128747795414e-02",
        "4.305831128747795414462081128747795414e-02",
        "4.305831128747795414462081128747795414e-02",
        "4.305831128747795414462081128747795414e-02",
        "4.305831128747795414462081128747795414e-02",
        "4.305831128747795414462081128747795414e-02",
        "4.305831128747795414462081128747795414e-02",
        "4.305831128747795414462081128747795414e-02",
        "4.305831128747795414462081128747795414e-02",
        "4.305831128747795414462081128747795414e-02",
        "4.305831128747795414462081128747795414e-02",
    ],
    order_high=7,
    weights_low=[
        "4.444444444444444444444444444444444444e-02",
        "-9.642857142857142857142857142857142857e-02",
        "-9.642857142857142857142857142857142857e-02",
        "-9.642857142857142857142857142857142857e-02",
        "-9.642857142857142857142857142857142857e-02",
        "5.079365079365079365079365079365079365e-02",
        "5.079365079365079365079365079365079365e-02",
        "5.079365079365079365079365079365079365e-02",
        "5.079365079365079365079365079365079365e-02",
        "5.079365079365079365079365079365079365e-02",
        "5.079365079365079365079365079365079365e-02",
        "5.079365079365079365079365079365079365e-02",
        "5.079365079365079365079365079365079365e-02",
        "5.079365079365079365079365079365079365e-02",
        "5.079365079365079365079365079365079365e-02",
    ],
    order_low=5,
)
