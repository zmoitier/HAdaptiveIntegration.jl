"""
Radon-Laurie with 19 nodes.
Radon-Laurie with 19 nodes.
"""
const TRIANGLE_RL19 = TabulatedEmbeddedCubature{Triangle}(;
    description="Radon-Laurie with 19 nodes (TRIANGLE_RL19)",
    reference="""https://www.math.unipd.it/~alvise/SETS_CUBATURE_TRIANGLE/radon/set_radon_standard.m
    https://www.math.unipd.it/~alvise/SETS_CUBATURE_TRIANGLE/laurie/set_laurie_standard.m""",
    precision=15,
    nodes=[
        ["3.333333333333333e-01", "3.333333333333333e-01"],
        ["7.974269853530872e-01", "1.012865073234563e-01"],
        ["1.012865073234563e-01", "7.974269853530872e-01"],
        ["1.012865073234563e-01", "1.012865073234563e-01"],
        ["5.971587178976980e-02", "4.701420641051151e-01"],
        ["4.701420641051151e-01", "5.971587178976980e-02"],
        ["4.701420641051151e-01", "4.701420641051151e-01"],
        ["5.357953464498992e-01", "2.321023267750504e-01"],
        ["2.321023267750504e-01", "5.357953464498992e-01"],
        ["2.321023267750504e-01", "2.321023267750504e-01"],
        ["9.410382782311209e-01", "2.948086088443960e-02"],
        ["2.948086088443960e-02", "9.410382782311209e-01"],
        ["2.948086088443960e-02", "2.948086088443960e-02"],
        ["7.384168123405100e-01", "2.321023267750504e-01"],
        ["7.384168123405100e-01", "2.948086088443960e-02"],
        ["2.321023267750504e-01", "7.384168123405100e-01"],
        ["2.321023267750504e-01", "2.948086088443960e-02"],
        ["2.948086088443960e-02", "7.384168123405100e-01"],
        ["2.948086088443960e-02", "2.321023267750504e-01"],
    ],
    weights_high=[
        "1.893054560015735e-02",
        "1.881021270659145e-02",
        "1.881021270659145e-02",
        "1.881021270659145e-02",
        "3.917867612205870e-02",
        "3.917867612205870e-02",
        "3.917867612205870e-02",
        "5.813573982848295e-02",
        "5.813573982848295e-02",
        "5.813573982848295e-02",
        "6.722133687582750e-03",
        "6.722133687582750e-03",
        "6.722133687582750e-03",
        "1.875486122761585e-02",
        "1.875486122761585e-02",
        "1.875486122761585e-02",
        "1.875486122761585e-02",
        "1.875486122761585e-02",
        "1.875486122761585e-02",
    ],
    order_high=8,
    weights_low=[
        "1.125000000000000e-01",
        "6.296959027241358e-02",
        "6.296959027241358e-02",
        "6.296959027241358e-02",
        "6.619707639425308e-02",
        "6.619707639425308e-02",
        "6.619707639425308e-02",
    ],
    order_low=5,
)

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
