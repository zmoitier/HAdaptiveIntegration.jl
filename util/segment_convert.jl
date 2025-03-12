using StaticArrays
import Printf: Format, format

import HAdaptiveIntegration as hai

setprecision(BigFloat, 40; base=10)

segment = hai.segment(big"-1.0", big"1.0")

# https://www.advanpix.com/2011/11/07/gauss-kronrod-quadrature-nodes-weights/

# 7 points Gauss rule of order 13
g7 = (
    nodes=[
        big" 0.000000000000000000000000000000000e+00",
        big"-4.058451513773971669066064120769615e-01",
        big" 4.058451513773971669066064120769615e-01",
        big"-7.415311855993944398638647732807884e-01",
        big" 7.415311855993944398638647732807884e-01",
        big"-9.491079123427585245261896840478513e-01",
        big" 9.491079123427585245261896840478513e-01",
    ],
    weights=[
        big"4.179591836734693877551020408163265e-01",
        big"3.818300505051189449503697754889751e-01",
        big"3.818300505051189449503697754889751e-01",
        big"2.797053914892766679014677714237796e-01",
        big"2.797053914892766679014677714237796e-01",
        big"1.294849661688696932706114326790820e-01",
        big"1.294849661688696932706114326790820e-01",
    ],
)

# 15 points Kronrod rule of order 23
k15 = (
    nodes=[
        big" 0.000000000000000000000000000000000e+00",
        big"-4.058451513773971669066064120769615e-01",
        big" 4.058451513773971669066064120769615e-01",
        big"-7.415311855993944398638647732807884e-01",
        big" 7.415311855993944398638647732807884e-01",
        big"-9.491079123427585245261896840478513e-01",
        big" 9.491079123427585245261896840478513e-01",
        big"-2.077849550078984676006894037732449e-01",
        big" 2.077849550078984676006894037732449e-01",
        big"-5.860872354676911302941448382587296e-01",
        big" 5.860872354676911302941448382587296e-01",
        big"-8.648644233597690727897127886409262e-01",
        big" 8.648644233597690727897127886409262e-01",
        big"-9.914553711208126392068546975263285e-01",
        big" 9.914553711208126392068546975263285e-01",
    ],
    weights=[
        big"2.094821410847278280129991748917143e-01",
        big"1.903505780647854099132564024210137e-01",
        big"1.903505780647854099132564024210137e-01",
        big"1.406532597155259187451895905102379e-01",
        big"1.406532597155259187451895905102379e-01",
        big"6.309209262997855329070066318920429e-02",
        big"6.309209262997855329070066318920429e-02",
        big"2.044329400752988924141619992346491e-01",
        big"2.044329400752988924141619992346491e-01",
        big"1.690047266392679028265834265985503e-01",
        big"1.690047266392679028265834265985503e-01",
        big"1.047900103222501838398763225415180e-01",
        big"1.047900103222501838398763225415180e-01",
        big"2.293532201052922496373200805896959e-02",
        big"2.293532201052922496373200805896959e-02",
    ],
)

# 15 points Gauss rule of order 29
g15 = (
    nodes=[
        big" 0.000000000000000000000000000000000e+00",
        big"-2.011940939974345223006283033945962e-01",
        big" 2.011940939974345223006283033945962e-01",
        big"-3.941513470775633698972073709810455e-01",
        big" 3.941513470775633698972073709810455e-01",
        big"-5.709721726085388475372267372539106e-01",
        big" 5.709721726085388475372267372539106e-01",
        big"-7.244177313601700474161860546139380e-01",
        big" 7.244177313601700474161860546139380e-01",
        big"-8.482065834104272162006483207742169e-01",
        big" 8.482065834104272162006483207742169e-01",
        big"-9.372733924007059043077589477102095e-01",
        big" 9.372733924007059043077589477102095e-01",
        big"-9.879925180204854284895657185866126e-01",
        big" 9.879925180204854284895657185866126e-01",
    ],
    weights=[
        big"2.025782419255612728806201999675193e-01",
        big"1.984314853271115764561183264438393e-01",
        big"1.984314853271115764561183264438393e-01",
        big"1.861610000155622110268005618664228e-01",
        big"1.861610000155622110268005618664228e-01",
        big"1.662692058169939335532008604812088e-01",
        big"1.662692058169939335532008604812088e-01",
        big"1.395706779261543144478047945110283e-01",
        big"1.395706779261543144478047945110283e-01",
        big"1.071592204671719350118695466858693e-01",
        big"1.071592204671719350118695466858693e-01",
        big"7.036604748810812470926741645066734e-02",
        big"7.036604748810812470926741645066734e-02",
        big"3.075324199611726835462839357720442e-02",
        big"3.075324199611726835462839357720442e-02",
    ],
)

# 31 points Kronrod rule of order 47
k31 = (
    nodes=[
        big" 0.000000000000000000000000000000000e+00",
        big"-2.011940939974345223006283033945962e-01",
        big" 2.011940939974345223006283033945962e-01",
        big"-3.941513470775633698972073709810455e-01",
        big" 3.941513470775633698972073709810455e-01",
        big"-5.709721726085388475372267372539106e-01",
        big" 5.709721726085388475372267372539106e-01",
        big"-7.244177313601700474161860546139380e-01",
        big" 7.244177313601700474161860546139380e-01",
        big"-8.482065834104272162006483207742169e-01",
        big" 8.482065834104272162006483207742169e-01",
        big"-9.372733924007059043077589477102095e-01",
        big" 9.372733924007059043077589477102095e-01",
        big"-9.879925180204854284895657185866126e-01",
        big" 9.879925180204854284895657185866126e-01",
        big"-1.011420669187174990270742314473923e-01",
        big" 1.011420669187174990270742314473923e-01",
        big"-2.991800071531688121667800242663890e-01",
        big" 2.991800071531688121667800242663890e-01",
        big"-4.850818636402396806936557402323506e-01",
        big" 4.850818636402396806936557402323506e-01",
        big"-6.509967412974169705337358953132747e-01",
        big" 6.509967412974169705337358953132747e-01",
        big"-7.904185014424659329676492948179473e-01",
        big" 7.904185014424659329676492948179473e-01",
        big"-8.972645323440819008825096564544959e-01",
        big" 8.972645323440819008825096564544959e-01",
        big"-9.677390756791391342573479787843372e-01",
        big" 9.677390756791391342573479787843372e-01",
        big"-9.980022986933970602851728401522712e-01",
        big" 9.980022986933970602851728401522712e-01",
    ],
    weights=[
        big"1.013300070147915490173747927674925e-01",
        big"9.917359872179195933239317348460313e-02",
        big"9.917359872179195933239317348460313e-02",
        big"9.312659817082532122548687274734572e-02",
        big"9.312659817082532122548687274734572e-02",
        big"8.308050282313302103828924728610379e-02",
        big"8.308050282313302103828924728610379e-02",
        big"6.985412131872825870952007709914748e-02",
        big"6.985412131872825870952007709914748e-02",
        big"5.348152469092808726534314723943030e-02",
        big"5.348152469092808726534314723943030e-02",
        big"3.534636079137584622203794847836005e-02",
        big"3.534636079137584622203794847836005e-02",
        big"1.500794732931612253837476307580727e-02",
        big"1.500794732931612253837476307580727e-02",
        big"1.007698455238755950449466626175697e-01",
        big"1.007698455238755950449466626175697e-01",
        big"9.664272698362367850517990762758934e-02",
        big"9.664272698362367850517990762758934e-02",
        big"8.856444305621177064727544369377430e-02",
        big"8.856444305621177064727544369377430e-02",
        big"7.684968075772037889443277748265901e-02",
        big"7.684968075772037889443277748265901e-02",
        big"6.200956780067064028513923096080293e-02",
        big"6.200956780067064028513923096080293e-02",
        big"4.458975132476487660822729937327969e-02",
        big"4.458975132476487660822729937327969e-02",
        big"2.546084732671532018687400101965336e-02",
        big"2.546084732671532018687400101965336e-02",
        big"5.377479872923348987792051430127650e-03",
        big"5.377479872923348987792051430127650e-03",
    ],
)

for (name, cbt, n) in [("G7", g7, 34), ("K15", k15, 34), ("G15", g15, 33), ("K31", k31, 33)]
    println(">> $name <<")

    local Φ = hai.map_to_reference(segment)
    local j =
        hai.abs_det_jac(hai.reference_orthotope(BigFloat, 1)) / hai.abs_det_jac(segment)

    fmt_node = Format("[\"%.$(n)e\"],")
    for x in cbt[:nodes]
        local v = Φ(SVector{1}([x]))[1]
        println(format(fmt_node, round(v; sigdigits=n + 1, base=10)))
    end
    println()

    fmt_weight = Format("\"%.$(n)e\",")
    for w in cbt[:weights]
        println(format(fmt_weight, round(j * w; sigdigits=n + 1, base=10)))
    end
    println()
end
