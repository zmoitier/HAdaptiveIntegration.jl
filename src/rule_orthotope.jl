"""
Gauss-Kronrod with 15 nodes.
"""
const SEGMENT_GK15 = TabulatedEmbeddedCubature{Segment}(;
    description="Gauss-Kronrod with 15 nodes (SEGMENT_GK15)",
    reference="https://www.advanpix.com/2011/11/07/gauss-kronrod-quadrature-nodes-weights/",
    nb_significant_digits=34,
    nodes=[
        ["5.0000000000000000000000000000000000e-01"],
        ["2.9707742431130141654669679396151925e-01"],
        ["7.0292257568869858345330320603848075e-01"],
        ["1.2923440720030278006806761335960580e-01"],
        ["8.7076559279969721993193238664039420e-01"],
        ["2.5446043828620737736905157976074350e-02"],
        ["9.7455395617137926226309484202392565e-01"],
        ["3.9610752249605076619965529811337755e-01"],
        ["6.0389247750394923380034470188662245e-01"],
        ["2.0695638226615443485292758087063520e-01"],
        ["7.9304361773384556514707241912936480e-01"],
        ["6.7567788320115463605143605679536900e-02"],
        ["9.3243221167988453639485639432046310e-01"],
        ["4.2723144395936803965726512368357500e-03"],
        ["9.9572768556040631960342734876316425e-01"],
    ],
    weights_high=[
        "1.0474107054236391400649958744585715e-01",
        "9.5175289032392704956628201210506850e-02",
        "9.5175289032392704956628201210506850e-02",
        "7.0326629857762959372594795255118950e-02",
        "7.0326629857762959372594795255118950e-02",
        "3.1546046314989276645350331594602145e-02",
        "3.1546046314989276645350331594602145e-02",
        "1.0221647003764944620708099961732455e-01",
        "1.0221647003764944620708099961732455e-01",
        "8.4502363319633951413291713299275150e-02",
        "8.4502363319633951413291713299275150e-02",
        "5.2395005161125091919938161270759000e-02",
        "5.2395005161125091919938161270759000e-02",
        "1.1467661005264612481866004029484795e-02",
        "1.1467661005264612481866004029484795e-02",
    ],
    order_high=23,
    weights_low=[
        "2.0897959183673469387755102040816325e-01",
        "1.9091502525255947247518488774448755e-01",
        "1.9091502525255947247518488774448755e-01",
        "1.3985269574463833395073388571188980e-01",
        "1.3985269574463833395073388571188980e-01",
        "6.4742483084434846635305716339541000e-02",
        "6.4742483084434846635305716339541000e-02",
    ],
    order_low=13,
)

"""
Gauss-Kronrod with 31 nodes.
"""
const SEGMENT_GK31 = TabulatedEmbeddedCubature{Segment}(;
    description="Gauss-Kronrod with 31 nodes",
    reference="https://www.advanpix.com/2011/11/07/gauss-kronrod-quadrature-nodes-weights/",
    nb_significant_digits=33,
    nodes=[
        ["5.000000000000000000000000000000000e-01"],
        ["3.994029530012827388496858483027019e-01"],
        ["6.005970469987172611503141516972981e-01"],
        ["3.029243264612183150513963145094772e-01"],
        ["6.970756735387816849486036854905228e-01"],
        ["2.145139136957305762313866313730447e-01"],
        ["7.854860863042694237686133686269553e-01"],
        ["1.377911343199149762919069726930310e-01"],
        ["8.622088656800850237080930273069690e-01"],
        ["7.589670829478639189967583961289155e-02"],
        ["9.241032917052136081003241603871084e-01"],
        ["3.136330379964704784612052614489525e-02"],
        ["9.686366962003529521538794738551048e-01"],
        ["6.003740989757285755217140706693700e-03"],
        ["9.939962590102427142447828592933063e-01"],
        ["4.494289665406412504864628842763038e-01"],
        ["5.505710334593587495135371157236962e-01"],
        ["3.504099964234155939166099878668055e-01"],
        ["6.495900035765844060833900121331945e-01"],
        ["2.574590681798801596531721298838247e-01"],
        ["7.425409318201198403468278701161753e-01"],
        ["1.745016293512915147331320523433627e-01"],
        ["8.254983706487084852668679476566373e-01"],
        ["1.047907492787670335161753525910263e-01"],
        ["8.952092507212329664838246474089736e-01"],
        ["5.136773382795904955874517177275205e-02"],
        ["9.486322661720409504412548282272480e-01"],
        ["1.613046216043043287132601060783140e-02"],
        ["9.838695378395695671286739893921686e-01"],
        ["9.988506533014698574135799238644000e-04"],
        ["9.990011493466985301425864200761356e-01"],
    ],
    weights_high=[
        "5.066500350739577450868739638374625e-02",
        "4.958679936089597966619658674230156e-02",
        "4.958679936089597966619658674230156e-02",
        "4.656329908541266061274343637367286e-02",
        "4.656329908541266061274343637367286e-02",
        "4.154025141156651051914462364305190e-02",
        "4.154025141156651051914462364305190e-02",
        "3.492706065936412935476003854957374e-02",
        "3.492706065936412935476003854957374e-02",
        "2.674076234546404363267157361971515e-02",
        "2.674076234546404363267157361971515e-02",
        "1.767318039568792311101897423918002e-02",
        "1.767318039568792311101897423918002e-02",
        "7.503973664658061269187381537903635e-03",
        "7.503973664658061269187381537903635e-03",
        "5.038492276193779752247333130878485e-02",
        "5.038492276193779752247333130878485e-02",
        "4.832136349181183925258995381379467e-02",
        "4.832136349181183925258995381379467e-02",
        "4.428222152810588532363772184688715e-02",
        "4.428222152810588532363772184688715e-02",
        "3.842484037886018944721638874132950e-02",
        "3.842484037886018944721638874132950e-02",
        "3.100478390033532014256961548040146e-02",
        "3.100478390033532014256961548040146e-02",
        "2.229487566238243830411364968663984e-02",
        "2.229487566238243830411364968663984e-02",
        "1.273042366335766009343700050982668e-02",
        "1.273042366335766009343700050982668e-02",
        "2.688739936461674493896025715063825e-03",
        "2.688739936461674493896025715063825e-03",
    ],
    order_high=47,
    weights_low=[
        "1.012891209627806364403100999837596e-01",
        "9.921574266355578822805916322191965e-02",
        "9.921574266355578822805916322191965e-02",
        "9.308050000778110551340028093321140e-02",
        "9.308050000778110551340028093321140e-02",
        "8.313460290849696677660043024060440e-02",
        "8.313460290849696677660043024060440e-02",
        "6.978533896307715722390239725551415e-02",
        "6.978533896307715722390239725551415e-02",
        "5.357961023358596750593477334293465e-02",
        "5.357961023358596750593477334293465e-02",
        "3.518302374405406235463370822533367e-02",
        "3.518302374405406235463370822533367e-02",
        "1.537662099805863417731419678860221e-02",
        "1.537662099805863417731419678860221e-02",
    ],
    order_low=29,
)

"""
Cools and Haegemans with 25 nodes.
"""
const SQUARE_CH25 = TabulatedEmbeddedCubature{Rectangle}(;
    description="Cools and Haegemans with 25 nodes (SQUARE_CH25)",
    reference="https://link.springer.com/article/10.1007/BF01389339",
    nb_significant_digits=19,
    nodes=[
        ["9.5308992296933198185e-01", "9.5308992296933198185e-01"],
        ["4.6910077030668018150e-02", "9.5308992296933198185e-01"],
        ["4.6910077030668018150e-02", "4.6910077030668018150e-02"],
        ["9.5308992296933198185e-01", "4.6910077030668018150e-02"],
        ["9.5308992296933198185e-01", "7.6923465505284155386e-01"],
        ["4.6910077030668018150e-02", "7.6923465505284155386e-01"],
        ["4.6910077030668018150e-02", "2.3076534494715844614e-01"],
        ["9.5308992296933198185e-01", "2.3076534494715844614e-01"],
        ["7.6923465505284155386e-01", "9.5308992296933198185e-01"],
        ["7.6923465505284155386e-01", "4.6910077030668018150e-02"],
        ["2.3076534494715844614e-01", "4.6910077030668018150e-02"],
        ["2.3076534494715844614e-01", "9.5308992296933198185e-01"],
        ["7.6923465505284155386e-01", "7.6923465505284155386e-01"],
        ["2.3076534494715844614e-01", "7.6923465505284155386e-01"],
        ["2.3076534494715844614e-01", "2.3076534494715844614e-01"],
        ["7.6923465505284155386e-01", "2.3076534494715844614e-01"],
        ["9.5308992296933198185e-01", "5.0000000000000000000e-01"],
        ["4.6910077030668018150e-02", "5.0000000000000000000e-01"],
        ["5.0000000000000000000e-01", "9.5308992296933198185e-01"],
        ["5.0000000000000000000e-01", "4.6910077030668018150e-02"],
        ["5.0000000000000000000e-01", "5.0000000000000000000e-01"],
        ["7.6923465505284155386e-01", "5.0000000000000000000e-01"],
        ["2.3076534494715844614e-01", "5.0000000000000000000e-01"],
        ["5.0000000000000000000e-01", "7.6923465505284155386e-01"],
        ["5.0000000000000000000e-01", "2.3076534494715844614e-01"],
    ],
    weights_high=[
        "1.4033587215607158679e-02",
        "1.4033587215607158679e-02",
        "1.4033587215607158679e-02",
        "1.4033587215607158679e-02",
        "2.8349999999999999855e-02",
        "2.8349999999999999855e-02",
        "2.8349999999999999855e-02",
        "2.8349999999999999855e-02",
        "2.8349999999999999855e-02",
        "2.8349999999999999855e-02",
        "2.8349999999999999855e-02",
        "2.8349999999999999855e-02",
        "5.7271351055997779968e-02",
        "5.7271351055997779968e-02",
        "5.7271351055997779968e-02",
        "5.7271351055997779968e-02",
        "3.3696268096880225377e-02",
        "3.3696268096880225377e-02",
        "3.3696268096880225377e-02",
        "3.3696268096880225377e-02",
        "8.0908641975308641832e-02",
        "6.8071633137687675802e-02",
        "6.8071633137687675802e-02",
        "6.8071633137687675802e-02",
        "6.8071633137687675802e-02",
    ],
    order_high=9,
    weights_low=[
        "2.0593268489242785292e-02",
        "2.0593268489242785292e-02",
        "2.0593268489242785292e-02",
        "2.0593268489242785292e-02",
        "9.7723992924246520540e-03",
        "9.7723992924246520540e-03",
        "9.7723992924246520540e-03",
        "9.7723992924246520540e-03",
        "9.7723992924246520540e-03",
        "9.7723992924246520540e-03",
        "9.7723992924246520540e-03",
        "9.7723992924246520540e-03",
        "1.0988476833241696496e-01",
        "1.0988476833241696496e-01",
        "1.0988476833241696496e-01",
        "1.0988476833241696496e-01",
        "5.7732106964759669065e-02",
        "5.7732106964759669065e-02",
        "5.7732106964759669065e-02",
        "5.7732106964759669065e-02",
        "1.6898023051492500631e-01",
    ],
    order_low=7,
)

"""
Cools and Haegemans with 21 nodes.
"""
const SQUARE_CH21 = TabulatedEmbeddedCubature{Rectangle}(;
    description="Cools and Haegemans with 21 nodes (SQUARE_CH21)",
    reference="https://link.springer.com/article/10.1007/BF01389339",
    nb_significant_digits=19,
    nodes=[
        ["5.0000000000000000000e-01", "5.0000000000000000000e-01"],
        ["9.5308992296933199640e-01", "5.0000000000000000000e-01"],
        ["5.0000000000000000000e-01", "9.5308992296933199640e-01"],
        ["4.6910077030668003600e-02", "5.0000000000000000000e-01"],
        ["5.0000000000000000000e-01", "4.6910077030668003600e-02"],
        ["7.6923465505284154552e-01", "7.6923465505284154552e-01"],
        ["2.3076534494715845448e-01", "7.6923465505284154552e-01"],
        ["7.6923465505284154552e-01", "2.3076534494715845448e-01"],
        ["2.3076534494715845448e-01", "2.3076534494715845448e-01"],
        ["9.5308992296933199640e-01", "9.5308992296933199640e-01"],
        ["4.6910077030668003600e-02", "9.5308992296933199640e-01"],
        ["9.5308992296933199640e-01", "4.6910077030668003600e-02"],
        ["4.6910077030668003600e-02", "4.6910077030668003600e-02"],
        ["9.5308992296933199640e-01", "7.6923465505284154552e-01"],
        ["7.6923465505284154552e-01", "9.5308992296933199640e-01"],
        ["4.6910077030668003600e-02", "7.6923465505284154552e-01"],
        ["7.6923465505284154552e-01", "4.6910077030668003600e-02"],
        ["9.5308992296933199640e-01", "2.3076534494715845448e-01"],
        ["2.3076534494715845448e-01", "9.5308992296933199640e-01"],
        ["4.6910077030668003600e-02", "2.3076534494715845448e-01"],
        ["2.3076534494715845448e-01", "4.6910077030668003600e-02"],
    ],
    weights_high=[
        "1.6898023051492500631e-01",
        "5.7732106964759669065e-02",
        "5.7732106964759669065e-02",
        "5.7732106964759669065e-02",
        "5.7732106964759669065e-02",
        "1.0988476833241696496e-01",
        "1.0988476833241696496e-01",
        "1.0988476833241696496e-01",
        "1.0988476833241696496e-01",
        "2.0593268489242785292e-02",
        "2.0593268489242785292e-02",
        "2.0593268489242785292e-02",
        "2.0593268489242785292e-02",
        "9.7723992924246520540e-03",
        "9.7723992924246520540e-03",
        "9.7723992924246520540e-03",
        "9.7723992924246520540e-03",
        "9.7723992924246520540e-03",
        "9.7723992924246520540e-03",
        "9.7723992924246520540e-03",
        "9.7723992924246520540e-03",
    ],
    order_high=7,
    weights_low=[
        "1.5262184183613067345e-01",
        "6.5911301304156885498e-02",
        "6.5911301304156885498e-02",
        "6.5911301304156885498e-02",
        "6.5911301304156885498e-02",
        "1.1965716762484161701e-01",
        "1.1965716762484161701e-01",
        "1.1965716762484161701e-01",
        "1.1965716762484161701e-01",
        "2.6276070611968829130e-02",
        "2.6276070611968829130e-02",
        "2.6276070611968829130e-02",
        "2.6276070611968829130e-02",
    ],
    order_low=5,
)

"""
Genz and Malik with 17 nodes.
"""
const SQUARE_GM17 = TabulatedEmbeddedCubature{Rectangle}(;
    description="Genz and Malik with 17 nodes (SQUARE_GM17)",
    reference="https://doi.org/10.1016/0771-050X(80)90039-X",
    nb_significant_digits=35,
    nodes=[
        [
            "5.000000000000000000000000000000000000e-01",
            "5.000000000000000000000000000000000000e-01",
        ],
        [
            "6.792842914001590459953225769539687477e-01",
            "5.000000000000000000000000000000000000e-01",
        ],
        [
            "3.207157085998409540046774230460312523e-01",
            "5.000000000000000000000000000000000000e-01",
        ],
        [
            "9.743416490252568997998340316649077801e-01",
            "5.000000000000000000000000000000000000e-01",
        ],
        [
            "2.565835097474310020016596833509221994e-02",
            "5.000000000000000000000000000000000000e-01",
        ],
        [
            "5.000000000000000000000000000000000000e-01",
            "6.792842914001590459953225769539687477e-01",
        ],
        [
            "5.000000000000000000000000000000000000e-01",
            "3.207157085998409540046774230460312523e-01",
        ],
        [
            "5.000000000000000000000000000000000000e-01",
            "9.743416490252568997998340316649077801e-01",
        ],
        [
            "5.000000000000000000000000000000000000e-01",
            "2.565835097474310020016596833509221994e-02",
        ],
        [
            "9.743416490252568997998340316649077801e-01",
            "9.743416490252568997998340316649077801e-01",
        ],
        [
            "2.565835097474310020016596833509221994e-02",
            "9.743416490252568997998340316649077801e-01",
        ],
        [
            "9.743416490252568997998340316649077801e-01",
            "2.565835097474310020016596833509221994e-02",
        ],
        [
            "2.565835097474310020016596833509221994e-02",
            "2.565835097474310020016596833509221994e-02",
        ],
        [
            "8.441236008058426488608143671468117626e-01",
            "8.441236008058426488608143671468117626e-01",
        ],
        [
            "1.558763991941573511391856328531882374e-01",
            "8.441236008058426488608143671468117626e-01",
        ],
        [
            "8.441236008058426488608143671468117626e-01",
            "1.558763991941573511391856328531882374e-01",
        ],
        [
            "1.558763991941573511391856328531882374e-01",
            "1.558763991941573511391856328531882374e-01",
        ],
    ],
    weights_high=[
        "-1.938728852309099222679469593049839963e-01",
        "1.493674744703551287913427831123304374e-01",
        "1.493674744703551287913427831123304374e-01",
        "5.182136869379667733577198597774729462e-02",
        "5.182136869379667733577198597774729462e-02",
        "1.493674744703551287913427831123304374e-01",
        "1.493674744703551287913427831123304374e-01",
        "5.182136869379667733577198597774729462e-02",
        "5.182136869379667733577198597774729462e-02",
        "1.016105268505817202662195803485241071e-02",
        "1.016105268505817202662195803485241071e-02",
        "1.016105268505817202662195803485241071e-02",
        "1.016105268505817202662195803485241071e-02",
        "8.711832545851750241325001270131585632e-02",
        "8.711832545851750241325001270131585632e-02",
        "8.711832545851750241325001270131585632e-02",
        "8.711832545851750241325001270131585632e-02",
    ],
    order_high=7,
    weights_low=[
        "-1.331961591220850480109739368998628258e+00",
        "5.041152263374485596707818930041152263e-01",
        "5.041152263374485596707818930041152263e-01",
        "4.458161865569272976680384087791495199e-02",
        "4.458161865569272976680384087791495199e-02",
        "5.041152263374485596707818930041152263e-01",
        "5.041152263374485596707818930041152263e-01",
        "4.458161865569272976680384087791495199e-02",
        "4.458161865569272976680384087791495199e-02",
        "3.429355281207133058984910836762688615e-02",
        "3.429355281207133058984910836762688615e-02",
        "3.429355281207133058984910836762688615e-02",
        "3.429355281207133058984910836762688615e-02",
    ],
    order_low=5,
)

"""
Berntsen and Espelid with 65 nodes.
"""
const CUBE_BE65 = TabulatedEmbeddedCubature{Cuboid}(;
    description="Berntsen and Espelid with 65 nodes (CUBE_BE65)",
    reference="https://epubs.siam.org/doi/10.1137/0725016",
    nb_significant_digits=15,
    nodes=[
        ["5.0000000000000000e-01", "5.0000000000000000e-01", "5.0000000000000000e-01"],
        ["7.9824398257170165e-01", "5.0000000000000000e-01", "5.0000000000000000e-01"],
        ["5.0000000000000000e-01", "7.9824398257170165e-01", "5.0000000000000000e-01"],
        ["5.0000000000000000e-01", "5.0000000000000000e-01", "7.9824398257170165e-01"],
        ["2.0175601742829835e-01", "5.0000000000000000e-01", "5.0000000000000000e-01"],
        ["5.0000000000000000e-01", "2.0175601742829835e-01", "5.0000000000000000e-01"],
        ["5.0000000000000000e-01", "5.0000000000000000e-01", "2.0175601742829835e-01"],
        ["9.5575373953655815e-01", "5.0000000000000000e-01", "5.0000000000000000e-01"],
        ["5.0000000000000000e-01", "9.5575373953655815e-01", "5.0000000000000000e-01"],
        ["5.0000000000000000e-01", "5.0000000000000000e-01", "9.5575373953655815e-01"],
        ["4.4246260463441850e-02", "5.0000000000000000e-01", "5.0000000000000000e-01"],
        ["5.0000000000000000e-01", "4.4246260463441850e-02", "5.0000000000000000e-01"],
        ["5.0000000000000000e-01", "5.0000000000000000e-01", "4.4246260463441850e-02"],
        ["9.2871014331657190e-01", "9.2871014331657190e-01", "5.0000000000000000e-01"],
        ["9.2871014331657190e-01", "5.0000000000000000e-01", "9.2871014331657190e-01"],
        ["5.0000000000000000e-01", "9.2871014331657190e-01", "9.2871014331657190e-01"],
        ["7.1289856683428100e-02", "9.2871014331657190e-01", "5.0000000000000000e-01"],
        ["7.1289856683428100e-02", "5.0000000000000000e-01", "9.2871014331657190e-01"],
        ["5.0000000000000000e-01", "7.1289856683428100e-02", "9.2871014331657190e-01"],
        ["9.2871014331657190e-01", "7.1289856683428100e-02", "5.0000000000000000e-01"],
        ["9.2871014331657190e-01", "5.0000000000000000e-01", "7.1289856683428100e-02"],
        ["5.0000000000000000e-01", "9.2871014331657190e-01", "7.1289856683428100e-02"],
        ["7.1289856683428100e-02", "7.1289856683428100e-02", "5.0000000000000000e-01"],
        ["7.1289856683428100e-02", "5.0000000000000000e-01", "7.1289856683428100e-02"],
        ["5.0000000000000000e-01", "7.1289856683428100e-02", "7.1289856683428100e-02"],
        ["7.5276599277131730e-01", "7.5276599277131730e-01", "7.5276599277131730e-01"],
        ["2.4723400722868270e-01", "7.5276599277131730e-01", "7.5276599277131730e-01"],
        ["7.5276599277131730e-01", "2.4723400722868270e-01", "7.5276599277131730e-01"],
        ["2.4723400722868270e-01", "2.4723400722868270e-01", "7.5276599277131730e-01"],
        ["7.5276599277131730e-01", "7.5276599277131730e-01", "2.4723400722868270e-01"],
        ["2.4723400722868270e-01", "7.5276599277131730e-01", "2.4723400722868270e-01"],
        ["7.5276599277131730e-01", "2.4723400722868270e-01", "2.4723400722868270e-01"],
        ["2.4723400722868270e-01", "2.4723400722868270e-01", "2.4723400722868270e-01"],
        ["9.5147762226420635e-01", "9.5147762226420635e-01", "9.5147762226420635e-01"],
        ["4.8522377735793650e-02", "9.5147762226420635e-01", "9.5147762226420635e-01"],
        ["9.5147762226420635e-01", "4.8522377735793650e-02", "9.5147762226420635e-01"],
        ["4.8522377735793650e-02", "4.8522377735793650e-02", "9.5147762226420635e-01"],
        ["9.5147762226420635e-01", "9.5147762226420635e-01", "4.8522377735793650e-02"],
        ["4.8522377735793650e-02", "9.5147762226420635e-01", "4.8522377735793650e-02"],
        ["9.5147762226420635e-01", "4.8522377735793650e-02", "4.8522377735793650e-02"],
        ["4.8522377735793650e-02", "4.8522377735793650e-02", "4.8522377735793650e-02"],
        ["7.6250000000000000e-01", "7.6250000000000000e-01", "9.6750000000000000e-01"],
        ["7.6250000000000000e-01", "9.6750000000000000e-01", "7.6250000000000000e-01"],
        ["9.6750000000000000e-01", "7.6250000000000000e-01", "7.6250000000000000e-01"],
        ["2.3750000000000000e-01", "7.6250000000000000e-01", "9.6750000000000000e-01"],
        ["2.3750000000000000e-01", "9.6750000000000000e-01", "7.6250000000000000e-01"],
        ["9.6750000000000000e-01", "2.3750000000000000e-01", "7.6250000000000000e-01"],
        ["7.6250000000000000e-01", "2.3750000000000000e-01", "9.6750000000000000e-01"],
        ["7.6250000000000000e-01", "9.6750000000000000e-01", "2.3750000000000000e-01"],
        ["9.6750000000000000e-01", "7.6250000000000000e-01", "2.3750000000000000e-01"],
        ["2.3750000000000000e-01", "2.3750000000000000e-01", "9.6750000000000000e-01"],
        ["2.3750000000000000e-01", "9.6750000000000000e-01", "2.3750000000000000e-01"],
        ["9.6750000000000000e-01", "2.3750000000000000e-01", "2.3750000000000000e-01"],
        ["7.6250000000000000e-01", "7.6250000000000000e-01", "3.2500000000000000e-02"],
        ["7.6250000000000000e-01", "3.2500000000000000e-02", "7.6250000000000000e-01"],
        ["3.2500000000000000e-02", "7.6250000000000000e-01", "7.6250000000000000e-01"],
        ["2.3750000000000000e-01", "7.6250000000000000e-01", "3.2500000000000000e-02"],
        ["2.3750000000000000e-01", "3.2500000000000000e-02", "7.6250000000000000e-01"],
        ["3.2500000000000000e-02", "2.3750000000000000e-01", "7.6250000000000000e-01"],
        ["7.6250000000000000e-01", "2.3750000000000000e-01", "3.2500000000000000e-02"],
        ["7.6250000000000000e-01", "3.2500000000000000e-02", "2.3750000000000000e-01"],
        ["3.2500000000000000e-02", "7.6250000000000000e-01", "2.3750000000000000e-01"],
        ["2.3750000000000000e-01", "2.3750000000000000e-01", "3.2500000000000000e-02"],
        ["2.3750000000000000e-01", "3.2500000000000000e-02", "2.3750000000000000e-01"],
        ["3.2500000000000000e-02", "2.3750000000000000e-01", "2.3750000000000000e-01"],
    ],
    weights_high=[
        "4.5340290436037275e-03",
        "4.1800060049505412e-02",
        "4.1800060049505412e-02",
        "4.1800060049505412e-02",
        "4.1800060049505412e-02",
        "4.1800060049505412e-02",
        "4.1800060049505412e-02",
        "1.3209778122026900e-02",
        "1.3209778122026900e-02",
        "1.3209778122026900e-02",
        "1.3209778122026900e-02",
        "1.3209778122026900e-02",
        "1.3209778122026900e-02",
        "1.3159017373052862e-02",
        "1.3159017373052862e-02",
        "1.3159017373052862e-02",
        "1.3159017373052862e-02",
        "1.3159017373052862e-02",
        "1.3159017373052862e-02",
        "1.3159017373052862e-02",
        "1.3159017373052862e-02",
        "1.3159017373052862e-02",
        "1.3159017373052862e-02",
        "1.3159017373052862e-02",
        "1.3159017373052862e-02",
        "2.6680584820591875e-02",
        "2.6680584820591875e-02",
        "2.6680584820591875e-02",
        "2.6680584820591875e-02",
        "2.6680584820591875e-02",
        "2.6680584820591875e-02",
        "2.6680584820591875e-02",
        "2.6680584820591875e-02",
        "3.6652379333158925e-03",
        "3.6652379333158925e-03",
        "3.6652379333158925e-03",
        "3.6652379333158925e-03",
        "3.6652379333158925e-03",
        "3.6652379333158925e-03",
        "3.6652379333158925e-03",
        "3.6652379333158925e-03",
        "1.1030506309137748e-02",
        "1.1030506309137748e-02",
        "1.1030506309137748e-02",
        "1.1030506309137748e-02",
        "1.1030506309137748e-02",
        "1.1030506309137748e-02",
        "1.1030506309137748e-02",
        "1.1030506309137748e-02",
        "1.1030506309137748e-02",
        "1.1030506309137748e-02",
        "1.1030506309137748e-02",
        "1.1030506309137748e-02",
        "1.1030506309137748e-02",
        "1.1030506309137748e-02",
        "1.1030506309137748e-02",
        "1.1030506309137748e-02",
        "1.1030506309137748e-02",
        "1.1030506309137748e-02",
        "1.1030506309137748e-02",
        "1.1030506309137748e-02",
        "1.1030506309137748e-02",
        "1.1030506309137748e-02",
        "1.1030506309137748e-02",
        "1.1030506309137748e-02",
    ],
    order_high=9,
    weights_low=[
        "2.0049410275506235e-01",
        "-5.1495158847434000e-02",
        "-5.1495158847434000e-02",
        "-5.1495158847434000e-02",
        "-5.1495158847434000e-02",
        "-5.1495158847434000e-02",
        "-5.1495158847434000e-02",
        "4.4635009008031388e-02",
        "4.4635009008031388e-02",
        "4.4635009008031388e-02",
        "4.4635009008031388e-02",
        "4.4635009008031388e-02",
        "4.4635009008031388e-02",
        "1.8642480739586644e-02",
        "1.8642480739586644e-02",
        "1.8642480739586644e-02",
        "1.8642480739586644e-02",
        "1.8642480739586644e-02",
        "1.8642480739586644e-02",
        "1.8642480739586644e-02",
        "1.8642480739586644e-02",
        "1.8642480739586644e-02",
        "1.8642480739586644e-02",
        "1.8642480739586644e-02",
        "1.8642480739586644e-02",
        "7.0756867601023362e-02",
        "7.0756867601023362e-02",
        "7.0756867601023362e-02",
        "7.0756867601023362e-02",
        "7.0756867601023362e-02",
        "7.0756867601023362e-02",
        "7.0756867601023362e-02",
        "7.0756867601023362e-02",
        "6.3627608247658450e-03",
        "6.3627608247658450e-03",
        "6.3627608247658450e-03",
        "6.3627608247658450e-03",
        "6.3627608247658450e-03",
        "6.3627608247658450e-03",
        "6.3627608247658450e-03",
        "6.3627608247658450e-03",
    ],
    order_low=7,
)

"""
Genz and Malik with 33 nodes.
"""
const CUBE_GM33 = TabulatedEmbeddedCubature{Cuboid}(;
    description="Genz and Malik with 33 nodes (CUBE_GM33)",
    reference="https://doi.org/10.1016/0771-050X(80)90039-X",
    nb_significant_digits=35,
    nodes=[
        [
            "5.000000000000000000000000000000000000e-01",
            "5.000000000000000000000000000000000000e-01",
            "5.000000000000000000000000000000000000e-01",
        ],
        [
            "6.792842914001590459953225769539687477e-01",
            "5.000000000000000000000000000000000000e-01",
            "5.000000000000000000000000000000000000e-01",
        ],
        [
            "3.207157085998409540046774230460312523e-01",
            "5.000000000000000000000000000000000000e-01",
            "5.000000000000000000000000000000000000e-01",
        ],
        [
            "9.743416490252568997998340316649077801e-01",
            "5.000000000000000000000000000000000000e-01",
            "5.000000000000000000000000000000000000e-01",
        ],
        [
            "2.565835097474310020016596833509221994e-02",
            "5.000000000000000000000000000000000000e-01",
            "5.000000000000000000000000000000000000e-01",
        ],
        [
            "5.000000000000000000000000000000000000e-01",
            "6.792842914001590459953225769539687477e-01",
            "5.000000000000000000000000000000000000e-01",
        ],
        [
            "5.000000000000000000000000000000000000e-01",
            "3.207157085998409540046774230460312523e-01",
            "5.000000000000000000000000000000000000e-01",
        ],
        [
            "5.000000000000000000000000000000000000e-01",
            "9.743416490252568997998340316649077801e-01",
            "5.000000000000000000000000000000000000e-01",
        ],
        [
            "5.000000000000000000000000000000000000e-01",
            "2.565835097474310020016596833509221994e-02",
            "5.000000000000000000000000000000000000e-01",
        ],
        [
            "5.000000000000000000000000000000000000e-01",
            "5.000000000000000000000000000000000000e-01",
            "6.792842914001590459953225769539687477e-01",
        ],
        [
            "5.000000000000000000000000000000000000e-01",
            "5.000000000000000000000000000000000000e-01",
            "3.207157085998409540046774230460312523e-01",
        ],
        [
            "5.000000000000000000000000000000000000e-01",
            "5.000000000000000000000000000000000000e-01",
            "9.743416490252568997998340316649077801e-01",
        ],
        [
            "5.000000000000000000000000000000000000e-01",
            "5.000000000000000000000000000000000000e-01",
            "2.565835097474310020016596833509221994e-02",
        ],
        [
            "9.743416490252568997998340316649077801e-01",
            "9.743416490252568997998340316649077801e-01",
            "5.000000000000000000000000000000000000e-01",
        ],
        [
            "2.565835097474310020016596833509221994e-02",
            "9.743416490252568997998340316649077801e-01",
            "5.000000000000000000000000000000000000e-01",
        ],
        [
            "9.743416490252568997998340316649077801e-01",
            "2.565835097474310020016596833509221994e-02",
            "5.000000000000000000000000000000000000e-01",
        ],
        [
            "2.565835097474310020016596833509221994e-02",
            "2.565835097474310020016596833509221994e-02",
            "5.000000000000000000000000000000000000e-01",
        ],
        [
            "9.743416490252568997998340316649077801e-01",
            "5.000000000000000000000000000000000000e-01",
            "9.743416490252568997998340316649077801e-01",
        ],
        [
            "2.565835097474310020016596833509221994e-02",
            "5.000000000000000000000000000000000000e-01",
            "9.743416490252568997998340316649077801e-01",
        ],
        [
            "9.743416490252568997998340316649077801e-01",
            "5.000000000000000000000000000000000000e-01",
            "2.565835097474310020016596833509221994e-02",
        ],
        [
            "2.565835097474310020016596833509221994e-02",
            "5.000000000000000000000000000000000000e-01",
            "2.565835097474310020016596833509221994e-02",
        ],
        [
            "5.000000000000000000000000000000000000e-01",
            "9.743416490252568997998340316649077801e-01",
            "9.743416490252568997998340316649077801e-01",
        ],
        [
            "5.000000000000000000000000000000000000e-01",
            "2.565835097474310020016596833509221994e-02",
            "9.743416490252568997998340316649077801e-01",
        ],
        [
            "5.000000000000000000000000000000000000e-01",
            "9.743416490252568997998340316649077801e-01",
            "2.565835097474310020016596833509221994e-02",
        ],
        [
            "5.000000000000000000000000000000000000e-01",
            "2.565835097474310020016596833509221994e-02",
            "2.565835097474310020016596833509221994e-02",
        ],
        [
            "8.441236008058426488608143671468117626e-01",
            "8.441236008058426488608143671468117626e-01",
            "8.441236008058426488608143671468117626e-01",
        ],
        [
            "1.558763991941573511391856328531882374e-01",
            "8.441236008058426488608143671468117626e-01",
            "8.441236008058426488608143671468117626e-01",
        ],
        [
            "8.441236008058426488608143671468117626e-01",
            "1.558763991941573511391856328531882374e-01",
            "8.441236008058426488608143671468117626e-01",
        ],
        [
            "1.558763991941573511391856328531882374e-01",
            "1.558763991941573511391856328531882374e-01",
            "8.441236008058426488608143671468117626e-01",
        ],
        [
            "8.441236008058426488608143671468117626e-01",
            "8.441236008058426488608143671468117626e-01",
            "1.558763991941573511391856328531882374e-01",
        ],
        [
            "1.558763991941573511391856328531882374e-01",
            "8.441236008058426488608143671468117626e-01",
            "1.558763991941573511391856328531882374e-01",
        ],
        [
            "8.441236008058426488608143671468117626e-01",
            "1.558763991941573511391856328531882374e-01",
            "1.558763991941573511391856328531882374e-01",
        ],
        [
            "1.558763991941573511391856328531882374e-01",
            "1.558763991941573511391856328531882374e-01",
            "1.558763991941573511391856328531882374e-01",
        ],
    ],
    weights_high=[
        "-5.556063608189808464156886653457298176e-01",
        "1.493674744703551287913427831123304374e-01",
        "1.493674744703551287913427831123304374e-01",
        "3.149926332368033328252806990804247320e-02",
        "3.149926332368033328252806990804247320e-02",
        "1.493674744703551287913427831123304374e-01",
        "1.493674744703551287913427831123304374e-01",
        "3.149926332368033328252806990804247320e-02",
        "3.149926332368033328252806990804247320e-02",
        "1.493674744703551287913427831123304374e-01",
        "1.493674744703551287913427831123304374e-01",
        "3.149926332368033328252806990804247320e-02",
        "3.149926332368033328252806990804247320e-02",
        "1.016105268505817202662195803485241071e-02",
        "1.016105268505817202662195803485241071e-02",
        "1.016105268505817202662195803485241071e-02",
        "1.016105268505817202662195803485241071e-02",
        "1.016105268505817202662195803485241071e-02",
        "1.016105268505817202662195803485241071e-02",
        "1.016105268505817202662195803485241071e-02",
        "1.016105268505817202662195803485241071e-02",
        "1.016105268505817202662195803485241071e-02",
        "1.016105268505817202662195803485241071e-02",
        "1.016105268505817202662195803485241071e-02",
        "1.016105268505817202662195803485241071e-02",
        "4.355916272925875120662500635065792816e-02",
        "4.355916272925875120662500635065792816e-02",
        "4.355916272925875120662500635065792816e-02",
        "4.355916272925875120662500635065792816e-02",
        "4.355916272925875120662500635065792816e-02",
        "4.355916272925875120662500635065792816e-02",
        "4.355916272925875120662500635065792816e-02",
        "4.355916272925875120662500635065792816e-02",
    ],
    order_high=7,
    weights_low=[
        "-2.292181069958847736625514403292181070e+00",
        "5.041152263374485596707818930041152263e-01",
        "5.041152263374485596707818930041152263e-01",
        "-2.400548696844993141289437585733882030e-02",
        "-2.400548696844993141289437585733882030e-02",
        "5.041152263374485596707818930041152263e-01",
        "5.041152263374485596707818930041152263e-01",
        "-2.400548696844993141289437585733882030e-02",
        "-2.400548696844993141289437585733882030e-02",
        "5.041152263374485596707818930041152263e-01",
        "5.041152263374485596707818930041152263e-01",
        "-2.400548696844993141289437585733882030e-02",
        "-2.400548696844993141289437585733882030e-02",
        "3.429355281207133058984910836762688615e-02",
        "3.429355281207133058984910836762688615e-02",
        "3.429355281207133058984910836762688615e-02",
        "3.429355281207133058984910836762688615e-02",
        "3.429355281207133058984910836762688615e-02",
        "3.429355281207133058984910836762688615e-02",
        "3.429355281207133058984910836762688615e-02",
        "3.429355281207133058984910836762688615e-02",
        "3.429355281207133058984910836762688615e-02",
        "3.429355281207133058984910836762688615e-02",
        "3.429355281207133058984910836762688615e-02",
        "3.429355281207133058984910836762688615e-02",
    ],
    order_low=5,
)
