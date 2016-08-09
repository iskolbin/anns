#ifndef ANNS_TASK
#define ANNS_TASK

#define PROBLEM_SWITCH  0
#define PROBLEM_SENSORS 160
#define PROBLEM_ERROR   0

#define NP_OMEGA     500
#define NP_INSULATED 20
#define NP_FLUX      20
#define NP_INITIAL   20

#if   (PROBLEM_SWITCH == 0)
    #if   ((PROBLEM_SENSORS == 160) && (PROBLEM_ERROR == 0))
    double Usensors[] = {1.000556,1.001755,1.003448,1.005462,1.007669,1.009986,1.012363,1.014773,1.017199,1.019632,1.022067,1.024503,1.026940,1.029377,1.031816,1.034256,1.036700,1.039149,1.041603,1.044063,1.046531,1.049008,1.051493,1.053988,1.056493,1.059010,1.061538,1.064078,1.066630,1.069196,1.071774,1.074366,1.076972,1.079592,1.082226,1.084875,1.087538,1.090216,1.092909,1.095618,1.098341,1.101079,1.103833,1.106602,1.109387,1.112187,1.115003,1.117835,1.120682,1.123544,1.126422,1.129316,1.132226,1.135151,1.138091,1.141048,1.144020,1.147007,1.150010,1.153029,1.156063,1.159113,1.162178,1.165258,1.168354,1.171465,1.174592,1.177733,1.180890,1.184063,1.187250,1.190452,1.193670,1.196902,1.200149,1.203412,1.206689,1.209981,1.213288,1.216609,1.219914,1.223166,1.226340,1.229419,1.232397,1.235269,1.238036,1.240700,1.243262,1.245726,1.248094,1.250370,1.252556,1.254656,1.256672,1.258607,1.260463,1.262242,1.263947,1.265580,1.267143,1.268637,1.270064,1.271427,1.272726,1.273964,1.275141,1.276259,1.277319,1.278323,1.279272,1.280166,1.281007,1.281796,1.282534,1.283222,1.283861,1.284452,1.284995,1.285491,1.285941,1.286347,1.286708,1.287025,1.287300,1.287532,1.287722,1.287871,1.287980,1.288049,1.288078,1.288069,1.288021,1.287936,1.287813,1.287654,1.287458,1.287226,1.286959,1.286657,1.286321,1.285950,1.285545,1.285107,1.284636,1.284133,1.283597,1.283029,1.282430,1.281799,1.281137,1.280445,1.279723,1.278970,1.278188,1.277376,1.276535,1.275666,1.274768,1.273841};
    #elif ((PROBLEM_SENSORS == 160) && (PROBLEM_ERROR == 1))
    double Usensors[] = {0.99532156251820258, 0.99833356200365952, 1.008842665541364, 1.0147193172688358, 1.00473464374809, 1.0069211816142356, 1.0099402940485549, 1.0199088286575344, 1.011384379046196, 1.0181573059833722, 1.0283811425806024, 1.0204025226608144, 1.035001357548426, 1.0373551662541929, 1.0279974394764952, 1.0426464398698563, 1.0344497125171517, 1.0286070450079214, 1.0388417887012105, 1.0435618391176058, 1.0500789623797235, 1.0568463413985585, 1.0502314970696758, 1.0555718827110936, 1.0706780137282099, 1.0652859626139048, 1.0601393602853726, 1.0617885956134685, 1.07918833000088, 1.0676752681974078, 1.0755092518177742, 1.0836379179077362, 1.0902532433297989, 1.087243384383686, 1.0777538262341091, 1.0913114377531228, 1.0869580415011635, 1.0908335806010541, 1.0982055250885079, 1.0811651717192337, 1.0891224497486909, 1.0951491494784933, 1.1066357424803881, 1.1079732929189499, 1.1143140292975988, 1.1090490336839665, 1.1111901053426989, 1.106858399790082, 1.1223009647008031, 1.1151725712387419, 1.1199384352376522, 1.137786431663528, 1.1220511943894498, 1.1317213281533107, 1.134154156337051, 1.1444185431370515, 1.1531017433908133, 1.1469477824974414, 1.1454705699476069, 1.1625602183622876, 1.1590378784417998, 1.1520626217244518, 1.1633051151106031, 1.1592659036401998, 1.1744231151495943, 1.1821007548254208, 1.1667831137395464, 1.1713348747264416, 1.1759686533888609, 1.1714035255796396, 1.1699147758031001, 1.1858931717289967, 1.1867294898428795, 1.202120000483091, 1.1974096582065019, 1.2043220271608992, 1.2073180144539402, 1.1954486125032857, 1.2156021046342045, 1.2181247268114526, 1.210909361239749, 1.2220835311317209, 1.2330022080655871, 1.2375814692519287, 1.2341653794531005, 1.2423128028365555, 1.2329714970080419, 1.2490189612487663, 1.2482899550080915, 1.2451808956134367, 1.254499528141958, 1.2438547896631145, 1.2517809331203908, 1.2445531405121484, 1.2612842522626908, 1.2636131818845666, 1.2438312977319996, 1.2698845746123746, 1.2630836413958311, 1.2459646142703875, 1.2625141536997575, 1.2762238508918449, 1.2695814475464469, 1.2814411798757803, 1.2767555072165253, 1.2702337737486376, 1.2796038193766961, 1.2826599410557458, 1.2785915582416882, 1.2905852994200453, 1.2680645480034465, 1.2898480392668887, 1.2742961104655521, 1.277814584824065, 1.2778250465150796, 1.2778101735274063, 1.2947338634434111, 1.2831803750940414, 1.2876153411023152, 1.2744340355961605, 1.3019659083669435, 1.2821308167588117, 1.2672242579294892, 1.2835082043748267, 1.3012208212012415, 1.295048194988347, 1.2858190491363639, 1.2883099582797826, 1.2918595885596955, 1.2949248164666129, 1.2869102327106321, 1.2835603767004833, 1.2867278442320418, 1.2898985104463041, 1.2952327952181841, 1.284311806558982, 1.2916834332912617, 1.2867185022964849, 1.2887051459954344, 1.2836383969019092, 1.2729034294553443, 1.2848044529452645, 1.291302529600979, 1.28507142419654, 1.2817313366781122, 1.2782002959870231, 1.2843256826879126, 1.288622145365476, 1.2721419162125753, 1.2879601205009115, 1.2801253792223077, 1.2857747488094842, 1.272529282246273, 1.2654426040510383, 1.2719030125429571, 1.2726649462783721, 1.2770887779469049, 1.276696079800653, 1.2650414662239171, 1.2699726921719565};
    #elif ((PROBLEM_SENSORS == 16)  && (PROBLEM_ERROR == 0))
    double Usensors[] = {1.0196320000000001, 1.044063, 1.069196, 1.095618, 1.1235440000000001, 1.1530290000000001, 1.1840630000000001, 1.2166090000000001, 1.2457260000000001, 1.2655799999999999, 1.2783230000000001, 1.2854909999999999, 1.288049, 1.2866569999999999, 1.2817989999999999, 1.273841};
    #elif ((PROBLEM_SENSORS == 16)  && (PROBLEM_ERROR == 1))
    double Usensors[] = {1.0187751492215169, 1.0456138294939867, 1.0721065372020087, 1.0962151674796248, 1.1231386421727596, 1.1496194280540988, 1.1845587003031999, 1.2116309155742309, 1.240664518389694, 1.2649929094843377, 1.2819340391076903, 1.280425083481155, 1.2874944613486867, 1.2901594596128307, 1.279274725180696, 1.2725630131484753};
    #endif
#elif (PROBLEM_SWITCH == 1)
    #if   ((PROBLEM_SENSORS == 160) && (PROBLEM_ERROR == 0))
    double Usensors[] = {1.001288,1.004053,1.007940,1.012540,1.017552,1.022779,1.028106,1.033465,1.038817,1.044141,1.049427,1.054668,1.059861,1.065006,1.070102,1.075151,1.080152,1.085105,1.090011,1.094869,1.099679,1.104439,1.109149,1.113807,1.118412,1.122961,1.127453,1.131886,1.136257,1.140565,1.144805,1.148977,1.153078,1.157105,1.161055,1.164926,1.168715,1.172421,1.176040,1.179570,1.183009,1.186355,1.189606,1.192759,1.195814,1.198768,1.201619,1.204367,1.207011,1.209549,1.211980,1.214304,1.216521,1.218629,1.220629,1.222522,1.224306,1.225984,1.227555,1.229020,1.230381,1.231638,1.232794,1.233850,1.234807,1.235669,1.236436,1.237112,1.237700,1.238201,1.238619,1.238957,1.239219,1.239407,1.239526,1.239578,1.239569,1.239501,1.239379,1.239206,1.238988,1.238728,1.238431,1.238102,1.237744,1.237362,1.236961,1.236545,1.236120,1.235689,1.235257,1.234828,1.234408,1.234001,1.233610,1.233241,1.232898,1.232585,1.232306,1.232064,1.231865,1.231711,1.231607,1.231555,1.231559,1.231623,1.231749,1.231941,1.232200,1.232530,1.232932,1.233410,1.233964,1.234597,1.235309,1.236104,1.236980,1.237940,1.238985,1.240113,1.241326,1.242625,1.244007,1.245474,1.247024,1.248657,1.250371,1.252166,1.254039,1.255990,1.258016,1.260115,1.262285,1.264523,1.266827,1.269194,1.271621,1.274105,1.276642,1.279230,1.281864,1.284541,1.287258,1.290010,1.292793,1.295604,1.298437,1.301290,1.304157,1.307034,1.309918,1.312803,1.315685,1.318560,1.321423,1.324270,1.327097,1.329900,1.332673,1.335414};
    #elif ((PROBLEM_SENSORS == 160) && (PROBLEM_ERROR == 1))
    double Usensors[] = {1.0046623344453949, 1.0040930243117223, 1.006153064399594, 1.0169001734720495, 1.0201794602100223, 1.0224566013427119, 1.035917466838163, 1.0304003227073102, 1.0408190412167613, 1.0515800164115223, 1.0451743443835133, 1.0553495682011131, 1.0633875811303248, 1.0678851691077176, 1.0746835921427589, 1.077479258815244, 1.0785699503082506, 1.0860063074363802, 1.0947717640914123, 1.0972314840209283, 1.1020794599916695, 1.1050216233619594, 1.1105846565585236, 1.1201265105918397, 1.1208782898796967, 1.1229893203149719, 1.1323989516601276, 1.1346107019686003, 1.1340354936061519, 1.1348144224277317, 1.1428023929216378, 1.1397068522972917, 1.1525317361874394, 1.1586095826442162, 1.1529862368610901, 1.1694123350287862, 1.1726092499582705, 1.1765386716339534, 1.1822346424591577, 1.1811794185312019, 1.181348944384734, 1.1820462635214626, 1.1878143615664276, 1.1932301284637525, 1.1963417622789811, 1.1965433508077516, 1.1949360506763376, 1.2016108100181713, 1.2050177225135812, 1.2087391872784237, 1.2102511200299564, 1.2157645077838495, 1.2174033510409745, 1.2176092176718902, 1.2187931654447122, 1.2245807814810066, 1.2258703244169562, 1.2289263649784743, 1.2282338444746439, 1.2251554441930192, 1.2343310439949184, 1.2308911602195025, 1.2321551598003002, 1.2373927257285455, 1.2358233357272843, 1.2405591132966667, 1.2335252915164836, 1.2404421145359912, 1.2383813945072197, 1.2370798579400355, 1.2422421703435453, 1.2344230247772516, 1.244095383536602, 1.237140315449154, 1.2456924567728742, 1.2427391388955427, 1.2355722395386959, 1.2382849456218854, 1.2351827130442865, 1.2338449730252148, 1.2440538681627811, 1.24073419278467, 1.2353773543938118, 1.2390751290088327, 1.2334264839373366, 1.2363757905285582, 1.2353365877779741, 1.2384207374615757, 1.2326620430452822, 1.2373347610672614, 1.2300156326687948, 1.2408156772490655, 1.2335278111735359, 1.2306980782371399, 1.2304129081001516, 1.2382987142798414, 1.228936293356883, 1.2318124090106941, 1.23664358574744, 1.2348986086441143, 1.236957829033464, 1.2323590724510813, 1.2310151762078057, 1.2286736565131713, 1.2364356837160615, 1.2289582403712302, 1.2323001760610905, 1.2338076460278271, 1.2359435931556384, 1.2307257918696441, 1.229295529555734, 1.2350691974577328, 1.2339992555475714, 1.2384280096760047, 1.2334041778845899, 1.23541547989567, 1.2348257596304149, 1.2406963312950097, 1.2377922038779958, 1.2380676255418481, 1.2402060370824242, 1.2404360896850768, 1.237524725884328, 1.2481373569328822, 1.2426702992121821, 1.2415755951934397, 1.2497583314162575, 1.252809034027178, 1.2545540343573067, 1.2573087909476119, 1.2595374272737783, 1.2622363045460179, 1.2601145575874122, 1.2573740794798007, 1.2684693759582448, 1.2698246498386605, 1.2691346983202032, 1.2698725077189954, 1.2757928590568108, 1.279026277551556, 1.2793626540450607, 1.2832043377384168, 1.2858276489312173, 1.2868888702902783, 1.291497660540377, 1.2949088262843669, 1.2947384882993513, 1.30136652949759, 1.3080668958913935, 1.3077054579480591, 1.312723308754365, 1.3072647345055879, 1.3163188221283579, 1.3268339203124624, 1.3192061100947021, 1.3264339026477288, 1.3251345450087022, 1.3276263969622835, 1.3358928298389756, 1.3317265219855066};
    #elif ((PROBLEM_SENSORS == 16)  && (PROBLEM_ERROR == 0))
    double Usensors[] = {1.044141, 1.0948690000000001, 1.1405650000000001, 1.17957, 1.209549, 1.22902, 1.2382010000000001, 1.239206, 1.235689, 1.232064, 1.2325299999999999, 1.240113, 1.2559899999999999, 1.2792300000000001, 1.307034, 1.3354140000000001};
    #elif ((PROBLEM_SENSORS == 16)  && (PROBLEM_ERROR == 1))
    double Usensors[] = {1.0457970109438164, 1.1054477226072723, 1.1417222753025, 1.1595755460684396, 1.2067907776339739, 1.2296617988162584, 1.2399588348978321, 1.2377570610569824, 1.2355633884628898, 1.238119677240213, 1.2419901690662787, 1.2445799313862442, 1.2513822915723283, 1.2806583917936045, 1.2898094361222339, 1.3395126248012141};
    #endif
#elif (PROBLEM_SWITCH == 2)
    #if   ((PROBLEM_SENSORS == 160) && (PROBLEM_ERROR == 0))
    #elif ((PROBLEM_SENSORS == 160) && (PROBLEM_ERROR == 1))
    #elif ((PROBLEM_SENSORS == 16)  && (PROBLEM_ERROR == 0))
    #elif ((PROBLEM_SENSORS == 16)  && (PROBLEM_ERROR == 1))
    #endif
#endif

#define X_m_ 10.0
#define T_m  1.0
#define U_m  1.0

#define X_o (X_m_*0.01)
#define T_o (T_m/PROBLEM_SENSORS)

#define X_m (X_m_-X_o)

#define TAU (T_m/(FD_STEPS-1))

#define X_1 (X_m*0.9)

#define C_alpha ( (T_m/(X_m*X_m)) * 202.33576642336 )
#define C_beta  ( (U_m/X_m) * 0.24912184549463 )

int data[] = {
    PJ_INVERSE_UNKNOWN_BOUNDARY, PJ_PARABOLIC_PROBLEM, 3, 0,

    FUNCTIONAL_SQUARE, 2, 5, 2, 0,

    NP_OMEGA, 1,   // ��������� �������� (v[0] - ��������������� �������� ��������)
    NP_INSULATED, 0,   // ������������� ������� (x = 0)
    NP_FLUX, 0,   // ����������� �������
    NP_INITIAL, 0,   // ��������� ������� (t = 0)

    PROBLEM_SENSORS, 1,

    NORMALIZED | RADIAL | GAUSSIAN, 16, 2, 0,      2,  1, 1,  0,  0,
    NORMALIZED | RADIAL | GAUSSIAN, 8, 1, 1, 1,   -1, -1, 0, -1, -1,
};

double bounds[] = {
    GP_UNIFORM_SOLID_CUBOID, X_o, X_m,  T_o, T_m,  // ��������� ��������
    GP_UNIFORM_SOLID_CUBOID, 0.0, 0.0,  T_o, T_m,  // ������������� �������
    GP_UNIFORM_SOLID_CUBOID, X_m_, X_m_,  T_o, T_m,  // ����������� �������
    GP_UNIFORM_SOLID_CUBOID, X_o, X_m,  0.0, 0.0,        // ��������� �������
    GP_GRID_2D,              1, X_1, X_1,   PROBLEM_SENSORS, T_o, T_m,  // ����������� �������
};

#if   (PROBLEM_SWITCH == 0)
    ANNS_FUN(F_q) {double t = x[1]/T_m; return C_beta * ((t <= 0.5) ? (3.0 + 14.0*t) : (10 - 14.0*(t-0.5)));}
#elif (PROBLEM_SWITCH == 1)
    ANNS_FUN(F_q) {double t = x[1]/T_m; return C_beta * (7.0 + 3.0*sin(2.5*M_PI*t));}
#elif (PROBLEM_SWITCH == 2)
    ANNS_FUN(F_q) {double t = x[1]/T_m; return C_beta * ((t <= 0.5) ? 10.0 : 5.0);}
#endif

ANNS_FUN(F_u) {double x_ = x[0]/X_m; return C_alpha * x_ * (1-x_);}

ANNS_CND(Ctransfer)  {return u_xx[0][0] - v[0]*u_x[0][1];}
ANNS_CND(Cinsulated) {return u_x[0][0];}
ANNS_CND(Cheatflux)  {return u_x[0][0] - u[1];}
ANNS_CND(Cinitial)   {return u[0] - U_m;}
ANNS_CND(Cs) {return u[0] - v[0];}

ANNS_CND_G(Ctransfer_a)  {return u_axx[0] - v[0]*u_ax[1];}
ANNS_CND_G(Cinsulated_a) {return u_ax[0];}
ANNS_CND_G(Cheatflux_a_1)  {return u_ax[0];}
ANNS_CND_G(Cheatflux_a_2)  {return -u_a;}
ANNS_CND_G(Cinitial_a)   {return u_a;}
ANNS_CND_G(Cs_a) {return u_a;}

fbvp2_t cond[] = {
    Ctransfer,
    Cinsulated,
    Cheatflux,
    Cinitial,
    Cs};

gbvp2_t cond_a[] = {
    Ctransfer_a,
    Cinsulated_a,
    Cheatflux_a_1,
    Cinitial_a,
    Cs_a,

    NULL, NULL, Cheatflux_a_2, NULL, NULL,};

double delta[] = {
    30./NP_OMEGA,
    1./NP_INSULATED,
    10./NP_FLUX,
    1./NP_INITIAL,
    20./PROBLEM_SENSORS};

void ai_init(anns_instance_t *ai) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

void ai_reload(anns_instance_t *ai) {
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 0, 0, F_u, NULL);
    ai_load_v(ai, 4, 0, Usensors);
}

void ai_initann(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;

    for (i = ai->I[0]; i < ai->I[1]; ) {
        x[i++] = uniform(-4, 4);               // ����
        x[i++] = uniform(0.05 * T_m, 1.05 * T_m);  // ������
        x[i++] = uniform(-0.05 * X_m, 1.05 * X_m); // ������
        x[i++] = uniform(-0.05 * T_m, 1.05 * T_m); // ������
    }

    for (i = ai->I[1]; i < ai->I[2]; ) {
        x[i++] = 0.25;               // ����
        x[i++] = uniform(0.05 * T_m, 0.25 * T_m);  // ������
        x[i++] = uniform(-0.05 * T_m, 1.05 * T_m); // ������
    }
}

double ai_mse(anns_instance_t *ai, double *x) {
    int npoints = 101, i;
    double mse = 0, tmp, hx = T_m/(npoints-1), point[] = {X_m, 0.0};

    for (i = 0; i < npoints; i++) {
        tmp = anns_evalsingle(1, point, x, ai);
        tmp -= F_q(point, 2, NULL);
        point[1] += hx;
        mse += tmp*tmp;
    }
    return sqrt(mse / (npoints-1));
}

#endif // ANNS_TASK