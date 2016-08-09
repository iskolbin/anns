#ifndef ANNS_TASK
#define ANNS_TASK

#define N_SENSORS 160

//#define TAU (T_m/N_SENSORS)

double Ub[] = {1.000000,1.081801,1.107909,1.127021,1.142858,1.156801,1.169500,1.181313,1.192457,1.203073,1.213255,1.223070,1.232565,1.241775,1.250725,1.259433,1.267912,1.276171,1.284216,1.292051,1.299678,1.307098,1.314310,1.321312,1.328103,1.334680,1.341041,1.347182,1.353099,1.358791,1.364253,1.369483,1.374477,1.379233,1.383748,1.388020,1.392046,1.395826,1.399357,1.402639,1.405671,1.408453,1.410985,1.413269,1.415304,1.417093,1.418638,1.419941,1.421005,1.421834,1.422431,1.422801,1.422949,1.422879,1.422598,1.422112,1.421426,1.420548,1.419485,1.418245,1.416836,1.415265,1.413543,1.411677,1.409677,1.407553,1.405314,1.402970,1.400531,1.398008,1.395412,1.392752,1.390040,1.387287,1.384503,1.381699,1.378886,1.376076,1.373279,1.370505,1.367766,1.365072,1.362433,1.359860,1.357363,1.354952,1.352636,1.350425,1.348327,1.346352,1.344509,1.342804,1.341246,1.339842,1.338600,1.337525,1.336625,1.335904,1.335367,1.335021,1.334868,1.334913,1.335158,1.335607,1.336263,1.337126,1.338198,1.339479,1.340971,1.342671,1.344581,1.346698,1.349020,1.351545,1.354269,1.357191,1.360305,1.363607,1.367092,1.370756,1.374592,1.378593,1.382754,1.387068,1.391526,1.396122,1.400847,1.405692,1.410649,1.415709,1.420863,1.426100,1.431412,1.436789,1.442219,1.447693,1.453201,1.458732,1.464276,1.469822,1.475359,1.480876,1.486365,1.491812,1.497210,1.502546,1.507811,1.512995,1.518088,1.523080,1.527963,1.532726,1.537361,1.541860,1.546215,1.550416,1.554458,1.558332,1.562032,1.565551,1.568884};
double Us[] = {1.001288,1.004053,1.007940,1.012540,1.017552,1.022779,1.028106,1.033465,1.038817,1.044141,1.049427,1.054668,1.059861,1.065006,1.070102,1.075151,1.080152,1.085105,1.090011,1.094869,1.099679,1.104439,1.109149,1.113807,1.118412,1.122961,1.127453,1.131886,1.136257,1.140565,1.144805,1.148977,1.153078,1.157105,1.161055,1.164926,1.168715,1.172421,1.176040,1.179570,1.183009,1.186355,1.189606,1.192759,1.195814,1.198768,1.201619,1.204367,1.207011,1.209549,1.211980,1.214304,1.216521,1.218629,1.220629,1.222522,1.224306,1.225984,1.227555,1.229020,1.230381,1.231638,1.232794,1.233850,1.234807,1.235669,1.236436,1.237112,1.237700,1.238201,1.238619,1.238957,1.239219,1.239407,1.239526,1.239578,1.239569,1.239501,1.239379,1.239206,1.238988,1.238728,1.238431,1.238102,1.237744,1.237362,1.236961,1.236545,1.236120,1.235689,1.235257,1.234828,1.234408,1.234001,1.233610,1.233241,1.232898,1.232585,1.232306,1.232064,1.231865,1.231711,1.231607,1.231555,1.231559,1.231623,1.231749,1.231941,1.232200,1.232530,1.232932,1.233410,1.233964,1.234597,1.235309,1.236104,1.236980,1.237940,1.238985,1.240113,1.241326,1.242625,1.244007,1.245474,1.247024,1.248657,1.250371,1.252166,1.254039,1.255990,1.258016,1.260115,1.262285,1.264523,1.266827,1.269194,1.271621,1.274105,1.276642,1.279230,1.281864,1.284541,1.287258,1.290010,1.292793,1.295604,1.298437,1.301290,1.304157,1.307034,1.309918,1.312803,1.315685,1.318560,1.321423,1.324270,1.327097,1.329900,1.332673,1.335414};
double Us_err[] = {1.0005093136810457, 1.0046623344453949, 1.0040930243117223, 1.006153064399594, 1.0169001734720495, 1.0201794602100223, 1.0224566013427119, 1.035917466838163, 1.0304003227073102, 1.0408190412167613, 1.0515800164115223, 1.0451743443835133, 1.0553495682011131, 1.0633875811303248, 1.0678851691077176, 1.0746835921427589, 1.077479258815244, 1.0785699503082506, 1.0860063074363802, 1.0947717640914123, 1.0972314840209283, 1.1020794599916695, 1.1050216233619594, 1.1105846565585236, 1.1201265105918397, 1.1208782898796967, 1.1229893203149719, 1.1323989516601276, 1.1346107019686003, 1.1340354936061519, 1.1348144224277317, 1.1428023929216378, 1.1397068522972917, 1.1525317361874394, 1.1586095826442162, 1.1529862368610901, 1.1694123350287862, 1.1726092499582705, 1.1765386716339534, 1.1822346424591577, 1.1811794185312019, 1.181348944384734, 1.1820462635214626, 1.1878143615664276, 1.1932301284637525, 1.1963417622789811, 1.1965433508077516, 1.1949360506763376, 1.2016108100181713, 1.2050177225135812, 1.2087391872784237, 1.2102511200299564, 1.2157645077838495, 1.2174033510409745, 1.2176092176718902, 1.2187931654447122, 1.2245807814810066, 1.2258703244169562, 1.2289263649784743, 1.2282338444746439, 1.2251554441930192, 1.2343310439949184, 1.2308911602195025, 1.2321551598003002, 1.2373927257285455, 1.2358233357272843, 1.2405591132966667, 1.2335252915164836, 1.2404421145359912, 1.2383813945072197, 1.2370798579400355, 1.2422421703435453, 1.2344230247772516, 1.244095383536602, 1.237140315449154, 1.2456924567728742, 1.2427391388955427, 1.2355722395386959, 1.2382849456218854, 1.2351827130442865, 1.2338449730252148, 1.2440538681627811, 1.24073419278467, 1.2353773543938118, 1.2390751290088327, 1.2334264839373366, 1.2363757905285582, 1.2353365877779741, 1.2384207374615757, 1.2326620430452822, 1.2373347610672614, 1.2300156326687948, 1.2408156772490655, 1.2335278111735359, 1.2306980782371399, 1.2304129081001516, 1.2382987142798414, 1.228936293356883, 1.2318124090106941, 1.23664358574744, 1.2348986086441143, 1.236957829033464, 1.2323590724510813, 1.2310151762078057, 1.2286736565131713, 1.2364356837160615, 1.2289582403712302, 1.2323001760610905, 1.2338076460278271, 1.2359435931556384, 1.2307257918696441, 1.229295529555734, 1.2350691974577328, 1.2339992555475714, 1.2384280096760047, 1.2334041778845899, 1.23541547989567, 1.2348257596304149, 1.2406963312950097, 1.2377922038779958, 1.2380676255418481, 1.2402060370824242, 1.2404360896850768, 1.237524725884328, 1.2481373569328822, 1.2426702992121821, 1.2415755951934397, 1.2497583314162575, 1.252809034027178, 1.2545540343573067, 1.2573087909476119, 1.2595374272737783, 1.2622363045460179, 1.2601145575874122, 1.2573740794798007, 1.2684693759582448, 1.2698246498386605, 1.2691346983202032, 1.2698725077189954, 1.2757928590568108, 1.279026277551556, 1.2793626540450607, 1.2832043377384168, 1.2858276489312173, 1.2868888702902783, 1.291497660540377, 1.2949088262843669, 1.2947384882993513, 1.30136652949759, 1.3080668958913935, 1.3077054579480591, 1.312723308754365, 1.3072647345055879, 1.3163188221283579, 1.3268339203124624, 1.3192061100947021, 1.3264339026477288, 1.3251345450087022, 1.3276263969622835, 1.3358928298389756, 1.3317265219855066};
double Us_1[] = {1.000000,1.000556,1.001755,1.003448,1.005462,1.007669,1.009986,1.012363,1.014773,1.017199,1.019632,1.022067,1.024503,1.026940,1.029377,1.031816,1.034256,1.036700,1.039149,1.041603,1.044063,1.046531,1.049008,1.051493,1.053988,1.056493,1.059010,1.061538,1.064078,1.066630,1.069196,1.071774,1.074366,1.076972,1.079592,1.082226,1.084875,1.087538,1.090216,1.092909,1.095618,1.098341,1.101079,1.103833,1.106602,1.109387,1.112187,1.115003,1.117835,1.120682,1.123544,1.126422,1.129316,1.132226,1.135151,1.138091,1.141048,1.144020,1.147007,1.150010,1.153029,1.156063,1.159113,1.162178,1.165258,1.168354,1.171465,1.174592,1.177733,1.180890,1.184063,1.187250,1.190452,1.193670,1.196902,1.200149,1.203412,1.206689,1.209981,1.213288,1.216609,1.219914,1.223166,1.226340,1.229419,1.232397,1.235269,1.238036,1.240700,1.243262,1.245726,1.248094,1.250370,1.252556,1.254656,1.256672,1.258607,1.260463,1.262242,1.263947,1.265580,1.267143,1.268637,1.270064,1.271427,1.272726,1.273964,1.275141,1.276259,1.277319,1.278323,1.279272,1.280166,1.281007,1.281796,1.282534,1.283222,1.283861,1.284452,1.284995,1.285491,1.285941,1.286347,1.286708,1.287025,1.287300,1.287532,1.287722,1.287871,1.287980,1.288049,1.288078,1.288069,1.288021,1.287936,1.287813,1.287654,1.287458,1.287226,1.286959,1.286657,1.286321,1.285950,1.285545,1.285107,1.284636,1.284133,1.283597,1.283029,1.282430,1.281799,1.281137,1.280445,1.279723,1.278970,1.278188,1.277376,1.276535,1.275666,1.274768,1.273841};

int data[] = {
    PJ_DIRECT, PJ_PARABOLIC_PROBLEM, 1, 0,

    FUNCTIOAL_SQUARE, 2, 10, 2, 0,

    // ������� 0 < t <= 0.5
    100, 1, // 0) ��������� �������� (v[0] - ��������������� �������� ��������) 0 < t <= 0.5
    20, 0,   // 1) ������������� ������� (x = 0)
    20, 0,   // 2) ����������� �������
    (N_SENSORS/2), 1,   // 3) ��������� ��������
    20, 0,  // 4) ��������� ������� (t = 0)

    100, 0,  // 5) ������������

    // ������� 0.5 < t <= 1
    100, 1, // 6) ��������� �������� (v[0] - ��������������� �������� ��������) 0 < t <= 0.5
    20, 0,   // 7) ������������� ������� (x = 0) 0 < t <= 0.5
    20, 0,   // 8) ����������� ������� 0 < t <= 0.5
    (N_SENSORS/2), 1,   // 9) ��������� ��������

    NORMALIZED | RADIAL | GAUSSIAN, 32, 2, 0,     2,  1,  1,  0,  0,   0,  -1, -1, -1, -1,
    NORMALIZED | RADIAL | GAUSSIAN, 32, 2, 0,    -1, -1, -1, -1, -1,   0,   2,  1,  1,  0,
    NORMALIZED | RADIAL | GAUSSIAN, 8, 1, 1, 1,  -1, -1, -1,  0, -1,  -1,  -1, -1,  0, -1,
};

#define X_m 10.0
#define T_m 1.0
#define U_m  1.0

#define TAU (T_m/N_SENSORS)

#define X_1 (X_m*0.9)

#define C_alpha ( (T_m/(X_m*X_m)) * 202.33576642336 )
#define C_beta  ( (U_m/X_m) * 0.24912184549463 )

double bounds[] = {
    GP_UNIFORM_SOLID_CUBOID, 0.0, X_m,  TAU, 0.5*T_m,          // 0) ��������� ��������
    GP_UNIFORM_SOLID_CUBOID, 0.0, 0.0,  TAU, 0.5*T_m,          // 1) ������������� �������
    GP_UNIFORM_SOLID_CUBOID, X_m, X_m,  TAU, 0.5*T_m,  // 2) ����������� �������
    GP_GRID_2D, 1, X_1, X_1, (N_SENSORS/2), TAU, 0.5*T_m,                 // 3) �������
    GP_UNIFORM_SOLID_CUBOID, 0.0, X_m,  0.0, 0.0,              // 4) ��������� �������

    GP_UNIFORM_SOLID_CUBOID, 0.0, X_m,  0.5*T_m, 0.5*T_m,      // 5) ������������

    GP_UNIFORM_SOLID_CUBOID, 0.0, X_m,  0.5*T_m + TAU, T_m,    // 6) ��������� ��������
    GP_UNIFORM_SOLID_CUBOID, 0.0, 0.0,  0.5*T_m + TAU, T_m,    // 7) ������������� �������
    GP_UNIFORM_SOLID_CUBOID, X_m, X_m,  0.5*T_m + TAU, T_m,    // 8) ����������� �������
    GP_GRID_2D, 1, X_1, X_1, (N_SENSORS/2), 0.5*T_m + TAU, T_m,           // 9) �������
};

ANNS_FUN(F_q) {double t = x[1]/T_m; return C_beta * ((t <= 0.5) ? (3.0 + 14.0*t) : (10 - 14.0*(t-0.5)));}
//ANNS_FUN(F_q) {double t = x[1]/T_m; return C_beta * (7.0 + 3.0*sin(2.5*M_PI*t));}
//ANNS_FUN(F_q) {double t = x[1]/T_m; return C_beta * ((t <= 0.5) ? 10.0 : 5.0);}

ANNS_FUN(F_u) {double x_ = x[0]/X_m; return C_alpha * x_ * (1-x_);}

ANNS_CND(Ctransfer_1)  {return u_xx[0][0] - v[0]*u_x[0][1];}    // 0
ANNS_CND(Cinsulated_1) {return u_x[0][0];}                      // 1
ANNS_CND(Cheatflux_1)  {return u_x[0][0] - u[2];}               // 2
ANNS_CND(Csensors_1) {return u[0] - v[0];}                      // 3
ANNS_CND(Cinitial)   {return u[0] - U_m;}                       // 4

ANNS_CND(Cconsistency) {return u[0] - u[1];}                    // 5

ANNS_CND(Ctransfer_2)  {return u_xx[1][0] - v[0]*u_x[1][1];}    // 6
ANNS_CND(Cinsulated_2) {return u_x[1][0];}                      // 7
ANNS_CND(Cheatflux_2)  {return u_x[1][0] - u[2];}               // 8
ANNS_CND(Csensors_2) {return u[1] - v[0];}                      // 9


ANNS_CND_G(Ctransfer_a)  {return u_axx[0] - v[0]*u_ax[1];}
ANNS_CND_G(Cinsulated_a) {return u_ax[0];}
ANNS_CND_G(Cheatflux_1_a)  {return u_ax[0];}
ANNS_CND_G(Cheatflux_2_a)  {return -u_a;}
ANNS_CND_G(Csensors_a) {return u_a;}
ANNS_CND_G(Cinitial_a)   {return u_a;}

ANNS_CND_G(Cconsistency_1) {return u_a;}
ANNS_CND_G(Cconsistency_2) {return -u_a;}


fbvp_t cond[] = {
    Ctransfer_1,     // 0
    Cinsulated_1,  // 1
    Cheatflux_1,   // 2
    Csensors_1,    // 3
    Cinitial,      // 4

    Cconsistency,  // 5

    Ctransfer_2,   // 6
    Cinsulated_2,  // 7
    Cheatflux_2,   // 8
    Csensors_2,    // 9
};

gbvp_t cond_a[] = {
    Ctransfer_a,    // 0
    Cinsulated_a,   // 1
    Cheatflux_1_a,  // 2
    Csensors_a,     // 3
    Cinitial_a,     // 4
    Cconsistency_1, // 5
    NULL,           // 6,7,8,9
    NULL,
    NULL,
    NULL,


    NULL,           // 1,2,3,4
    NULL,
    NULL,
    NULL,
    Cconsistency_2, // 5
    Ctransfer_a,    // 6
    Cinsulated_a,   // 7
    Cheatflux_1_a,  // 8
    Csensors_a,     // 9


    NULL,          // 0,1
    NULL,
    Cheatflux_2_a, // 2
    NULL,          // 3,4,5,6,7
    NULL,
    NULL,
    NULL,
    NULL,
    Cheatflux_2_a, // 8
    NULL,          // 9
};

double delta[] = {
    1.0,   // 0
    1.0,   // 1
    1.0,   // 2
    10.0,  // 3
    100.0, // 4
    1.0,   // 5
    1.0,   // 6
    1.0,   // 7
    1.0,   // 8
    10.0,  // 9
};

void ai_init(anns_instance_t *ai) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

void ai_reload(anns_instance_t *ai) {
    ai_genpoints(ai, bounds);
    ai_eval_v(ai, 0, 0, F_u, NULL);
    ai_load_v(ai, 3, 0, Us_1);
    ai_eval_v(ai, 6, 0, F_u, NULL);
    ai_load_v(ai, 9, 0, Us_1+80);
}

void ai_initann(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;

    for (i = ai->I[0]; i < ai->I[1]; ) {
        x[i++] = U_m;               // ����
        x[i++] = uniform(0.05 * T_m, 1.05 * T_m);  // ������
        x[i++] = uniform(0.75 * X_m, 1.05 * X_m); // ������
        x[i++] = uniform(-0.05 * T_m, 0.5 * T_m); // ������
    }

    for (i = ai->I[1]; i < ai->I[2]; ) {
        x[i++] = U_m;               // ����
        x[i++] = uniform(0.05 * T_m, 1.05 * T_m);  // ������
        x[i++] = uniform(0.75 * X_m, 1.05 * X_m); // ������
        x[i++] = uniform(0.5 * T_m , 1.05 * T_m); // ������
    }

    for (i = ai->I[2]; i < ai->I[3]; ) {
        x[i++] = 0.25;               // ����
        x[i++] = uniform(0.05 * T_m, 0.25 * T_m);  // ������
        x[i++] = uniform(-0.05 * T_m, 1.05 * T_m); // ������
    }
}

double ai_mse(anns_instance_t *bi, double *x) {
    return 0.0;
}

#endif // ANNS_TASK
