#ifndef ANNS_TASK
#define ANNS_TASK

#define PROBLEM_SENSORS 32
#define G_M 100

//exact
//double usensors[] = {0.000000,0.030969,0.061226,0.090192,0.117406,0.142513,0.165248,0.185418,0.202897,0.217612,0.229532,0.238668,0.245059,0.248768,0.249883,0.248505,0.244751,0.238746,0.230625,0.220528,0.208601,0.194993,0.179855,0.163338,0.145596,0.126783,0.107052,0.086558,0.065452,0.043889,0.022021,0.000000};
////0.001
//double usensors[] = {0.000000,0.028565,0.060877,0.089527,0.118056,0.142083,0.164743,0.185504,0.203039,0.216634,0.229958,0.238054,0.246073,0.248544,0.248799,0.247449,0.245596,0.239774,0.229936,0.219963,0.208880,0.192961,0.180141,0.163350,0.145468,0.127073,0.107613,0.085385,0.065385,0.042729,0.021287,0.000000};
////0.01
//double usensors[] = {0.000000,0.006935,0.057734,0.083544,0.123908,0.138213,0.160204,0.186282,0.204312,0.207841,0.233788,0.232525,0.255205,0.246523,0.239042,0.237945,0.253198,0.249027,0.223734,0.214877,0.211387,0.174673,0.182714,0.163460,0.144315,0.129685,0.112662,0.074832,0.064781,0.032294,0.014680,0.000000};
////0.1
//double usensors[] = {0.000000,-0.209363,0.026303,0.023714,0.182426,0.099513,0.114811,0.194060,0.217044,0.119901,0.272088,0.177237,0.346525,0.226318,0.141476,0.142905,0.329221,0.341560,0.161722,0.164019,0.236453,-0.008211,0.208443,0.164557,0.132786,0.155805,0.163144,-0.030695,0.058743,-0.072064,-0.051391,0.000000};
// ������
double usensors[] = {0,0.0309685390250146,0.0612264960007404,0.0901922036669068,0.11740621726668,0.142513455227179,0.165247670402037,0.185417970154539,0.202897136293787,0.217611524915719,0.229532351952742,0.238668193062393,0.245058546710018,0.248768327210078,0.24988317034153,0.248505448172131,0.244750902116977,0.238745814197909,0.23062464612205,0.2205280843015,0.208601436417187,0.194993331699245,0.179854682851749,0.163337872577566,0.145596132035115,0.126783082349203,0.107052413560524,0.0865576781826515,0.0654521788841493,0.0438889317629586,0.0220206882606651,0.000000000000000};
// 0.001
//double usensors[] = {0.000000000000000,0.031520011885503,0.060123261095788,0.091308566157251,0.117257684709312,0.144127390849381,0.164729865690680,0.184317892996862,0.202118581966265,0.217345426089344,0.230418602363026,0.238623373735714,0.245678896747534,0.249019899224309,0.250228026878119,0.247951618903617,0.243674805990902,0.238755734744094,0.230164771731013,0.219571083508054,0.209842500771629,0.194436928420838,0.177842353091552,0.162983680313322,0.146512159214164,0.124961105526915,0.107542101820890,0.085018178513787,0.064837322810122,0.045078583172488,0.020291968655025,0.000000000000000};
// 0.01
//double usensors[] = {0.000000000000000,0.036483267629900,0.050194146951221,0.101355828570344,0.115920891693003,0.158652811449195,0.160069623288467,0.174417198577767,0.195111593018566,0.214950536651977,0.238394856055581,0.238219999795599,0.251262047085178,0.251284047352382,0.253331735707421,0.242967155486989,0.233989940856227,0.238845019659759,0.226025902211684,0.210958076367043,0.221012079961601,0.189429298915169,0.159731385249787,0.159795949935123,0.154756403825597,0.108563314126318,0.111949296164191,0.071162681494010,0.059303618143877,0.055785445858251,0.004733492204262,0.000000000000000};
// 0.1
//double usensors[] = {0.000000000000000,0.086115825073868,-0.049096994494458,0.201828452701281,0.102552961529918,0.303907017447334,0.113467199266333,0.075410254386816,0.125041703541579,0.191001642278299,0.318157392981131,0.234186260394450,0.307093550461613,0.273925528633113,0.284368824000433,0.193122521320714,0.137141289509480,0.239737868816411,0.184637207018389,0.124828004956930,0.332707871861327,0.139353003858482,-0.021378293167868,0.127918646153137,0.237198849939931,-0.055414599879649,0.156021239597194,-0.067392288703765,0.003966571481427,0.162854072715884,-0.150851272303365,0.000000000000000};


int data[] = {
    PJ_INVERSE_UNKNOWN_SOURCE, PJ_ELLIPTIC_PROBLEM, 4, 2,

    FUNCTIONAL_SQUARE, 2, 5, 2, 0,

    256, 1,   // ��������� ��������
    240, 0,   // x = 0
    240, 0,   // x = pi
    220, 0,   // y = 0
    PROBLEM_SENSORS, 1,   // y = 1
//        2, 0,
//        1, 1,

    NORMALIZED | RADIAL | GAUSSIAN, 24, 2, 0,   2, 0,  0,  0,  0,
//     -1, 2,
    NORMALIZED | RADIAL | GAUSSIAN, 8,  1, 0,   0, -1, -1, -1, -1,
//    0, -1,
};

double bounds[] = {
//    GP_UNIFORM_SOLID_CUBOID, 0.0, M_PI,   0.0, 2.5,  // ��������� ��������
    GP_GRID_2D,              16, 0.0, M_PI,   16, 0.0, 5.5,  // y = 1
    GP_UNIFORM_SOLID_CUBOID, 0.0, 0.0,    0.0, 5.5,  // x = 0
    GP_UNIFORM_SOLID_CUBOID, M_PI, M_PI,  0.0, 5.5,  // x = pi

    GP_UNIFORM_SOLID_CUBOID, 0.0, M_PI,   0.0, 0.0,        // y = 0
    GP_GRID_2D,              PROBLEM_SENSORS, 0.0, M_PI,   1, 1.0, 1.0,  // y = 1
//        GP_POINTS, 0.0, 0.0, M_PI, 0.0,

//     GP_UNIFORM_SOLID_CUBOID, 0.0, M_PI,   0.0, 1.650,  // ��������� ��������

};

ANNS_CND(A)  {return u_xx[0][0] + u_xx[0][1] + u[1];}
ANNS_CND(B) {return u[0];}
ANNS_CND(Btop)  {return u[0] - v[0];}
ANNS_CND(C) {return u[1];}
ANNS_CND(D) {return u_xx[0][0] + u_xx[0][1] + v[0];}

ANNS_CND_G(A_a_1)  {return u_axx[0] + u_axx[1];}
ANNS_CND_G(A_a_2)  {return u_a;}
ANNS_CND_G(B_a) {return u_a;}
ANNS_CND_G(C_a) {return u_a;}

ANNS_FUN(g) {int i;
    double s = 0;
    for (i = 1; i < G_M; i++) {
        s += sin(i*x[0]) / ( (1+i*i) * (1+i*i) * i);
    }
    return s;}

ANNS_FUN(f) {int i;
    double s = 0;
    for (i = 1; i < G_M; i++) {
        s += sin(i*x[0]) / ( (1+i*i) * (1+i*i) * i);
    }
    return s;
}

fbvp2_t cond[] = {
    A,
    B, B, B, Btop,
//    C,
//    D,
    };

gbvp2_t cond_a[] = {
     A_a_1,
     B_a, B_a, B_a, B_a,
//     NULL, A_a_1,

    A_a_2, NULL, NULL, NULL, NULL,
//    C_a, NULL
    };

double delta[] = {
    1.0, 10.0, 10.0, 10.0, 100.0, 1.0, 1.0};

void ai_init(anns_instance_t *ai) {
    ai_conditions(ai, cond, cond_a, delta);
    ai_loadfuncs(ai);
}

void ai_reload(anns_instance_t *ai) {
    ai_genpoints(ai, bounds);
    ai_load_v(ai, 4, 0, usensors);
//    ai_eval_v(ai, 6, 0, f, NULL);
}

void ai_initann(anns_instance_t *ai, double *x) {
    int i, j, t, n, s;

    for (i = ai->I[0]; i < ai->I[1]; ) {
        x[i++] = uniform(0.0,0.0);               // ����
        x[i++] = uniform(0.05, 1.05);  // ������
        x[i++] = uniform(-0.05*M_PI, 1.05*M_PI); // ������
        x[i++] = uniform(-0.05, 3.05); // ������
    }
//
    for (i = ai->I[1]; i < ai->I[2]; ) {
        x[i++] = uniform(0.0,0.0);                // ����
        x[i++] = uniform(0.05, 0.5);  // ������
        x[i++] = uniform(-0.05*M_PI, 1.05*M_PI); // ������
    }
}

double ai_mse(anns_instance_t *bi, double *x) {
    double hx = 0.001;
    int i, n = (int)(M_PI/hx);
    double point[] = {0.0, 0.0};
    double mse = 0, tmp;

    for (i = 0; i < n; i++) {
        point[0] = uniform(0.0, M_PI);
        point[1] = uniform(1.0, 1.0);

        bi->ann_val[0][2](point, x, 0, 0, bi);
        tmp = bi->cache->u_xx[0][0] + bi->cache->u_xx[0][1] + f(point, 1, NULL);
        mse += tmp*tmp;
    }

    return sqrt(mse/(n-1));
}

#endif // ANNS_TASK
