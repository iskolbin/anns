#include "anns_tablecalc.h"

anns_tablecalc_t *at_new(double from, double to, double step, double (*func)(double)) {
    anns_tablecalc_t *at = malloc(sizeof *at);
    int len = (int)((to-from)/step) + 2;
    int i;
    double f1, f2, x1, x2, k;

    assert(at);
    assert(len > 0);

    at->from = from;
    at->to = to;
    at->step = step;
    at->invstep = 1. / step;
    at->func = func;
    at->len = len;
    at->k = calloc(len, sizeof *at->k);
    at->b = calloc(len, sizeof *at->b);

    assert(at->k);
    assert(at->b);

    x2 = from;
    f2 = func(x2);
    for (i = 0; i < len; i++) {
        x1 = x2;
        x2 += step;
        f1 = f2;
        f2 = func(x2);
        k = (f2-f1) / (x2-x1);
        at->k[i] = k;
        at->b[i] = f1 - k * x1;
    }

    return at;
}

void at_delete(anns_tablecalc_t *at) {
    if (at) {
        free(at->k);
        free(at->b);
        free(at);
    }
}

void at_fprintf(anns_tablecalc_t *at, FILE *fout, char *format) {
    int i;

    assert(at);
    assert(fout);
    assert(format);

    fprintf(fout, "[%g,%g] step=%g len=%d\n", at->from, at->to, at->step, at->len);
    fprintf(fout, "k = {");
    for (i = 0; i < at->len; i++) {
        fprintf(fout, format, at->k[i]);
    }
    fprintf(fout, "}\nb = {");
    for (i = 0; i < at->len; i++) {
        fprintf(fout, format, at->b[i]);
    }
    fprintf(fout, "}");
}

double at_eval(anns_tablecalc_t *at, double x) {
    int index;

//    assert(0);
    assert(at);
    assert(x >= at->from);
    assert(x <= at->to);

    index = (int)((x-at->from) * at->invstep);
    return at->k[index] * x + at->b[index];
}
