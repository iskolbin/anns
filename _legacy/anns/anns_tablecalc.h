#ifndef ANNS_TABLECALC_H_INCLUDED
#define ANNS_TABLECALC_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

typedef struct {
    int len;
    double from;
    double to;
    double step;
    double invstep;
    double *k;
    double *b;
    double (*func)(double);

    int tmpindex;
    //    double tmpx;
} anns_tablecalc_t;

anns_tablecalc_t *at_new(double from, double to, double step, double (*func)(double));
void at_delete(anns_tablecalc_t *at);
void at_fprintf(anns_tablecalc_t *at, FILE *fout, char *format);
double at_eval(anns_tablecalc_t *at, double x);

#define at_size(at)          ((at)->len*(sizeof(*(at)->k)+sizeof(*(at)->b)))
#define at_idx(at,x)         ((int)(((x)-(at)->from) * ((at)->invstep)))
//#define at_fasteval(at,x)    ((at)->k[(at)->tmpindex = at_idx((at),((at)->tmpx = (x)))]*((at)->tmpx) + (at)->b[(at)->tmpindex])
#define at_fasteval(at,x)    ((at)->k[(at)->tmpindex = at_idx((at),(x))]*(x) + (at)->b[(at)->tmpindex])
#define at_printf(at,format) (at_fprintf((at),stdout,(format)))

#endif // ANNS_TABLECALC_H_INCLUDED
