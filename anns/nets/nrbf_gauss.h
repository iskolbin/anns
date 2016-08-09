#ifndef ANNS_NETS_NRBF_GAUSS_H_INCLUDED
#define ANNS_NETS_NRBF_GAUSS_H_INCLUDED

// НРБС-Г (NRBF-G)
// NORMALIZED | RADIAL | GAUSSIAN

#include <math.h>
#include <stdio.h>
#include "../anns_types.h"
//#include "../anns_tablecalc.h"

// Вычисление НРБС на гауссианах с нормировкой (тавтология, но тем не менее)
//
// В отличие от обычного вычисления экспонент, тут предварительно вычисляется самый большой
// показатель, и вычитается из всех. С математической точки зрения это просто деление
// числителя и знаменателя на константу, т.е. полученные нормированные выражения будут тождествено
// равны ненормированным. Нормировка значительно повышает точность машинных вычислений, при этом практически
// не меняя алгоритм (само уравнение не меняется). Мотивация:
// 1) Экспонента очень быстро убывает, и даже на не очень значительном
// отдалении от центров обнуляется (из-за ограничений встроенного вещественномго типа)
// Для ненормализованных сетей это не принципиально (выход сети обнулится), а в случае НРБС возникает
// деление 0/0.
// 2) Возможны ситуации, когда рядом стоящиие нейроны: часть обнуляется, а часть просто
// даёт очень малый выход, тогда вычисление отношения может давать совершенно неточные
// результаты.
// 3) Таким образом, вычислять НРБС без нормировки просто не имеет смысла, т.к. дополнительная вычислительная
// нагрузка крайне незначительна (nmax операций сравнения и вычитания), а эффект существенен

//void nrbf_gauss_val(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);
//void nrbf_gauss_val_d1(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);
//void nrbf_gauss_val_d2(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);
//
//void nrbf_gauss_valgrad(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);
//void nrbf_gauss_valgrad_d1(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);
//void nrbf_gauss_valgrad_d2(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);
//
//void nrbf_gauss_valgrad_wonly(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);
//void nrbf_gauss_valgrad_d1_wonly(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);
//void nrbf_gauss_valgrad_d2_wonly(double *point, double *x, int t, int mixed_derivatives, void *anns_instance);

void nrbf_gauss_val(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);
void nrbf_gauss_val_d1(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);
void nrbf_gauss_val_d2(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);

void nrbf_gauss_valgrad(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);
void nrbf_gauss_valgrad_d1(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);
void nrbf_gauss_valgrad_d2(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);

void nrbf_gauss_val_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);
void nrbf_gauss_val_d1_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);
void nrbf_gauss_val_d2_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);

void nrbf_gauss_valgrad_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);
void nrbf_gauss_valgrad_d1_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);
void nrbf_gauss_valgrad_d2_wonly(double *x, int c, int p, int t, int mixed_derivatives, void *anns_instance);



//void nrbf_gauss_loadfuncs(int t, void *anns_instance);
void nrbf_gauss_loadfuncs(anns_net_t *net);

#endif // ANNS_NETS_NRBF_GAUSS_H_INCLUDED
