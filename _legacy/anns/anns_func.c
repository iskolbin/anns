#include "anns_func.h"

//anns_arrfuncs_t *aa_new(int n) {
//    anns_arrfuncs_t *aa = NULL;
//
//    assert(n > 0);
//
//    aa = malloc(sizeof *aa);
//
//    assert(aa);
//
//    aa->n = n;
//    aa->funcs = calloc(n, sizeof *aa->funcs);
//
//    assert(aa->funcs);
//
//    return aa;
//}
//
//anns_arrfuncs_t *aa_fromlist(int n, ...) {
//    va_list args;
//    anns_arrfuncs_t *aa = NULL;
//    int i;
//
//    aa = af_new(n);
//
//    va_start(args, n);
//    for (i = 0; i < n; i ++) {
//        aa->funcs[i] = va_arg(args, anns_func_t);
//    }
//    va_end(args);
//}
//
//void aa_delete(anns_arrfuncs_t *aa) {
//    assert(aa);
//
//    free(af->funcs);
//    free(af);
//}
