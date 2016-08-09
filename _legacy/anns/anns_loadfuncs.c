#include "anns_loadfuncs.h"

void ai_loadfuncs(anns_instance_t *ai) {
    int t;

    assert(ai);

    for (t = 0; t < ai->nnets; t++) {
        switch(ai->nets[t]->T) {
            case NORMALIZED | RADIAL | GAUSSIAN: nrbf_gauss_loadfuncs(ai->nets[t]); break;
//            case NORMALIZED | RADIAL | GAUSSIAN: nrbf_gauss_loadfuncs(t, ai); break;
            case CLASSIC | RADIAL | GAUSSIAN: rbf_gauss_loadfuncs(ai->nets[t]); break;
//            case NORMALIZED | ELLIPTIC | GAUSSIAN: nebf_gauss_loadfuncs(t, ai); break;
//            case CLASSIC | ELLIPTIC | GAUSSIAN: ebf_gauss_loadfuncs(t, ai); break;
//
//            case NORMALIZED | RADIAL | CAUCHY_FUNCTION: nrbf_cauchy_loadfuncs(t, ai); break;
//            case CLASSIC | RADIAL | CAUCHY_FUNCTION: rbf_cauchy_loadfuncs(t, ai); break;
        }
    }
}
