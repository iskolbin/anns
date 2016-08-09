#ifndef ANNS_INSTANCE_H_INCLUDED
#define ANNS_INSTANCE_H_INCLUDED

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

//#include "anns_net.h"
//#include "anns_cond.h"
#include "anns_cache.h"
#include "anns_types.h"
#include "anns_const.h"


anns_instance_t *ai_new(int data[], void *instance);
void ai_delete(anns_instance_t *ai);
void ai_conditions(anns_instance_t *ai, fbvp_t *fbvp, gbvp_t *gbvp, double *delta);
void ai_updatenets(anns_instance_t *ai, double *newz, int *nc);

#define ai_eval_v(ai,c,iv,f,inst) (ac_eval_v((ai->cache),(c),(iv),(f),(inst)))
#define ai_load_v(ai,c,iv,darr) (ac_load_v((ai->cache),(c),(iv),(darr)))
#define ai_set_v(ai,c,iv,val) (ac_set_v((ai->cache),(c),(iv),(val)))
#define ai_push(ai,x) (ac_put((ai)->cache, (x)))
#define ai_put(ai,x,i) (ac_put((ai)->cache,(x),(i)))
#define ai_get(ai,i) (ac_get((ai)->cache,(i)))
#define ai_genpoints(ai,bounds) (ac_genpoints((ai)->cache, (bounds)))
//#define ai_savez(ai) {int i; for (i = 0; i < ai->A; i++) ai->zbest[i] = ai->z[i];}
#define ai_ischeckpoints(ai) (ac_ischeckpoints((ai)->cache))
#define ai_switchpoints(ai)  (ac_switchpoints((ai)->cache))
#define ai_isbestz(ai)       ((ai)->z == (ai)->zbest)

void ai_refine(anns_instance_t *ai);
void ai_fprintf(anns_instance_t *ai, FILE *fout, char *valf, char *valsep, char *nsep, int aux);
void ai_savez(anns_instance_t *ai);
void ai_switchz(anns_instance_t *ai);
double *ai_collectweights(anns_instance_t *ai);
double *ai_expandweights(anns_instance_t *ai);

#define ai_printf(ai) (ai_fprintf(ai, stdout, "%g", " ", "\n", 1))

#endif // ANNS_INSTANCE_H_INCLUDED
