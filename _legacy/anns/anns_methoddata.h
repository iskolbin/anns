#ifndef ANNS_METHODDATA_H_INCLUDED
#define ANNS_METHODDATA_H_INCLUDED

#include <time.h>
//#include <stdio.h>
//#include "anns_types.h"
#include "anns.h"
#include "../optimize/fmin_callback.h"

typedef struct {
    int nofileout;
    int randomseed;
    long int seed;
    double tol;
    double gtol;

    int *data;
    double *bounds;
    fbvp_t *cond;
    gbvp_t *cond_a;
    double *delta;

    fbvp_t *cond_i;
    gbvp_t *cond_i_a;
    double tau;
    double T0;
    double Tm;
    int nsteps;

    void (*ai_init)(anns_instance_t *, fbvp_t *, gbvp_t *, double *);
    void (*ai_update)(anns_instance_t *, double *, int ,double );
    void (*ai_newpoints)(anns_instance_t *, double *);
    void (*ai_initann)(anns_instance_t *, double *);
    double (*ai_mse)( double *, int, void *);

    double rho;
    double lambda;

    int fdsteps;  // Число шагов для конечно-разностной схемы
    double h;   // Величина шага

    int n0;           // Начальное число нейронов
    int nmax;         // Максимальное число нейронов

    int maxTrials;
    int maxIters;

    int Nopt;              // Максимальное число шагов оптимизации
    int ticksGrow;         // Число тактов перед ростом сети
    int ticksReduce;       // Число тактов перед редукцией сети
    int ticksRegeneration; // Число тактов перед перегенерацией точек
    int ticksTrial;        // Число тактов перед перезапуском
    int ticksMax;          // Максимальное число тактов

    void (*preeval)(double *z, void *anns_instance);
    void (*updatepoints)(double *z, void *anns_instance);
} anns_methoddata_t;

anns_methoddata_t *am_new(void);

anns_methoddata_t *am_new2(void (*ai_init)(anns_instance_t *, fbvp_t *, gbvp_t *, double *),
                            void (*ai_newpoints)(anns_instance_t *, double *),
                            void (*ai_initann)(anns_instance_t *, double *),
                            double (*ai_mse)(double *, int, void *),
                            void (*ai_update)(anns_instance_t *, double *, int ,double ),
                            int *data, double *bounds,
                            fbvp_t *cond, gbvp_t *cond_a, double *delta);

anns_methoddata_t *am_new3(void (*ai_init)(anns_instance_t *, fbvp_t *, gbvp_t *, double *),
                            void (*ai_newpoints)(anns_instance_t *, double *),
                            void (*ai_initann)(anns_instance_t *, double *),
                            double (*ai_mse)(double *, int, void *),
                            void (*ai_update)(anns_instance_t *, double *, int , double ),
                            int *data, double *bounds,
                            fbvp_t *cond, gbvp_t *cond_a, double *delta,
                            fbvp_t *cond_i, gbvp_t *cond_i_a, double T0, double Tm, int nsteps);

typedef struct {
    int key;
    void *value;
} anns_methoddata_pointer_t;

typedef struct {
    int key;
    double value;
} anns_methoddata_number_t;

anns_methoddata_t *am_new_n(anns_methoddata_pointer_t* pointers, anns_methoddata_number_t* numbers);

#define am_delete(am) {if (am) free(am);}

#endif // ANNS_METHODDATA_H_INCLUDED
