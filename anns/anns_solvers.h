#ifndef ANNS_SOLVERS_H_INCLUDED
#define ANNS_SOLVERS_H_INCLUDED

#include <time.h>
#include "anns.h"
#include "anns_cache.h"
#include "../optimize/fmin_cg_descent.h"

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

    void (*ai_init)(anns_instance_t *, fbvp_t *, gbvp_t *, double *);
    void (*ai_update)(anns_instance_t *, double *, int );
    void (*ai_newpoints)(anns_instance_t *, double *);
    void (*ai_initann)(anns_instance_t *, double *);
    double (*ai_mse)(anns_instance_t *, double *);

    double rho;
    double lambda;

    int fdsteps;  // ����� ����� ��� �������-���������� �����
    double h;   // �������� ����

    int n0;           // ��������� ����� ��������
    int nmax;         // ������������ ����� ��������

    int maxTrials;

    int Nopt;              // ������������ ����� ����� �����������
    int ticksGrow;         // ����� ������ ����� ������ ����
    int ticksReduce;       // ����� ������ ����� ��������� ����
    int ticksRegeneration; // ����� ������ ����� �������������� �����
    int ticksTrial;        // ����� ������ ����� ������������
    int ticksMax;          // ������������ ����� ������
} anns_solverdata_t;

anns_solverdata_t *as_new(void);
anns_solverdata_t *as_new2(void (*ai_init)(anns_instance_t *, fbvp_t *, gbvp_t *, double *),
                            void (*ai_newpoints)(anns_instance_t *, double *),
                            void (*ai_initann)(anns_instance_t *, double *),
                            double (*ai_mse)(anns_instance_t *, double *),
                            int *data, double *bounds, fbvp_t *cond, gbvp_t *cond_a, double *delta);

#define as_delete(as) {if (as) free(as);}

int anns_solve_restarts(anns_solverdata_t *as);

#endif // ANNS_SOLVERS_H_INCLUDED
