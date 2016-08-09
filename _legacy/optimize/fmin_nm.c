#include "fmin_nm.h"

fmin_nm_args_t fmin_nm_default_args = {
    1.0, 2.0, -0.5, 0.5,
    0xFFFF, NULL, -1.0, 1.0, NULL,
    FMIN_NM_CONV_SIMPLEX_SIZE
};

// A Simplex Method for Function Minimization
// J. A. Nelder and R. Mead
// The Computer Journal (1965) 7 (4): 308-313.
// doi: 10.1093/comjnl/7.4.308
double *fmin_nm(double *x,
               int n,
               double (*val)(double *, int, void *),
               double tol,
               void *instance,

               fmin_nm_args_t *args,
               fmin_nm_stats_t *stats,
               double *work_array) {

    void *preallocated = work_array;
    int iw, is, ib;             // ������� ������, ���������� � ������ �������
    double *xw, *xs, *xb, *xo;  // ������, ����������, ������, ����������� �������
    double *xr, *xe, *xc;       // ����������, ��������� � ������ �������
    double *xt;
    double fw, fs, fb, fo;      // �������� ������� � ������, ����������, ������, ����������� �������
    double fr, fe, fc;          // �������� ������� � ����������, ���������� � ������ �������
    double *simplex;            // ������� ��������� ���������� � ���� ������ � ����� n
    double *f;                  // �������� ������� � ��������
    double a, b, s;
    int current_iteration, max_iterations;
    double reflection, expansion, contraction, reduction, current_tol;
    int i, j;
    FILE *log = NULL;

    // ���� ������ �� ���� ����� - ���������� ���� (����� ���� n * (n + 8) + 1)
    if (!preallocated) {
        work_array = calloc(n * (n + 8) + 1, sizeof *work_array);
    }

    // ���� ���� ���������� ��������� - ��������� �� (����� ���������� �� ���������)
    args = args ? args : &fmin_nm_default_args;
    reflection = args->reflection;
    expansion = args->expansion;
    contraction = args->contraction;
    reduction = args->reduction;
    max_iterations = args->max_iterations;
    a = args->simplex_min;
    b = args->simplex_max;
    log = args->log;

    assert(work_array);

    // ������������ ������
    simplex = work_array;
    xo = simplex + n*(n+1);
    xw = xo + n;
    xs = xw + n;
    xb = xs + n;
    xr = xb + n;
    xe = xr + n;
    xc = xe;
    f = xc + n;

    if (log) {
        fprintf(log, "[Begin Nelder-Mead]\n[simplex of %d vertices]\n", n+1);
        if (!preallocated) fprintf(log, "[memory allocated %d * %d]", n*(n+8)+1, sizeof(double));
    }

    // �������������� �������� - ��������� ������������ ����������� ����� ��� ������ �� ����������������
    for (i = 0; i < n + 1; i++) {
        for (j = 0; j < n; j++) {
            simplex[i*n + j] = args->init ? args->init[i*n + j] : x[j] +  rand()/(RAND_MAX + 1.0)*(b-a) + a;
        }
    }

    // ������� ����
    for (current_iteration = 0; current_iteration < max_iterations; current_iteration++) {
        // ��������� ������� � �������� ���������
        // ���������� ������, �������������� ������� � ������ ������� (�������)
        xw = xs = xb = simplex;
        iw = is = ib = 0;
        fw = fs = fb = fc = f[0] = val(simplex, n, instance);
        for  (i = 1; i < n + 1; i++) {
            f[i] = val(simplex + i*n, n, instance);
            if (f[i] < fb) {
                fb = f[i];
                ib = i;
                xb = simplex + ib * n;
            } else if (f[i] > fw) {
                fs = fw;
                is = iw;
                xs = xw;

                fw = f[i];
                iw = i;
                xw = simplex + iw * n;
            } else if (f[i] > fs) {
                fs = f[i];
                is = i;
                xs = simplex + is * n;
            }
        }

        // ��������� ����� ���� (����� ������ �������)
        for (j = 0; j < n; j++) xo[j] = 0;

        for (i = 0; i < iw; i++) {
            for (j = 0; j < n; j++) xo[j] += simplex[i*n + j];
        }

        for (i = iw+1; i < n+1; i++) {
            for (j = 0; j < n; j++) xo[j] += simplex[i*n + j];
        }

        for (j = 0; j < n; j++) {
            xo[j] /= n;
        }
        fo = val(xo, n, instance);

        // ��������� �������� ����������
        if (args->convergence_criterion == FMIN_NM_CONV_FUNCTION) {
            current_tol = fb;
        } else if (args->convergence_criterion == FMIN_NM_CONV_SIMPLEX_SIZE) {
            s = 0;
            for (i = 0; i < n+1; i++) {
                s += (f[i] - fo)*(f[i] - fo);
            }
            current_tol = sqrt(s/(n+1));
        }

        if (log) fprintf(log, "\n[current tolerance %g]", current_tol);

        if (current_tol <= tol) {
            if (log) fprintf(log, "\n[converged (tol=%g, needed=%g)]", current_tol, tol);
            break;
        }

        // ��������� ���������� �����
        for (j = 0; j < n; j++) xr[j] = xo[j] + reflection * (xo[j] - xw[j]);
        fr = val(xr, n, instance);

        if (log) {
            fprintf(log, "\n\n[Iteration %d]\n[simplex]\n", current_iteration);

            for (i = 0; i < n + 1; i++) {
                fprintf(log, "(%g) ",f[i]);
                for (j = 0; j < n; j++) fprintf(log, "%g ", simplex[i*n + j]);
                fprintf(log, "\n");
            }

            fprintf(log, "[best(%g)] ", fb);
            for (j = 0; j < n; j++) fprintf(log, "%g ", xb[j]);

            fprintf(log, "\n[worst(%g)] ", fw);
            for (j = 0; j < n; j++) fprintf(log, "%g ", xw[j]);

            fprintf(log, "\n[preworst(%g)] ", fs);
            for (j = 0; j < n; j++) fprintf(log, "%g ", xs[j]);

            fprintf(log, "\n[center(%g)] ", fo);
            for (j = 0; j < n; j++) fprintf(log, "%g ", xo[j]);

            fprintf(log, "\n[reflected(%g)] ", fr);
            for (j = 0; j < n; j++) fprintf(log, "%g ", xr[j]);
        }

        // ���� ������ ������� ����� ����������, ���������� ����� ����������
        if (fb <= fr && fr < fs) {
            // ��������� - ���������� � ������ ������� ���������� �����
            if (log) fprintf(log, "\n[reflection]");

            for (j = 0; j < n; j++) xw[j] = xr[j];
            f[iw] = fr;

        // �����, ���� ���������� ����� ������, �� ������� ���� � ��� �� �����������
        } else if (fr < fb) {
            // ��������� ����� ����������
//            for (j = 0; j < n; j++) xe[j] = xo[j] + expansion*(xo[j] - xw[j]);
            for (j = 0; j < n; j++) xe[j] = xo[j] + expansion*(xr[j] - xo[j]);
            fe = val(xe, n, instance);

            if (log) {
                fprintf(log, "\n[expanded(%g)] ", fe);
                for (j = 0; j < n; j++) fprintf(log, "%g ", xe[j]);
            }

            // ���� ����� ���������� ����� ���������� (�.�. ���� ������ ���������)
            if (fe < fr) {
                // ���������� - ���������� � ������ ������� ����� ����������

                if (log) fprintf(log, "\n[expansion]");

                for (j = 0; j < n; j++) xw[j] = xe[j];
                f[iw] = fe;

            // ����� (���������� ����� �� ����, ����� �����������)
            } else {
                // ��������� - ���������� � ������ ������� ���������� �����
                if (log) fprintf(log, "\n[reflection]");

                for (j = 0; j < n; j++) xw[j] = xr[j];
                f[iw] = fr;
            }

        // ����� (���������� �� ����� ����������)
        } else {
            // ��������� ����� ������
            for (j = 0; j < n; j++) xc[j] = xo[j] + contraction*(xo[j] - xw[j]);
            fc = val(xc, n, instance);

            if (log) {
                fprintf(log, "\n[contracted(%g)] ", fc);
                for (j = 0; j < n; j++) fprintf(log, "%g ", xc[j]);

            }

//            if (fc < fr) {
            if (fc < fw && fc < fr) {
                // ������ - ���������� � ������ ������� ����� ������
                if (log) fprintf(log, "\n[contraction]");

                for (j = 0; j < n; j++) xw[j] = xc[j];
                f[iw] = fc;

            } else {
                // �������� - ���������� ������ ������������ ������ �������
                if (log) fprintf(log, "\n[reduction]");
                for (i = 0; i < ib; i++) {
                    xt = simplex + i*n;
                    for (j = 0; j < n; j++) xt[j] = xb[j] + reduction*(xt[j] - xb[j]);
                    f[i] = val(xt, n, instance);
                }

                for (i = ib+1; i < n+1; i++) {
                    xt = simplex + i*n;
                    for (j = 0; j < n; j++) xt[j] = xb[j] + reduction*(xt[j] - xb[j]);
                    f[i] = val(xt, n, instance);
                }
            }
        }
    }
    printf("Nelder-Mead ended after %d iterataions.\nObjective function value is %g.\nTolerance is %g.\n", current_iteration, fb, current_tol);

    // ���������� � ��������� ������ ���������
    for (j = 0; j < n; j++) {
        x[j] = xb[j];
    }

    // ������������ ������, ���� ��������
    if (!preallocated) {
        free(work_array);
    }

    if (log) {
        fprintf(log, "\n[release memory]\n[best vertix(%g)] ", fb);
        for (j = 0; j < n; j++) {
            fprintf(log, "%g ", x[j]);
        }
        fprintf(log, "\n[End Nelder-Mead]");
    }

    return x;
}
