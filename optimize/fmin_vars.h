#ifndef FMIN_VARS_H_INCLUDED
#define FMIN_VARS_H_INCLUDED

// ����� ������ ������� ���������� ����������
// val, grad, valgrad, x, xtemp, g, gtemp, d,
// ��������� f0, fa = 0, ��������� g0, ga = 0,
// ���������� n, instance, nval = 0, ngrad = 0

// ����� ������ � xtemp ����� ����� �����������, � gtemp - ����� ���,
// � fa - �������� � ����� ����������� �������,
// � nval ������������� ����� ���������� �������,
// � ngrad - ��������� � ������� ��������� ������
typedef struct {
    double (*val)(double *, int, void *);               // ������� �������
    void (*grad)(double *, double *, int, void *);      // ������� ��� ���������� ���������
    double (*valgrad)(double *, double *, int, void *); // ������� ��� ���������� ��������� � ��������

    double *x;      // ������� ����������� (�� �������� � �������� ������)
    double *xtemp;  // ������ ��� ������ ����������� (�����������)
    double *g;      // �������� � ������� ����������� (�� �����������)
    double *gtemp;  // ������ ��� ��������� � ����� ����������� (�����������)
    double *d;      // ������ ����������� (�� �����������)

    double old_fa;

    double a;   // ���
    double f0;  // �������� ������� � ������� ����������� f(x)
    double fa;  // �������� ������� � ����� ����������� f(x + ad)
    double g0;  // �������� (d, g(x))
    double ga;  // ��������� (d, g(x + ad))

    int n;          // ����������� ������
    void *instance; // ������������ ����������

    int niterations; // ����� ��������
    int nval;   // ����� ���������� �������
    int ngrad;  // ����� ���������� ���������
} fmin_vars_t;

#endif // FMIN_VARS_H_INCLUDED
