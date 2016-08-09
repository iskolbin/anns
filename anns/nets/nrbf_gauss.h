#ifndef ANNS_NETS_NRBF_GAUSS_H_INCLUDED
#define ANNS_NETS_NRBF_GAUSS_H_INCLUDED

// ����-� (NRBF-G)
// NORMALIZED | RADIAL | GAUSSIAN

#include <math.h>
#include <stdio.h>
#include "../anns_types.h"
//#include "../anns_tablecalc.h"

// ���������� ���� �� ���������� � ����������� (����������, �� ��� �� �����)
//
// � ������� �� �������� ���������� ���������, ��� �������������� ����������� ����� �������
// ����������, � ���������� �� ����. � �������������� ����� ������ ��� ������ �������
// ��������� � ����������� �� ���������, �.�. ���������� ������������� ��������� ����� �����������
// ����� ���������������. ���������� ����������� �������� �������� �������� ����������, ��� ���� �����������
// �� ����� �������� (���� ��������� �� ��������). ���������:
// 1) ���������� ����� ������ �������, � ���� �� �� ����� ������������
// ��������� �� ������� ���������� (��-�� ����������� ����������� �������������� ����)
// ��� ����������������� ����� ��� �� ������������� (����� ���� ���������), � � ������ ���� ���������
// ������� 0/0.
// 2) �������� ��������, ����� ����� �������� �������: ����� ����������, � ����� ������
// ��� ����� ����� �����, ����� ���������� ��������� ����� ������ ���������� ��������
// ����������.
// 3) ����� �������, ��������� ���� ��� ���������� ������ �� ����� ������, �.�. �������������� ��������������
// �������� ������ ������������� (nmax �������� ��������� � ���������), � ������ �����������

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
