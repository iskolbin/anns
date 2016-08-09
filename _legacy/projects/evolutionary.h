#ifndef EVOLUTIONARY_H_INCLUDED
#define EVOLUTIONARY_H_INCLUDED

#include "../anns/methods/method_hybrid.h"
#include "../anns/anns_solvers.h"

void evolutionary_1d(double Tm, int fdsteps);
void evolutionary_2d(double Tm, int fdsteps);
void evolutionary_1d_nl(double Tm, int fdsteps);


#endif // EVOLUTIONARY_H_INCLUDED
