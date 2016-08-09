#ifndef STEFAN_H_INCLUDED
#define STEFAN_H_INCLUDED

#include "../anns/methods/method_restarts.h"
#include "../anns/methods/method_hybrid.h"
#include "../anns/methods/method_updatepoints.h"
#include "../anns/anns_solvers.h"

void stefan1phase(void);
void stefan1phasefd(double Tmax, int fdsteps);
void stefan2phase(void);

#endif // STEFAN_H_INCLUDED
