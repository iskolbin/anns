#include "lbfgs_geterror.h"

char lbfgs_errors_strings[][40] = {
    "SUCCESS",
    "LBFGSERR_UNKNOWNERROR",
    "LBFGSERR_LOGICERROR",
    "LBFGSERR_OUTOFMEMORY",
    "LBFGSERR_CANCELED",
    "LBFGSERR_INVALID_N",
    "LBFGSERR_INVALID_N_SSE",
    "LBFGSERR_INVALID_X_SSE",
    "LBFGSERR_INVALID_EPSILON",
    "LBFGSERR_INVALID_TESTPERIOD",
    "LBFGSERR_INVALID_DELTA",
    "LBFGSERR_INVALID_LINESEARCH",
    "LBFGSERR_INVALID_MINSTEP",
    "LBFGSERR_INVALID_MAXSTEP",
    "LBFGSERR_INVALID_FTOL",
    "LBFGSERR_INVALID_WOLFE",
    "LBFGSERR_INVALID_GTOL",
    "LBFGSERR_INVALID_XTOL",
    "LBFGSERR_INVALID_MAXLINESEARCH",
    "LBFGSERR_INVALID_ORTHANTWISE",
    "LBFGSERR_INVALID_ORTHANTWISE_START",
    "LBFGSERR_INVALID_ORTHANTWISE_END",
    "LBFGSERR_OUTOFINTERVAL",
    "LBFGSERR_INCORRECT_TMINMAX",
    "LBFGSERR_ROUNDING_ERROR",
    "LBFGSERR_MINIMUMSTEP",
    "LBFGSERR_MAXIMUMSTEP",
    "LBFGSERR_MAXIMUMLINESEARCH",
    "LBFGSERR_MAXIMUMITERATION",
    "LBFGSERR_WIDTHTOOSMALL",
    "LBFGSERR_INVALIDPARAMETERS",
    "LBFGSERR_INCREASEGRADIENT"};

char* lbfgs_geterror(int code){
    if (code < 0) {
        return lbfgs_errors_strings[1024+code+1];
    } else {
        return lbfgs_errors_strings[0];
    }
}