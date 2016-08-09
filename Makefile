CC = gcc
CFLAGS = -g -Wall
LDFLAGS = -lm
INCLUDES = -Ianns -Ianns/methods -Ianns/nets -Ioptimize -Ioptmize/cg_descent -Iprojects

ANNS_SRCS = anns/anns.c anns/anns_instance.c anns/anns_solvers.c anns/anns_cache.c anns/anns_loadfuncs.c anns/anns_tablecalc.c anns/anns_cond.c anns/anns_methoddata.c anns/anns_tests.c anns/anns_func.c anns/anns_net.c anns/anns_utils.c anns/anns_genpoints.c anns/anns_numeric.c 

METHODS_SRCS = anns/methods/method_hybrid.c anns/methods/method_updatepoints.c anns/methods/method_restarts.c

NETS_SRCS = anns/nets/nrbf_gauss.c anns/nets/rbf_gauss.c

OPTIMIZE_SRCS = optimize/fmin_callback.c optimize/fmin_hj.c optimize/fmin_nm.c optimize/fmin_cg.c optimize/fmin_linesearch.c optimize/cg_descent/cg_descent.c

SRCS = main.c $(ANNS_SRCS) $(NETS_SRCS) $(METHODS_SRCS) $(OPTIMIZE_SRCS)

OBJS = $(SRCS:.c=.o)

MAIN = auroc

.PHONY: depend clean

all: $(MAIN)

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) -c $< -o $@ $(LDFLAGS)

clean:
	$(RM) *.o *~ $(MAIN)

restart:
	@(make clean)
	@(make)
	@(mkdir output)
	@(./auroc)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE
