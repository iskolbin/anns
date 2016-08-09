#include "anns_tests.h"

void anns_test_net_funcs(void) {
    int nn1 = 5, nn2 = 10, dim = 3, len = (nn1 + nn2) * (2 + dim), i, j, nearest, farthest;
    double *z = calloc(len, sizeof *z), *neuron, beta, d;
    anns_net_t *net1 = an_new(z, nn1, nn1, 5, 3, NULL);

    srand(123L);
    for (i = 0; i < nn1; i++) {
        neuron = an_neuron(net1, i);
        neuron[0] = uniform(-2,2);
        neuron[1] = 1;
        neuron[2] = uniform(-1, 1);
        neuron[3] = uniform(-1, 1);
        neuron[4] = uniform(-1, 1);
    }
    an_printf(net1);
    printf("\nSwapping 1 and 3\n");
    an_swap(net1, 1, 3);
    an_printf(net1);
    printf("\nDistances list\n");
    for (i = 0; i < nn1; i++) {
        for (j = 0; j < nn1; j++) {
            if (i == j) continue;
            printf("%d to %d distance=%f\n", i, j, an_distance(net1, i, j));
        }
    }
    printf("\nDistances nearest and farthest\n");
    for (i = 0; i < nn1; i++) {
        printf("%d min=%f ", i, an_min_distance(net1, i, &nearest));
        printf("nearest=%d max=%f ", nearest, an_max_distance(net1, i, &farthest));
        printf("farhest=%d\n", farthest);
    }
    printf("\nSet width to min distances\n");
    an_setwidth_nearest(net1, 1.0);
    an_printf(net1);

    free(z);
    an_delete(net1);
}
