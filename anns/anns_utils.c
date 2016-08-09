#include "anns_utils.h"

double **d2alloc(int rows, int cols) {
    double **d2 = NULL;
    int i;

    d2 = calloc(rows, sizeof *d2);
    assert(d2);

    for (i = 0; i < rows; i++) {
        d2[i] = calloc(cols, sizeof **d2);
        assert(d2[i]);
    }

    return d2;
}

double **d2allocs(int rows, int *cols) {
    double **d2 = NULL;
    int i;

    d2 = calloc(rows, sizeof **d2);
    assert(d2);

    for (i = 0; i < rows; i++) {
        d2[i] = calloc(cols[i], sizeof **d2);
        assert(d2[i]);
    }

    return d2;
}

double ***d3alloc(int rows, int cols, int depth) {
    double ***d3 = NULL;
    int i, j;

    d3 = calloc(rows, sizeof *d3);
    assert(d3);

    for (i = 0; i < rows; i++) {
        d3[i] = calloc(cols, sizeof **d3);
        assert(d3[i]);
        for (j = 0; j < cols; j++) {
            d3[i][j] = calloc(depth, sizeof ***d3);
            assert(d3[i][j]);
        }
    }

    return d3;
}

void d2free(double **d2, int rows) {
    int i;

    assert(d2);

    for (i = 0; i < rows; i++) {
        assert(d2[i]);
        free(d2[i]);
    }
}

void d3free(double ***d3, int rows, int cols) {
    int i, j;

    assert(d3);

    for (i = 0; i < rows; i++) {
        assert(d3[i]);
        for (j = 0; j < cols; j++) {
            assert(d3[i][j]);
            free(d3[i][j]);
        }
        free(d3[i]);
    }
    free(d3);
}
