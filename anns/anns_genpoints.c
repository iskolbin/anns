#include "anns_cache.h"
#include <math.h>

void ac_genpoints(anns_cache_t *ac, double bounds[]) {
    int j, k, m, cursor = 0, p,  c, T, nx, ny, nz, npts;
    double ax, bx, hx, ay, by, hy, az, bz, hz, r_in, r_out, r, xc;
    double from[64], to[64], pts[64];

    assert(ac);

    for (c = 0; c < ac->nconds; c++) {
        T = (int) bounds[cursor++];

        switch (T) {
            case GP_MANUAL: break;

            case GP_UNIFORM_SOLID_SPHERE: {
                r_in = bounds[cursor++];
                r_out = bounds[cursor++];

                assert(r_in >= 0);
                assert(r_in <= r_out);
                assert(r_out > 0);

//                for (p = 0; p < ac->np[c]; p++) {
                for (p = 0; p < ac->conds[c]->m; p++) {
                    /* Заполняем точку случайными углами */
                    for (j = 0; j < ac->dim-1; j++) {
                        ac->points[c][p][j] = uniform(0, 2*M_PI);
                    }

                    /* Находим точку из интервала, при этом
                            r_in = 0 - шар,
                            r_in == r_out - сфера,
                            r_in < r_out - полый шар */
                    r = uniform(r_in, r_out);

                    for (j = ac->dim-1; j > 1; j--) {
                        ac->points[c][p][j] = r * cos(ac->points[c][p][j-1]);
                        r *= sin(ac->points[c][p][j-1]);
                    }

                    ac->points[c][p][1] = r * sin(ac->points[c][p][0]); /* y */
                    ac->points[c][p][0] = r * cos(ac->points[c][p][0]); /* x */
                }

                /* Перемещаем точки с учётом координат центра */
                for (j = 0; j < ac->dim; j++) {
                    xc = bounds[cursor++];
//                    for (p = 0; p < ac->np[c]; p++) {
                    for (p = 0; p < ac->conds[c]->m; p++) {
                        ac->points[c][p][j] += xc;
                    }
                }

                break;
            }

            case GP_UNIFORM_SPHERE: {
//                r_out = bounds[cursor++];
            }

            case GP_POINTS: {
//                for (p = 0; p < ac->np[c]; p++) {
                for (p = 0; p < ac->conds[c]->m; p++) {
                    for (j = 0; j < ac->dim; j++) {
                        ac->points[c][p][j] = bounds[cursor++];
                    }
                }

                break;
            }

            case GP_UNIFORM_SOLID_CUBOID: {
                for (j = 0; j < ac->dim; j++) {
                    ax = bounds[cursor++]; bx = bounds[cursor++];
//                    for (p = 0; p < ac->np[c]; p++) {
                    for (p = 0; p < ac->conds[c]->m; p++) {
                        ac->points[c][p][j] = uniform(ax, bx);
                    }

                    for (p = 0; p < ac->conds[c]->m_; p++) {
                        ac->points_[c][p][j] = uniform(ax, bx);
                    }
                }
                break;
            }

            case GP_UNIFORM_CUBOID: {
                break;
            }

            case GP_GRID_1D: {
                nx = (int) bounds[cursor++]; ax = bounds[cursor++]; bx = bounds[cursor++];

                hx = nx > 1 ? (bx-ax) / (nx-1) : 0;

                for (p = 0; p < ac->conds[c]->m; p++) {
                    ac->points[c][p][0] = ax + p*hx;
                    printf("%g ",ac->points[c][p][0]);
                }
                break;
            }

            case GP_GRID_2D: {
                nx = (int)bounds[cursor++]; ax = bounds[cursor++]; bx = bounds[cursor++];
                ny = (int)bounds[cursor++]; ay = bounds[cursor++]; by = bounds[cursor++];

                hx = nx > 1 ? (bx-ax) / (nx-1) : 0;
                hy = ny > 1 ? (by-ay) / (ny-1) : 0;

                for (j = 0; j < nx; j++) {
                    for (k = 0; k < ny; k++) {
                        p = j*ny+k;
                        ac->points[c][p][0] = ax + j*hx;
                        ac->points[c][p][1] = ay + k*hy;
                    }
                }

                break;
            }

            case GP_GRID_3D: {
                nx = (int) bounds[cursor++];ax = bounds[cursor++]; bx = bounds[cursor++];
                ny = (int) bounds[cursor++];ay = bounds[cursor++]; by = bounds[cursor++];
                nz = (int) bounds[cursor++];az = bounds[cursor++]; bz = bounds[cursor++];

                hx = nx > 1 ? (bx-ax) / (nx-1) : 0;
                hy = ny > 1 ? (by-ay) / (ny-1) : 0;
                hz = nz > 1 ? (bz-az) / (nz-1) : 0;

                for (j = 0; j < nx; j++) {
                    for (k = 0; k < ny; k++) {
                        for (m = 0; m < nz; m++) {
                            p = j*(ny*nz) + k*nz + m;
                            ac->points[c][p][0] = ax + j*hx;
                            ac->points[c][p][1] = ay + k*hy;
                            ac->points[c][p][2] = az + m*hz;
                        }
                    }
                }
                break;
            }


            case GP_UNIFORM_COIL_SECTOR: {
                double x0, y0, r, R, len, startangle, endangle, angle;
                x0 = bounds[cursor++];
                y0 = bounds[cursor++];
                r = bounds[cursor++];
                R = bounds[cursor++];
                startangle = bounds[cursor++];
                endangle = bounds[cursor++];


                for (p = 0; p < ac->conds[c]->m; p++) {
                    angle = uniform(startangle, endangle);
                    len = uniform(r, R);
                    ac->points[c][p][0] = len * cos(angle) + x0;
                    ac->points[c][p][1] = len * sin(angle) + y0;
                }

                break;
            }

            case GP_UNIFORM_ARC: {
                double x0, y0, R, startangle, endangle, angle;
                x0 = bounds[cursor++];
                y0 = bounds[cursor++];
                R = bounds[cursor++];
                startangle = bounds[cursor++];
                endangle = bounds[cursor++];

                for (p = 0; p < ac->conds[c]->m; p++) {
                    angle = uniform(startangle, endangle);
                    ac->points[c][p][0] = R * cos(angle) + x0;
                    ac->points[c][p][1] = R * sin(angle) + y0;
                }

                break;
            }

            case GP_UNIFORM_CUT_SECTOR: {
                double x0, y0, R, startangle, endangle, angle, r, b, tmp;

                x0 = bounds[cursor++];
                y0 = bounds[cursor++];
                R = bounds[cursor++];
                startangle = bounds[cursor++];
                endangle = bounds[cursor++];

                for (p = 0; p < ac->conds[c]->m; p++) {
                    angle = uniform(startangle, endangle);
//                    b = R * (1 - cos(endangle - startangle)/cos(angle) );
                     b = R * cos(endangle - startangle)/cos(angle) ;
                    r = uniform(b, R);

                    ac->points[c][p][0] = r * cos(angle) + x0;
                    ac->points[c][p][1] = r * sin(angle) + y0;
                    tmp = 0;
                }

                break;
            }

            case GP_UNIFORM_LINE: {
                break;
            }

            case GP_UNIFORM_FRAME: {
                break;
            }
        }
    }
}

void ac_fprintfpoints(FILE *fout, anns_cache_t *ac, int mode) {
    int c, p, j;

    if (fout) {
        switch(mode) {
            case 0: {

                for (c = 0; c < ac->nconds; c++) {
                    fprintf(fout, "Condition %d points (%d):\n[", c, ac->conds[c]->m);
                    for (p = 0; p < ac->conds[c]->m; p++) {
                        fprintf(fout, "[%g", ac->points[c][p][0]);
                        for (j = 1; j < ac->dim; j++) {
                            fprintf(fout, ", %g", ac->points[c][p][j]);
                        }
                        fprintf(fout, "],");
                    }
                    fprintf(fout, "]\n");
                }

                break;
            }
            case 1: {

                for (c = 0; c < ac->nconds; c++) {
                    fprintf(fout, "Condition %d points (%d):\n", c, ac->conds[c]->m);
                    for (j = 0; j < ac->dim; j++) {
                        fprintf(fout, "[%g", ac->points[c][0][j]);
                        for (p = 1; p < ac->conds[c]->m; p++) {
                            fprintf(fout, ", %g", ac->points[c][p][j]);
                        }
                        fprintf(fout, "]\n");
                    }
                    fprintf(fout, "\n");
                }
            }
        }
    }
}
