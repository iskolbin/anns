#ifndef ANNS_CONST_H_INCLUDED
#define ANNS_CONST_H_INCLUDED

// Описание задачи (влияет только на название файла
#define PJ_DIRECT                      100
#define PJ_INVERSE_UNKNOWN_SOURCE      101
#define PJ_INVERSE_UNKNOWN_BOUNDARY    102
#define PJ_INVERSE_UNKNOWN_INITIAL     103
#define PJ_INVERSE_UNKNOWN_COEFFITIENT 104

#define PJ_ELLIPTIC_PROBLEM   200
#define PJ_PARABOLIC_PROBLEM  201
#define PJ_HYPERBOLIC_PROBLEM 202
#define PJ_APPROXIMATION      203

// Метод генерации выборочных точек
#define GP_MANUAL 1100 /* Ручной */
#define GP_POINTS 1101 /* Потеченый */

#define GP_UNIFORM_CUBOID     1210 /* Равномерное распределение на поверхности кубоида */
#define GP_UNIFORM_SPHERE     1211 /* Равномерное распределение на сфере */
#define GP_UNIFORM_TORUS      1212 /* Равномерное распределение на поверхности тора */
#define GP_UNIFORM_CYLINDER   1213 /* Равномерное распределение на поверхности цилиндра */

#define GP_UNIFORM_SOLID_CUBOID   1300 /* Равномерное распределение в кубоиде */
#define GP_UNIFORM_SOLID_SPHERE   1301 /* Равномерное распределение в шаре */
#define GP_UNIFORM_SOLID_TORUS    1302 /* Равномерное распределение в торе */
#define GP_UNIFORM_SOLID_CYLINDER 1303 /* Равномерное распределение в цилиндре */

#define GP_UNIFORM_LINE  1250 /* Отрезок: координаты концов */
#define GP_UNIFORM_FRAME 1251 /* Рамка: число точек, координаты концов */
#define GP_UNIFORM_ARC   1252 /* Арка: координаты центра, радиус, начальный угол относительно Oy, конечный угол относительно Oy */
#define GP_UNIFORM_COIL_SECTOR 1253 /* Сектор кольца: координаты центра, внутренний радиус, внещний радиус, начальный угол относительно Oy, конечный угол относительно Oy */
#define GP_UNIFORM_CUT_SECTOR  1254

#define GP_GRID_MD 1400 /* Многомерная регулярная прямоугольная сетка */
#define GP_GRID_1D 1401 /* Одномерная регулярная сетка */
#define GP_GRID_2D 1402 /* Двумерная регулярная сетка */
#define GP_GRID_3D 1403 /* Трехмерная регулярная сетка */

// Тип нейронной сети
// Радиально-базисные сети
#define RADIAL   256
#define ELLIPTIC 512

#define CLASSIC    1024
#define NORMALIZED 2048

#define GAUSSIAN              1
#define MULTIQUADRIC          2
#define INVERSE_MULTIQUADRIC  4
#define THIN_PLATE_SPLINE     8
#define CAUCHY_FUNCTION       16

// Прочие глобальные константы
#define MIXED_DERIVATIVES 0x7000  /* Маска для степени производной, указывающая на смешаныне производные */
#define MAX_DERIVATIVES 2         /* Максимальная степень производной */
#define OUTPUT_DIR "output"       /* Папка для выходных данных */
#define EVAL_MANUAL -1            /* Псевдо-степень производной - вычисления не проводятся */

// Типы функционалов
#define FUNCTIONAL_SQUARE       0 /* Квадратичный f^2 */
#define FUNCTIONAL_ABSOLUTE     1 /* Абсолютный |f|, не особо работает*/
#define FUNCTIONAL_LINEAR       2 /* Линейный f */
#define FUNCTIONAL_SQUAREROOT   3 /* Квадратичный с последующим извлечением корня */

// Значения поправочных коэффициентов
#define DELTA_MANUAL     0 /* Ручная установка (или все единичные) */
#define DELTA_INV_COND   1 /* Обратное число точек условия (1/m) */
#define DELTA_INV_ALL    2 /* Обратное число всех точек (1/M) */

// Ключи для параметров методов
#define AM_NULL       0 // Нужно ставить в конце
#define AM_SEED       1
#define AM_TOL        2
#define AM_GTOL       3
#define AM_NOFILEOUT  4

#define AM_DATA       50
#define AM_COND       51
#define AM_COND_A     52
#define AM_DELTA      53
#define AM_BOUNDS     54
#define AM_COND_I     55
#define AM_COND_I_A   56

#define AM_UPDATE     101
#define AM_MSE        102
#define AM_INIT       103
#define AM_NEWPOINTS  104
#define AM_INITANN    105
#define AM_PREEVAL    106
#define AM_UPDATEPOINTS 107

#define AM_OPTSTEPS   150
#define AM_GROW       151
#define AM_REDUCE     152
#define AM_REGENERATE 153
#define AM_TRIAL      154
#define AM_MAXTICKS   155
#define AM_INITIAL_N  156
#define AM_MAXIMAL_N  157
#define AM_MAXTRIALS  158
#define AM_MAXITERS   159

#define AM_TBEGIN     201
#define AM_TEND       202
#define AM_TSTEPS     203

#endif
