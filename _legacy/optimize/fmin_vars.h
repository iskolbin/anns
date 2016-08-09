#ifndef FMIN_VARS_H_INCLUDED
#define FMIN_VARS_H_INCLUDED

// Перед первым вызовом необходимо установить
// val, grad, valgrad, x, xtemp, g, gtemp, d,
// Вычислить f0, fa = 0, вычислить g0, ga = 0,
// установить n, instance, nval = 0, ngrad = 0

// После поиска в xtemp будет новое приближенеи, в gtemp - новый шаг,
// в fa - значение в новом приближении функции,
// в nval накапливается число вычислений функции,
// в ngrad - градиента в течении линейного поиска
typedef struct {
    double (*val)(double *, int, void *);               // Целевая функция
    void (*grad)(double *, double *, int, void *);      // Функция для вычисления градиента
    double (*valgrad)(double *, double *, int, void *); // Функция для вычисления градиента и значения

    double *x;      // Текущее приближение (не меняется в процессе поиска)
    double *xtemp;  // Вектор для нового приближения (обновляется)
    double *g;      // Градиент в текущем приближении (не обновляется)
    double *gtemp;  // Вектор для градиента в новом приближении (обновляется)
    double *d;      // Вектор направлений (не обновляется)

    double old_fa;

    double a;   // Шаг
    double f0;  // Значение функции в текущем приближении f(x)
    double fa;  // Значение функции в новом приближении f(x + ad)
    double g0;  // Значение (d, g(x))
    double ga;  // Значнение (d, g(x + ad))

    int n;          // Размерность задачи
    void *instance; // Произвольная переменная

    int niterations; // Число итераций
    int nval;   // Число вычислений функции
    int ngrad;  // Число вычислений градиента
} fmin_vars_t;

#endif // FMIN_VARS_H_INCLUDED
