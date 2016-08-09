// В данном файле решается уравнение Пуассона на плоскости:
// u_xx + u_yy = f(x,y), f(x,y) = sin(x)sin(y)
// удвлетворяющее граничному условию:
// u_Г = 0,
// где u_xx, u_yy -- вторые производные по x и у соответственно,
// u_Г -- значение функции на границе.
// Граница: квадрат 0 ~ PI
// Аналитическое решение: -0.5sin(x)sin(y)
//
// Решаться задача будет НРБС из 8 нейронов, число точек 30 внутри области и по 10 на границе

// Сначала объявляем прототипы функций, сначала локальных
static void ai_init(anns_instance_t *, fbvp_t *, gbvp_t *, double *);
static void ai_newpoints(anns_instance_t *, double *bounds);
static void ai_initann(anns_instance_t *, double *);
static double ai_mse(anns_instance_t *, double *);

// Потом внешеней, которую будем вызывать в main.c
void poisson_2d_sines(void);

// Теперь объявляем условия задач.
// A: u_xx + u_yy = sin(x) * sin(y)
ANNS_CND(A) {return u_xx[0][0] + u_xx[0][1] - v[0];}
// B: u_Г = 0
ANNS_CND(B) {return u[0];}
// f(x,y) = sin(x)sin(y)
ANNS_FUN(A_f) {return sin(x[0])*sin(x[1]);}

// u_xx[0][0] -- в данной записи означает у первой сети(нумерация с нуля идет) вторая производная по x,
// предполагается, что 0-я переменная это "x", 1-я -- это "y". В данной задаче используется только одна
// сеть с индексом ноль. u_xx[0][1] -- это вторая производная по y. u[0] -- непосредственно значение.
// v[0] -- значение источника. Строго говоря, в v может передаваться что угодно, но в данной задаче это
// именно значение источника.
//
// Далее объявляются формулы-производные сети по параметру ( параметр -- это для НРБС веса, ширины и
// координаты центров ).
ANNS_CND_G(A_a) {return u_axx[0] + u_axx[1];}
ANNS_CND_G(B_a) {return u_a;}

// Инициализация алгоритма
void ai_init(anns_instance_t *ai, fbvp_t *cond, gbvp_t *cond_a, double *delta) {
    // Загружаем настройки, условия, производные условий по параметру, сглаживающие коэффициенты
		ai_conditions(ai, cond, cond_a, delta);
		// Загружаем нейросетевые функции
		ai_loadfuncs(ai);
}

// Генерация точек и заполнение расчет значений источника в них
void ai_newpoints(anns_instance_t *ai, double *bounds) {
    // Генерация точек в указаных рамках
		ai_genpoints(ai, bounds);
		// Заполнение "v" -- в данном случае это значение источника в точках
    ai_eval_v(ai, 0, 0, A_f, NULL);
}

// Инициализация нейронных сетей
void ai_initann(anns_instance_t *ai, double *x) {
    int i, t;

		// Цикл по сетям (на самом деле в данной задаче она одна)
    for (t = 0; t < ai->nnets; t++) {
				// Цикл по нейронам
        for (i = ai->nets[t]->I; i < ai->nets[t]->lenc; ) {
            x[i++] = uniform(-3., 3.);             // Инициализация весов
            x[i++] = M_PI/sqrt(ai->nets[t]->nc);   // Инициализация ширин
            x[i++] = uniform(-M_PI/4,M_PI+M_PI/4); // Инициализация координаты "x" центров
            x[i++] = uniform(-M_PI/4,M_PI+M_PI/4); // Инициализация координаты "y" центров
        }
    }
}

// Расчет среднеквадратичного отклонения для оценки реальной точности полученной модели.
// Алгоритм простой: генерируется 10000 точек в области и вычисляется среднеквадратичная сумма
// разницы значения сети и аналитического решения в точках.
double ai_mse(anns_instance_t *ai, double *x) {
    double s, tmp, point[2];
    int p, M = 10000;

    s = 0;
    for (p = 0; p < M; p++) {
        point[0] = uniform(0, M_PI);
        point[1] = uniform(0, M_PI);

        tmp = 0;
        tmp += anns_eval(point, x, ai);
        tmp += 0.5*sin(point[0])*sin(point[1]);

        s += tmp * tmp;
    }

    return sqrt(s / (M-1) );
}

// Точка входа в задачу
void poisson_2d_sines(void) {
    // Ввод настроек
		int data[] = {
				// Описание задачи -- прямая, эллиптическая, индекс 1, под-индекс 0 --
				// это нужно для формирования имени выходного файла
        PJ_DIRECT, PJ_ELLIPTIC_PROBLEM, 1, 0,
				
				// Функционал -- квадратичный
        FUNCTIONAL_SQUARE, 
				
				// Эвристический расчет сглаживающих коэффициентов -- обратное число всех точек
				DELTA_INV_ALL, 
				
				// Размерность пространства
				2, 
				
				// Число уравнений (5 -- 4 на границе и 1 внутри)
				5,
				
				// Число сетей
				1, 
				
				// Число запоминаемых наборов сетей
				0,

				// Параметры генерации контрольных и проверочных точек; число внешних переменных в уравнении
        
				// Условие в области, генерируем 30 контрольных и 15 проверочных точек; 
				// в уравнении 1 внешняя переменная (источник f(x,y)
				30, 15,  1,         
				// Условия на границе
				10, 5,   0,         
				10, 5,   0, 
        10, 5,   0,
        10, 5,   0, 

				// Параметры нейронной сети
        NORMALIZED | RADIAL | GAUSSIAN, 
				// Число нейронов -- 8
				8, 
				// Размерность пространства для сети
				2, 
				
				0,  
				
				// Расчет производных выхода сети для условий задачи
				// Внутри области требуется 2-я производная
				2, 
				// На границе требуется только рассчет значения сети
				0, 
				0, 
				0, 
				0,};

		// Параметры для генерции точек
    double bounds[] = {
        // в области
				GP_UNIFORM_SOLID_CUBOID, 0, M_PI, 0, M_PI,
        // на границах
				GP_UNIFORM_SOLID_CUBOID, 0, M_PI, M_PI, M_PI,
        GP_UNIFORM_SOLID_CUBOID, 0, 0, 0, M_PI,
        GP_UNIFORM_SOLID_CUBOID, M_PI, M_PI, 0, M_PI,
        GP_UNIFORM_SOLID_CUBOID, 0, M_PI, 0, 0};

    fbvp_t cond[] = {A, B, B, B, B};
    gbvp_t cond_a[] = {A_a, B_a, B_a, B_a, B_a};
    //double delta[] = {1./30, 1./10, 1./10, 1./10, 1./10};

    anns_solverdata_t *as = as_new2(ai_init, ai_newpoints, ai_initann, ai_mse, data, bounds, cond, cond_a, NULL);
    as->Nopt = 5;
    as->gtol = 1e-5;
    as->tol = 1e-5;
    as->n0 = 8;
    anns_solve_restarts(as);
    as_delete(as);
}
