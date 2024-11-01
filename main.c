#include <stdio.h>
#include <math.h>
#include <omp.h>

#define M 160
#define N 180

// это задаем область П - прямоугольник, содержащий нашу трапецию
#define A1 -3.0
#define A2 0.0
#define B1 3.0
#define B2 3.0


// #define EPS 1e-3 // TODO сделать = h^2
#define ACC 8*1e-7 // это точность метода


typedef struct {
    double x, y;
} Point;

typedef struct {
    Point A, B, C, D; // Вершины трапеции
} Trapezoid;

typedef struct {
    Point F, G, H, P; // Вершины прямоугольника TODO переименовать
} Rectangle;


double left_line(Point a){
    return 3.0*a.x + 9 - a.y;
}

double right_line(Point a){
    return -3.0*a.x + 9 - a.y;
}

Point calc_left_point(Point a){
    Point res = {(a.y - 9.0) / 3.0, a.y};
    return res;
}

Point calc_right_point(Point a){
    Point res = {(a.y - 9.0) / -3.0, a.y};
    return res;
}

// Основная функция для нахождения площади пересечения
double intersection_area(Rectangle rect) {

    // мне стыдно за этот код. честно. но он работает :)
    
    Trapezoid trap = {{-3, 0}, {3, 0}, {-2, 3}, {2, 3}};
    
    Point left_point_bottom;
    Point right_point_bottom;

    Point left_point_top;
    Point right_point_top;

    if ((left_line(rect.F) <= 0 && left_line(rect.G) <= 0) || (right_line(rect.F) < 0 && right_line(rect.G) < 0)){
        return 0.0;
    }
    else if (left_line(rect.F) < 0 && left_line(rect.G) > 0 && left_line(rect.H) < 0 && left_line(rect.P) < 0) {
        left_point_top.x = rect.G.x;
        right_point_top.x = rect.G.x;

        left_point_top.y = 3.0 * rect.G.x + 9.0;
        right_point_top.y = 3.0 * rect.G.x + 9.0;

        left_point_bottom = calc_left_point(rect.F);
        right_point_bottom = rect.G;
    }
    else if (right_line(rect.F) > 0 && right_line(rect.G) < 0 && right_line(rect.H) < 0 && right_line(rect.P) < 0) {
        left_point_top.x = rect.F.x;
        right_point_top.x = rect.F.x;

        left_point_top.y = -3.0 * rect.F.x + 9.0;
        right_point_top.y = -3.0 * rect.F.x + 9.0;

        left_point_bottom = rect.F;
        right_point_bottom = calc_right_point(rect.G);
    } else {
        // внутри трапеции
        if(left_line(rect.F) >= 0){
            left_point_bottom = rect.F;
        }
        else {
            left_point_bottom = calc_left_point(rect.F);
        }
        if(right_line(rect.G) >= 0){
            right_point_bottom = rect.G;
        }
        else{
            right_point_bottom = calc_right_point(rect.G);
        }
        if(left_line(rect.H) >= 0){
            left_point_top = rect.H;
        }
        else {
            left_point_top = calc_left_point(rect.H);
        }

        if(right_line(rect.P) >= 0){
            right_point_top = rect.P;
        }
        else {
            right_point_top = calc_right_point(rect.P);
        }
    }

    return (right_point_bottom.x - left_point_bottom.x + right_point_top.x - left_point_top.x) / 2.0 * (left_point_top.y - left_point_bottom.y);
}

// Находит определенный интеграл, когда линия вертикальная (внутри трапеции)
double lengh_of_vert_line(Point top, Point bottom){

    Point top_point, bottom_point;

    if(left_line(bottom) <= 0 || right_line(bottom) <= 0) {
        return 0.0;
    } else if (left_line(bottom) > 0 && left_line(top) < 0){
        top_point.y = 3.0 * bottom.x + 9.0;  // TODO вынести в функцию
        top_point.x = top.x;
        bottom_point = bottom;
    } else if (right_line(bottom) > 0 && left_line(top) < 0){
        top_point.y = -3.0 * bottom.x + 9.0;
        top_point.x = top.x;
        bottom_point = bottom;
    }
    else {
        top_point = top;
        bottom_point = bottom;
    }

    return (top_point.y - bottom_point.y);
}

// Находит определенный интеграл, когда линия горизонтальная (внутри трапеции)
double lengh_of_horz_line(Point left, Point right){

    Point left_point, right_point;

    if(left_line(right) <= 0 || right_line(left) <= 0) {
        return 0.0;
    } else if (left_line(right) > 0 && left_line(left) < 0){
        left_point = calc_left_point(right);
        right_point = right;
    } else if (right_line(left) > 0 && left_line(right) < 0){
        left_point = left;
        right_point = calc_right_point(left);
    }
    else {
        left_point = left;
        right_point = right;
    }

    return (right_point.x - left_point.x);
}

// Возвращает значение правой части (10)
double f(double x_l, double y_l, double x_r, double y_r, double h1, double h2) {
    Rectangle rect;

    rect.F = (Point){x_l, y_l};
    rect.G = (Point){x_r, y_l};

    rect.H = (Point){x_l, y_r};
    rect.P = (Point){x_r, y_r};

    return intersection_area(rect) / (h1 * h2);

}

int main() {
    double h1 = (B1 - A1) / M;
    double h2 = (B2 - A2) / N;

    double EPS = h1 > h2 ? h1*h1 : h2*h2;

    double w[M][N];
    double r[M][N];
    double a[M][N], b[M][N];
    double F[M][N];
    double Diff[M][N];
    double Ar[M][N];
    int i, j, k;
    int err_cnt = 0;

    double err_by_iter[20];
    int step = 12000;
    int step_2 = 24000;

    for (i = 0; i < 20; i++){
        err_by_iter[i] = 0.0;
    }

    // Инициализация
    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            w[i][j] = 0.0;
            r[i][j] = 0.0;
            a[i][j] = 0.0;
            b[i][j] = 0.0;
            F[i][j] = 0.0;
            Diff[i][j] = 0.0;  // TODO убрать
            Ar[i][j] = 0.0;
        }
    }

    double start = omp_get_wtime();

    // Вычисление коэффициентов a, b и правой части F
    for (i = 1; i < M-1; i++) {
        for (j = 1; j < N-1; j++) {
            double xi = A1 + i*h1;
            double yj = A2 + j*h2;

            double lij = lengh_of_vert_line((Point){xi-0.5*h1, yj+0.5*h2}, (Point){xi-0.5*h1, yj-0.5*h2});
            double gij = lengh_of_horz_line((Point){xi-0.5*h1, yj-0.5*h2}, (Point){xi+0.5*h1, yj-0.5*h2});

            a[i][j] = lij / h2 + (1.0 - lij / h2) / EPS;
            b[i][j] = gij / h1 + (1.0 - gij / h1) / EPS;

            F[i][j] = f(xi-0.5*h1, yj-0.5*h2, xi+0.5*h1, yj+0.5*h2, h1, h2);
        }
    }

    // Метод скорейшего спуска
    double tau, norm_r, norm_dr, norm;
    int iter = 0;
    do {
        norm_r = 0.0;
        norm_dr = 0.0;
        norm = 0.0;
        for (i = 1; i < M-1; i++) {
            for (j = 1; j < N-1; j++) {
                // тут вычисляем невязку
                r[i][j] = -(a[i+1][j] * (w[i+1][j] - w[i][j]) / h1 - a[i][j] * (w[i][j] - w[i-1][j]) / h1) / h1
                        - (b[i][j+1] * (w[i][j+1] - w[i][j]) / h2 - b[i][j] * (w[i][j] - w[i][j-1]) / h2) / h2
                        - F[i][j];
            }
        }

        for (i = 1; i < M-1; i++) {
            for (j = 1; j < N-1; j++) {
                // тут - матрицу A*r
                Ar[i][j] = -(a[i+1][j] * (r[i+1][j] - r[i][j]) / h1 - a[i][j] * (r[i][j] - r[i-1][j]) / h1) / h1
                        - (b[i][j+1] * (r[i][j+1] - r[i][j]) / h2 - b[i][j] * (r[i][j] - r[i][j-1]) / h2) / h2;
                
                norm_r += r[i][j] * r[i][j] * h1 * h2;
                norm_dr += Ar[i][j] * r[i][j] * h1 * h2;
            }
        }

        tau = norm_r / norm_dr;

        // тут немного неэффективно считаем норму для проверки условия останова + уравнение (15)
        for (i = 1; i < M; i++) {
            for (j = 1; j < N; j++) {
                double tmp = w[i][j];
                w[i][j] -= tau * r[i][j];
                Diff[i][j] = w[i][j] - tmp;
                norm += Diff[i][j] * Diff[i][j] * h1 * h2;
            }
        }

        norm = sqrt(norm);

        if(!(iter % step) && err_cnt < 19){
            err_by_iter[err_cnt++] = norm;
        }

        // if(!(iter % step_2)) {
        //     for (i = 0; i < M; i++) {
        //         for (j = 0; j < N; j++) {
        //             printf("%f ", r[i][j]);
        //         }
        //         printf("\n");
        //     }
        //     printf("==========================================================================================================\n");
        // }

        printf("Тау: %lf, norm_r: %lf, norm_dr%lf, norm: %lf\n", tau, norm_r, norm_dr, norm);
        iter++;
    } while (norm > ACC);

    printf("Число итераций: %d\n", iter);

    double end = omp_get_wtime();
    printf("Затраченное время: %f\n", end - start);


    // for (i = 0; i < 20; i++){
    //     printf("%f, ", err_by_iter[i]);
    // }


    // // Вывод решения
    // for (i = 0; i < M; i++) {
    //     for (j = 0; j < N; j++) {
    //         printf("%f ", w[i][j]);
    //     }
    //     printf("\n");
    // }

    return 0;
}
