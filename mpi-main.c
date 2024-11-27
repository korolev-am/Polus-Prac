#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"

// это задаем область П - прямоугольник, содержащий нашу трапецию
#define A1 -3.0
#define A2 0.0
#define B1 3.0
#define B2 3.0

#define M 41
#define N 41

// #define ACC 1e-6 // это точность метода
#define ACC 1e-6 // это точность метода


typedef struct {
    double x, y;
} Point;

typedef struct {
    Point A, B, C, D; // Вершины трапеции
} Trapezoid;

typedef struct {
    Point F, G, H, P; // Вершины прямоугольника
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
        top_point.y = 3.0 * bottom.x + 9.0;
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

double w[M][N];
double r[M][N];
// double a[M][N], b[M][N];
// double F[M][N];
double Ar[M][N];


int obl_index_to_rank(int i, int j, int n, int k){
    return i * k + j % n;
}


void initialize(double **a, double **b, double **F, double **w, double **r, double **Ar){
    // Инициализация

    int i = 0;
    int j = 0;

    // #pragma omp parallel for collapse(2) private(i, j)
    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            w[i][j] = 0.0;
            r[i][j] = 0.0;
            a[i][j] = 0.0;
            b[i][j] = 0.0;
            F[i][j] = 0.0;
            Ar[i][j] = 0.0;
        }
    }
}


void calc_coefs(double h1, double h2, double EPS, double **F, double **a, double **b, int row_start, int row_end, int col_start, int col_end){
    double xi, yj, lij, gij, tmp;
    int i, j;
    // Вычисление коэффициентов a, b и правой части F
    // #pragma omp parallel for collapse(2) private(i, j, xi, yj, lij, gij)
    for (i = row_start; i <= row_end; i++) {
        for (j = col_start; j <= col_end; j++) {
            if(i == 0 || j == 0){
                continue;
            }

            xi = A1 + i*h1;
            yj = A2 + j*h2;

            lij = lengh_of_vert_line((Point){xi-0.5*h1, yj+0.5*h2}, (Point){xi-0.5*h1, yj-0.5*h2});
            gij = lengh_of_horz_line((Point){xi-0.5*h1, yj-0.5*h2}, (Point){xi+0.5*h1, yj-0.5*h2});

            a[i][j] = lij / h2 + (1.0 - lij / h2) / EPS;
            b[i][j] = gij / h1 + (1.0 - gij / h1) / EPS;

            if(i != M-1 && j != N-1) {
                F[i][j] = f(xi-0.5*h1, yj-0.5*h2, xi+0.5*h1, yj+0.5*h2, h1, h2);
            }
        }
    }
}


int MRD(double h1, double h2, double **a, double **b, double **F, double **w, double **r, double **Ar, 
        int row_start, int row_end, int col_start, int col_end, int obj_i, int obj_j, int n, int k){
    // n - по вертикали
    // k - по горизонтали
    // Метод скорейшего спуска
    double tau, norm_r, norm_dr, norm, tmp;
    double global_norm, global_norm_dr, global_norm_r;
    int i, j;
    int iter = 0, msg_iter = 0;

    double *right_send_buf, *right_recv_buf, *left_send_buf, *left_recv_buf, *top_send_buf, *top_recv_buf, *down_send_buf, *down_recv_buf;

    if(obj_j + 1 < k){
        // если есть сосед справа
        // будем накапливать для него элементы
        right_send_buf = (double *) malloc ((row_end - row_start - 1) * sizeof(double));
        right_recv_buf = (double *) malloc ((row_end - row_start - 1) * sizeof(double));
    }
    if(obj_j - 1 >= 0){
        // если есть сосед слева
        left_send_buf = (double *) malloc ((row_end - row_start - 1) * sizeof(double));
        left_recv_buf = (double *) malloc ((row_end - row_start - 1) * sizeof(double));
    }
    // if(obj_i + 1 < n){
    if(obj_i - 1 >= 0){
        // если есть сосед сверху
        // а тут достаточно указателя на строку
        // top_send_buf;
        top_recv_buf = (double *) malloc ((col_end - col_start - 1) * sizeof(double));
    }
    // if(obj_i - 1 >= 0){
    if(obj_i + 1 < n){
        // если есть сосед снизу
        // down_send_buf;
        down_recv_buf = (double *) malloc ((col_end - col_start - 1) * sizeof(double));
    }

    do {
        norm_r = 0.0;
        norm_dr = 0.0;
        norm = 0.0;
        global_norm = 0.0;
        global_norm_dr = 0.0;
        global_norm_r = 0.0;

        MPI_Request recv_right_r, recv_left_r, recv_top_r, recv_down_r, recv_right_Ar, recv_left_Ar, recv_top_Ar, recv_down_Ar;

        int right_send_iter = 0, left_send_iter = 0, right_recv_iter = 0, left_recv_iter = 0, top_recv_iter = 0, down_recv_iter = 0;
        if(obj_j + 1 < k){
            // если есть сосед справа
            // будем накапливать для него элементы
            MPI_Irecv(right_recv_buf, row_end - row_start - 1, MPI_DOUBLE, obl_index_to_rank(obj_i, obj_j + 1, n, k), msg_iter, MPI_COMM_WORLD, &recv_right_r);
        }
        if(obj_j - 1 >= 0){
            // если есть сосед слева
            MPI_Irecv(left_recv_buf, row_end - row_start - 1, MPI_DOUBLE, obl_index_to_rank(obj_i, obj_j - 1, n, k), msg_iter, MPI_COMM_WORLD, &recv_left_r);
        }
        if(obj_i - 1 >= 0){
            // если есть сосед сверху
            MPI_Irecv(top_recv_buf, col_end - col_start - 1, MPI_DOUBLE, obl_index_to_rank(obj_i - 1, obj_j, n, k), msg_iter, MPI_COMM_WORLD, &recv_top_r);
        }
        if(obj_i + 1 < n){
            // если есть сосед снизу
            MPI_Irecv(down_recv_buf, col_end - col_start - 1, MPI_DOUBLE, obl_index_to_rank(obj_i + 1, obj_j, n, k), msg_iter, MPI_COMM_WORLD, &recv_down_r);
        }


        for (i = row_start + 1; i < row_end; i++) {
            for (j = col_start + 1; j < col_end; j++) {

                if(i == 0 || j == 0 || i == M-1 || j == N-1){
                    continue;
                }

                // тут вычисляем невязку
                r[i][j] = -(a[i+1][j] * (w[i+1][j] - w[i][j]) / h1 - a[i][j] * (w[i][j] - w[i-1][j]) / h1) / h1
                        - (b[i][j+1] * (w[i][j+1] - w[i][j]) / h2 - b[i][j] * (w[i][j] - w[i][j-1]) / h2) / h2
                        - F[i][j];
                // printf("%f ", r[i][j]);
                norm_r += r[i][j] * r[i][j] * h1 * h2;


                if(i == row_start + 1 && j == col_end - 1 && obj_i - 1 >= 0){
                    // отправляем соседу сверху
                    MPI_Request send_up;
                    MPI_Isend(r[i], col_end - col_start - 1, MPI_DOUBLE, obl_index_to_rank(obj_i - 1, obj_j, n, k), msg_iter, MPI_COMM_WORLD, &send_up);
                }
                if(i == row_end - 1 && j == col_end - 1 && obj_i + 1 < n){
                    // отправляем соседу снизу
                    MPI_Request send_down;
                    MPI_Isend(r[i], col_end - col_start - 1, MPI_DOUBLE, obl_index_to_rank(obj_i + 1, obj_j, n, k), msg_iter, MPI_COMM_WORLD, &send_down);
                }
                if(j == col_start + 1 && obj_j - 1 >= 0){
                    left_send_buf[left_send_iter++] = r[i][j];
                }
                if(j == col_end - 1 && obj_j + 1 < k){
                    right_send_buf[right_send_iter++] = r[i][j];
                }
                // printf("\n");
            }
        }

        // printf("After r calc (%d, %d)\n", obj_i, obj_j);

        if(obj_j + 1 < k){
            // отправляем соседу справа
            MPI_Request send_right;
            MPI_Isend(right_send_buf, row_end - row_start - 1, MPI_DOUBLE, obl_index_to_rank(obj_i, obj_j + 1, n, k), msg_iter, MPI_COMM_WORLD, &send_right);

            MPI_Wait(&recv_right_r, MPI_STATUS_IGNORE);

            for(i = row_start + 1; i < row_end; i++){
                r[i][col_end] = right_recv_buf[right_recv_iter++];
            }

            MPI_Irecv(right_recv_buf, row_end - row_start - 1, MPI_DOUBLE, obl_index_to_rank(obj_i, obj_j + 1, n, k), msg_iter+1, MPI_COMM_WORLD, &recv_right_Ar);
        }
        if(obj_j - 1 >= 0){
            // отправляем соседу слева
            MPI_Request send_left;
            MPI_Isend(left_send_buf, row_end - row_start - 1, MPI_DOUBLE, obl_index_to_rank(obj_i, obj_j - 1, n, k), msg_iter, MPI_COMM_WORLD, &send_left);

            MPI_Wait(&recv_left_r, MPI_STATUS_IGNORE);

            for(i = row_start + 1; i < row_end; i++){
                r[i][col_start] = left_recv_buf[left_recv_iter++];
            }

            MPI_Irecv(left_recv_buf, row_end - row_start - 1, MPI_DOUBLE, obl_index_to_rank(obj_i, obj_j - 1, n, k), msg_iter+1, MPI_COMM_WORLD, &recv_left_Ar);
        }
        if(obj_i - 1 >= 0){
            // получаем от соседа сверху
            MPI_Wait(&recv_top_r, MPI_STATUS_IGNORE);

            for(j = col_start + 1; j < col_end; j++){
                r[row_start][j] = top_recv_buf[top_recv_iter++];
            }

            MPI_Irecv(top_recv_buf, col_end - col_start - 1, MPI_DOUBLE, obl_index_to_rank(obj_i - 1, obj_j, n, k), msg_iter+1, MPI_COMM_WORLD, &recv_top_Ar);
        }
        if(obj_i + 1 < n){
            // получаем от соседа снизу
            MPI_Wait(&recv_down_r, MPI_STATUS_IGNORE);

            for(j = col_start + 1; j < col_end; j++){
                r[row_end][j] = down_recv_buf[down_recv_iter++];
            }

            MPI_Irecv(down_recv_buf, col_end - col_start - 1, MPI_DOUBLE, obl_index_to_rank(obj_i + 1, obj_j, n, k), msg_iter+1, MPI_COMM_WORLD, &recv_down_Ar);
        }

        left_send_iter = 0;
        right_send_iter = 0;
        right_recv_iter = 0;
        left_recv_iter = 0;
        top_recv_iter = 0;
        down_recv_iter = 0;

        msg_iter++;

        // printf("After 1 wait");


        for (i = row_start + 1; i < row_end; i++) {
            for (j = col_start + 1; j < col_end; j++) {

                if(i == 0 || j == 0 || i == M-1 || j == N-1){
                    continue;
                }
                // тут - матрицу A*r
                Ar[i][j] = -(a[i+1][j] * (r[i+1][j] - r[i][j]) / h1 - a[i][j] * (r[i][j] - r[i-1][j]) / h1) / h1
                        - (b[i][j+1] * (r[i][j+1] - r[i][j]) / h2 - b[i][j] * (r[i][j] - r[i][j-1]) / h2) / h2;
                norm_dr += Ar[i][j] * r[i][j] * h1 * h2;


                if(i == row_start + 1 && j == col_end - 1 && obj_i - 1 >= 0){
                    // отправляем соседу сверху
                    MPI_Request send_up;
                    MPI_Isend(Ar[i], col_end - col_start - 1, MPI_DOUBLE, obl_index_to_rank(obj_i - 1, obj_j, n, k), msg_iter, MPI_COMM_WORLD, &send_up);
                }
                if(i == row_end - 1 && j == col_end - 1 && obj_i + 1 < n){
                    // отправляем соседу снизу
                    MPI_Request send_down;
                    MPI_Isend(Ar[i], col_end - col_start - 1, MPI_DOUBLE, obl_index_to_rank(obj_i + 1, obj_j, n, k), msg_iter, MPI_COMM_WORLD, &send_down);
                }
                if(j == col_start + 1 && obj_j - 1 >= 0){
                    left_send_buf[left_send_iter++] = Ar[i][j];
                }
                if(j == col_end - 1 && obj_j + 1 < k){
                    right_send_buf[right_send_iter++] = Ar[i][j];
                }
            }
        }

        // printf("After Ar calc");

        if(obj_j + 1 < k){
            // отправляем соседу справа
            // printf("After Ar sending to right (%d, %d)\n", obj_i, obj_j);
            MPI_Request send_right;
            MPI_Isend(right_send_buf, row_end - row_start - 1, MPI_DOUBLE, obl_index_to_rank(obj_i, obj_j + 1, n, k), msg_iter, MPI_COMM_WORLD, &send_right);

            MPI_Wait(&recv_right_Ar, MPI_STATUS_IGNORE);

            for(i = row_start; i < row_end; i++){
                Ar[i][col_end] = right_recv_buf[right_recv_iter++];
            }
        }
        if(obj_j - 1 >= 0){
            // отправляем соседу слева
            MPI_Request send_left;
            MPI_Isend(left_send_buf, row_end - row_start - 1, MPI_DOUBLE, obl_index_to_rank(obj_i, obj_j - 1, n, k), msg_iter, MPI_COMM_WORLD, &send_left);

            MPI_Wait(&recv_left_Ar, MPI_STATUS_IGNORE);

            for(i = row_start; i < row_end; i++){
                Ar[i][col_start] = left_recv_buf[left_recv_iter++];
            }
            // printf("After Ar sending to left (%d, %d)\n", obj_i, obj_j);
        }
        if(obj_i - 1 >= 0){
            // printf("After Ar sending to top (%d, %d)\n", obj_i, obj_j);
            // получаем от соседа сверху
            MPI_Wait(&recv_top_Ar, MPI_STATUS_IGNORE);

            for(j = col_start; j < col_end; j++){
                Ar[row_start][j] = top_recv_buf[top_recv_iter++];
            }

        }
        if(obj_i + 1 < n){
            // получаем от соседа снизу
            // printf("After Ar sending to down (%d, %d)\n", obj_i, obj_j);
            MPI_Wait(&recv_down_Ar, MPI_STATUS_IGNORE);

            for(j = col_start; j < col_end; j++){
                Ar[row_end][j] = down_recv_buf[down_recv_iter++];
            }
        }

        // printf("After 2 wait");


        MPI_Allreduce(&norm_r, &global_norm_r, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&norm_dr, &global_norm_dr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        tau = global_norm_r / global_norm_dr;

        for (i = row_start; i < row_end; i++) {
            for (j = col_start; j < col_end; j++) {

                if(i == 0 || j == 0){
                    continue;
                }
                tmp = w[i][j];
                w[i][j] -= tau * r[i][j];
                tmp = w[i][j] - tmp;
                norm += tmp * tmp * h1 * h2;
            }
        }

        // printf("Before Allreduce Tau: %f, norm_r: %f, norm_dr: %f, global_norm: %f\n", tau, norm_r, norm_dr, global_norm);


        // printf("iter %d\n", iter);

        MPI_Allreduce(&norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        global_norm = sqrt(global_norm);

        printf("global_norm: %f, norm_r: %f, norm_dr: %f, tau: %f, iter: %d\n", global_norm, norm_r, norm_dr, tau, iter);
        iter++;
        msg_iter++;

    } while (global_norm > ACC);
    // } while (iter < 10);

    printf("Free memory in MRD");

    if(obj_j + 1 < k){
        free(right_send_buf);
        free(right_recv_buf);
    }
    if(obj_j - 1 >= 0){
        free(left_send_buf);
        free(left_recv_buf);
    }
    if(obj_i - 1 >= 0){
        // free(top_send_buf);
        free(top_recv_buf);
    }
    if(obj_i + 1 < n){
        // free(down_send_buf);
        free(down_recv_buf);
    }

    return iter;
}


void make_experiment(int num_treads) {

    // omp_set_num_threads(num_treads);
    // printf("Число потоков: %d", omp_get_num_threads());

    double start;
    double end;
    double h1 = (B1 - A1) / M;
    double h2 = (B2 - A2) / N;
    double EPS = h1 > h2 ? h1*h1 : h2*h2;

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Определяем количество подобластей по горизонтали и вертикали
    int k = size / 2 + size % 2;
    int n = size / k;

    int obj_i = rank / k;
    int obj_j = rank % k;

    // Определяем границы подматрицы для текущего rank
    int row_start = N/n * obj_i;
    int row_end = N/n + 1 + obj_i*(N/n + 1) - n*obj_i;
    int col_start = N/k * obj_j;
    int col_end = N/k + 1 + obj_j*(N/k + 1) - k*obj_j;

    // хехе костыль
    if(row_end > N){
        row_end = N-1;
    }
    if(col_end > N){
        col_end = N-1;
    }

    // start = omp_get_wtime();
    double  *data1   = (double  *) malloc (M * N * sizeof(double));
    double **F = (double **) malloc (M * sizeof(double *));
    for (int i = 0; i < M; i++)
        F[i] = & (data1[i * N]);
    double  *data2   = (double  *) malloc (M * N * sizeof(double));
    double **a = (double **) malloc (M * sizeof(double *));
    for (int i = 0; i < M; i++)
        a[i] = & (data2[i * N]);    
    double  *data3   = (double  *) malloc (M * N * sizeof(double));
    double **b = (double **) malloc (M * sizeof(double *));
    for (int i = 0; i < M; i++)
        b[i] = & (data3[i * N]);
    double  *data4   = (double  *) malloc (M * N * sizeof(double));
    double **Ar = (double **) malloc (M * sizeof(double *));
    for (int i = 0; i < M; i++)
        Ar[i] = & (data4[i * N]);
    double  *data5   = (double  *) malloc (M * N * sizeof(double));
    double **r = (double **) malloc (M * sizeof(double *));
    for (int i = 0; i < M; i++)
        r[i] = & (data5[i * N]);
    double  *data6   = (double  *) malloc (M * N * sizeof(double));
    double **w = (double **) malloc (M * sizeof(double *));
    for (int i = 0; i < M; i++)
        w[i] = & (data6[i * N]);


    

    initialize(a, b, F, w, r, Ar);
    calc_coefs(h1, h2, EPS, F, a, b, row_start, row_end, col_start, col_end);

    int i, j;
    // if(rank == 0){
    //     for (i = 0; i < M; i++) {
    //         for (j = 0; j < N; j++) {
    //             printf("%f ", b[i][j]);
    //         }
    //         printf("\n");
    //     }
    // }

    printf("Processing submatrix: [%d:%d, %d:%d]\n", row_start, row_end, col_start, col_end);
    int iterations = MRD(h1, h2, a, b, F, w, r, Ar, row_start, row_end, col_start, col_end, obj_i, obj_j, n, k);

    // end = omp_get_wtime();
    // printf("Time spent: %f\n", end - start);

    if(rank == 0) {
        printf("Iterations number: %d\n", iterations);
    }

    // Вывод решения
    // printf("rank %d\n", rank);
    // if(rank == 0){
    //     for (i = 0; i < M; i++) {
    //         for (j = 0; j < N; j++) {
    //             printf("%f ", r[i][j]);
    //         }
    //         printf("\n");
    //     }
    // }

    free (data1);
    free (a);
    free (data2);
    free (b);
    free (data3);
    free (F);
    free (data4);
    free (w);
    free (data5);
    free (r);
    free (data6);
    free (Ar);

    MPI_Finalize();
}

int main(int argc, char *argv[]){

    // if(argc > 1){
    //     int i;
    //     for(i = 1; i < argc; i++){
    //         int num_threads = atoi(argv[i]);
    //         printf("Проводим эксперимент с %d нитями\n", num_threads);
    //         make_experiment(num_threads);
    //     }
    // }
    MPI_Init(&argc, &argv);
    make_experiment(1);

    return 0;
}