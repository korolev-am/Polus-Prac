#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include "mpi.h"

// это задаем область П - прямоугольник, содержащий нашу трапецию
#define A1 -3.5
#define A2 -0.5
#define B1 3.5
#define B2 3.5


#define M 41
#define N 41

#define ACC 1e-5 // это точность метода


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
    } else if (right_line(bottom) > 0 && right_line(top) < 0){
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
    } else if (right_line(left) > 0 && right_line(right) < 0){
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
double Ar[M][N];


int obl_index_to_rank(int i, int j, int n, int k){
    return i * k + j % n;
}


void initialize(double **a, double **b, double **F, double **w, double **r, double **Ar){
    // Инициализация

    int i = 0;
    int j = 0;

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
    for (i = 1; i < M; i++) {
        for (j = 1; j < N; j++) {
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
    double global_norm = 1.0, global_norm_dr, global_norm_r;
    int i, j;
    int iter = 0, msg_iter = 0;
    double max_value;

    double *right_send_buf, *right_recv_buf, *left_send_buf, *left_recv_buf, *top_send_buf, *top_recv_buf, *down_send_buf, *down_recv_buf;
    right_send_buf = (double *) malloc ((row_end - row_start - 1) * sizeof(double));
    right_recv_buf = (double *) malloc ((row_end - row_start - 1) * sizeof(double));
    left_send_buf = (double *) malloc ((row_end - row_start - 1) * sizeof(double));
    left_recv_buf = (double *) malloc ((row_end - row_start - 1) * sizeof(double));
    top_send_buf = (double *) malloc ((col_end - col_start - 1) * sizeof(double));
    top_recv_buf = (double *) malloc ((col_end - col_start - 1) * sizeof(double));
    down_send_buf = (double *) malloc ((col_end - col_start - 1) * sizeof(double));
    down_recv_buf = (double *) malloc ((col_end - col_start - 1) * sizeof(double));

    double** sending_col_data = (double **) malloc (2 * sizeof(double *));
    double** sending_row_data = (double **) malloc (2 * sizeof(double *));
    double** receiving_col_data = (double **) malloc (2 * sizeof(double *));
    double** receiving_row_data = (double **) malloc (2 * sizeof(double *));

    int row_neighbour_ranks[2] = {-1, -1};
    int col_neighbour_ranks[2] = {-1, -1};
    int fixed_row[2] = {-1, -1};
    int fixed_col[2] = {-1, -1};
    int row_n_len = 0;
    int col_n_len = 0;

    MPI_Request r_col[2];
    MPI_Request r_row[2];
    MPI_Request Ar_col[2];
    MPI_Request Ar_row[2];

    if(obj_j + 1 < k){
        // если есть сосед справа
        // будем накапливать для него элементы
        sending_col_data[row_n_len] = right_send_buf;
        receiving_col_data[row_n_len] = right_recv_buf;
        fixed_col[row_n_len] = col_end;
        row_neighbour_ranks[row_n_len++] = obl_index_to_rank(obj_i, obj_j + 1, n, k);
    }
    if(obj_j - 1 >= 0){
        // если есть сосед слева
        sending_col_data[row_n_len] = left_send_buf;
        receiving_col_data[row_n_len] = left_recv_buf;
        fixed_col[row_n_len] = col_start;
        row_neighbour_ranks[row_n_len++] = obl_index_to_rank(obj_i, obj_j - 1, n, k);
    }
    if(obj_i - 1 >= 0){
        // если есть сосед сверху
        sending_row_data[col_n_len] = top_send_buf;
        receiving_row_data[col_n_len] = top_recv_buf;
        fixed_row[col_n_len] = row_start;
        col_neighbour_ranks[col_n_len++] = obl_index_to_rank(obj_i - 1, obj_j, n, k);
    }
    if(obj_i + 1 < n){
        // если есть сосед снизу
        sending_row_data[col_n_len] = down_send_buf;
        receiving_row_data[col_n_len] = down_recv_buf;
        fixed_row[col_n_len] = row_end;
        col_neighbour_ranks[col_n_len++] = obl_index_to_rank(obj_i + 1, obj_j, n, k);
    }

    for(iter = 0; iter < INT_MAX; iter++)
    {

        if (iter % 1000 == 0 && global_norm <= ACC) {
            break;
        }

        norm_r = 0.0;
        norm_dr = 0.0;
        norm = 0.0;
        global_norm = 0.0;
        global_norm_dr = 0.0;
        global_norm_r = 0.0;
        max_value = -1.0;

        int row_n_id, col_n_id;

        for(row_n_id = 0; row_n_id < row_n_len; ++row_n_id) {
            MPI_Irecv(receiving_col_data[row_n_id], row_end - row_start - 1, MPI_DOUBLE, row_neighbour_ranks[row_n_id], msg_iter, MPI_COMM_WORLD, &(r_col[row_n_id]));
        }
        for(col_n_id = 0; col_n_id < col_n_len; ++col_n_id) {
            MPI_Irecv(receiving_row_data[col_n_id], col_end - col_start - 1, MPI_DOUBLE, col_neighbour_ranks[col_n_id], msg_iter, MPI_COMM_WORLD, &(r_row[col_n_id]));
        }

        #pragma omp parallel for collapse(2) reduction(+:norm_r) reduction(max:max_value) private(i, j) shared(a, b, F, r, w)
        for (i = row_start + 1; i < row_end; i++) {
            for (j = col_start + 1; j < col_end; j++) {

                // тут вычисляем невязку
                r[i][j] = -(a[i+1][j] * (w[i+1][j] - w[i][j]) / h1 - a[i][j] * (w[i][j] - w[i-1][j]) / h1) / h1
                        - (b[i][j+1] * (w[i][j+1] - w[i][j]) / h2 - b[i][j] * (w[i][j] - w[i][j-1]) / h2) / h2
                        - F[i][j];
                // printf("%f ", r[i][j]);
                norm_r += r[i][j] * r[i][j] * h1 * h2;

                if(fabs(r[i][j]) > max_value){
                    max_value = fabs(r[i][j]);
                }

                // граничные условия
                top_send_buf[j - col_start - 1] = r[row_start + 1][j];
                down_send_buf[j - col_start - 1] = r[row_end - 1][j];
                left_send_buf[i - row_start - 1] = r[i][col_start + 1];
                right_send_buf[i - row_start - 1] = r[i][col_end - 1];
                
            }
            
        }

        for(row_n_id = 0; row_n_id < row_n_len; ++row_n_id) {
            MPI_Request send_col;
            MPI_Isend(sending_col_data[row_n_id], row_end - row_start - 1, MPI_DOUBLE, row_neighbour_ranks[row_n_id], msg_iter, MPI_COMM_WORLD, &send_col);

            MPI_Wait(&(r_col[row_n_id]), MPI_STATUS_IGNORE);

            for(i = row_start + 1; i < row_end; i++){
                r[i][fixed_col[row_n_id]] = receiving_col_data[row_n_id][i - row_start - 1];
            }

            MPI_Irecv(receiving_col_data[row_n_id], row_end - row_start - 1, MPI_DOUBLE, row_neighbour_ranks[row_n_id], msg_iter+1, MPI_COMM_WORLD, &(Ar_col[row_n_id]));
        }

        for(col_n_id = 0; col_n_id < col_n_len; ++col_n_id) {
            MPI_Request send_row;
            MPI_Isend(sending_row_data[col_n_id], col_end - col_start - 1, MPI_DOUBLE, col_neighbour_ranks[col_n_id], msg_iter, MPI_COMM_WORLD, &send_row);

            MPI_Wait(&(r_row[col_n_id]), MPI_STATUS_IGNORE);

            for(j = col_start + 1; j < col_end; j++){
                r[fixed_row[col_n_id]][j] = receiving_row_data[col_n_id][j - col_start - 1];
            }

            MPI_Irecv(receiving_row_data[col_n_id], col_end - col_start - 1, MPI_DOUBLE, col_neighbour_ranks[col_n_id], msg_iter+1, MPI_COMM_WORLD, &(Ar_row[col_n_id]));
        }

        msg_iter++;

        #pragma omp parallel for collapse(2) reduction(+:norm_dr) private(i, j) shared(a, b, r, Ar)
        for (i = row_start + 1; i < row_end; i++) {
            for (j = col_start + 1; j < col_end; j++) {

                // тут - матрицу A*r
                Ar[i][j] = -(a[i+1][j] * (r[i+1][j] - r[i][j]) / h1 - a[i][j] * (r[i][j] - r[i-1][j]) / h1) / h1
                        - (b[i][j+1] * (r[i][j+1] - r[i][j]) / h2 - b[i][j] * (r[i][j] - r[i][j-1]) / h2) / h2;
                norm_dr += Ar[i][j] * r[i][j] * h1 * h2;
            }
        }

        MPI_Allreduce(&norm_r, &global_norm_r, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&norm_dr, &global_norm_dr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        tau = global_norm_r / global_norm_dr;

        #pragma omp parallel for collapse(2) private(i, j, tmp) shared(w, r)
        for (i = row_start + 1; i < row_end; i++) {
            for (j = col_start + 1; j < col_end; j++) {

                tmp = w[i][j];
                w[i][j] -= tau * r[i][j];
                tmp = w[i][j] - tmp;

                top_send_buf[j - col_start - 1] = w[row_start + 1][j];
                down_send_buf[j - col_start - 1] = w[row_end - 1][j];
                left_send_buf[i - row_start - 1] = w[i][col_start + 1];
                right_send_buf[i - row_start - 1] = w[i][col_end - 1];
            }
        }

        for(row_n_id = 0; row_n_id < row_n_len; ++row_n_id) {
            MPI_Request send_col;
            MPI_Isend(sending_col_data[row_n_id], row_end - row_start - 1, MPI_DOUBLE, row_neighbour_ranks[row_n_id], msg_iter, MPI_COMM_WORLD, &send_col);

            MPI_Wait(&(Ar_col[row_n_id]), MPI_STATUS_IGNORE);

            for(i = row_start + 1; i < row_end; i++){
                w[i][fixed_col[row_n_id]] = receiving_col_data[row_n_id][i - row_start - 1];
            }
        }
        for(col_n_id = 0; col_n_id < col_n_len; ++col_n_id) {
            MPI_Request send_row;
            MPI_Isend(sending_row_data[col_n_id], col_end - col_start - 1, MPI_DOUBLE, col_neighbour_ranks[col_n_id], msg_iter, MPI_COMM_WORLD, &send_row);

            MPI_Wait(&(Ar_row[col_n_id]), MPI_STATUS_IGNORE);

            for(j = col_start + 1; j < col_end; j++){
                w[fixed_row[col_n_id]][j] = receiving_row_data[col_n_id][j - col_start - 1];
            }
        }

        MPI_Allreduce(&max_value, &global_norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        // printf("global_norm: %f, norm_r: %f, norm_dr: %f, tau: %f, iter: %d\n", global_norm, norm_r, norm_dr, tau, iter);

        msg_iter++;
    }
    
    free(right_send_buf);
    free(right_recv_buf);
    free(left_send_buf);
    free(left_recv_buf);
    free(top_send_buf);
    free(top_recv_buf);
    free(down_send_buf);
    free(down_recv_buf);

    free(sending_row_data);
    free(sending_col_data);
    free(receiving_row_data);
    free(receiving_col_data);
    

    return iter;
}


void make_experiment() {

    double h1 = (B1 - A1) / M;
    double h2 = (B2 - A2) / N;
    double EPS = h1 > h2 ? h1*h1 : h2*h2;

    double start = MPI_Wtime();

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

    if(row_end > N) {
        row_end = N-1;
    }
    if(col_end > N) {
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

    printf("Затраченное время %lf", MPI_Wtime()-start);
    MPI_Finalize();

}

int main(int argc, char *argv[]){

    MPI_Init(&argc, &argv);
    make_experiment();

    return 0;
}
