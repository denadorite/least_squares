#include "stdio.h"
#include "math.h"

// Функция gauss принимает на вход размерность матрицы N и саму матрицу mat
void gauss(int N, double mat[N][N+1]) {
    int i, j, k;
    double ratio, pivot;

    // Прямой ход метода Гаусса
    for (i = 0; i < N; i++) {
        pivot = mat[i][i];
        
        for (j = i + 1; j < N; j++) {
            ratio = mat[j][i] / pivot;
            
            for (k = 0; k <= N; k++) {
                mat[j][k] -= ratio * mat[i][k];
            }
        }
    }

    // Обратный ход метода Гаусса
    for (i = N - 1; i >= 0; i--) {
        mat[i][N] /= mat[i][i];
        mat[i][i] = 1;

        for (j = i - 1; j >= 0; j--) {
            mat[j][N] -= mat[j][i] * mat[i][N];
            mat[j][i] = 0;
        }
    }

    // Вывод решений системы уравнений
    if (N == 2)
    {
	printf("a = %.20f\n", mat[0][N]);
	printf("b = %.20f\n", mat[1][N]);
    }

    if (N == 3)
    {
	printf("a0 = %.20f\n", mat[0][N]);
	printf("a1 = %.20f\n", mat[1][N]);
	printf("a2 = %.20f\n", mat[2][N]);
    }
        
    printf("\n");
}

int main() {
    
    double x[] = {0.161, 0.118, 0.926, 0.967, 0.129, 0.765, 0.643, 0.081, 0.182, 0.563};
    double y[] = {3.243, 3.398, 1.287, 0.835, 3.448, 1.497, 1.935, 3.851, 3.103, 2.041};
    int n = sizeof(x) / sizeof(x[0]);

    // Считаем необходимые коэффициенты для составления СЛАУ
    double summ_xi_quad = 0, summ_xi = 0, summ_xi_yi = 0, summ_yi = 0, summ_xi_quad_yi = 0,
  summ_xi_cube = 0, summ_xi_quadro = 0;
    
    for (int i = 0; i < n; i++)
    {
        summ_xi_quad += pow(x[i], 2);
        summ_xi += x[i];
        summ_xi_yi += (x[i] * y[i]);
        summ_yi += y[i];
        summ_xi_quad_yi += pow(x[i], 2) * y[i];
        summ_xi_cube += pow(x[i], 3);
        summ_xi_quadro += pow(x[i], 4);
    }

    printf("xi:\n");

    for (int i = 0; i < n; i++)
    {
	printf("%f, ", x[i]);
    }
    
    printf("\nyi:\n");
    
    for (int i = 0; i < n; i++)
    {
	printf("%f, ", y[i]);
    }
    printf("\n");


    #define N 2
    #define K 3

    double mat[N][N+1] = {
        {summ_xi_quad, summ_xi, summ_xi_yi},
        {summ_xi, n, summ_yi}
    };

    printf("\nСЛАУ 2x2:\n");
    printf("%f*a + %f*b = %f\n%f*a + %d*b = %f\n", summ_xi_quad, summ_xi, summ_xi_yi, summ_xi, n, summ_yi);
    printf("\nРешение СЛАУ 2x2:\n");

    gauss(N, mat);

    // Вычисление погрешностей для линейной зависимости
    double a;
    double b;

    a = mat[0][N];
    b = mat[1][N];

    double fab = 0.0, fab2 = 0.0;

    // Вычисление погрешностей для линейной зависимости
    for (int i = 0; i < n; i++)
    {
        fab += pow((y[i] - (a * x[i] + b)), 2);
    }

    double mat2[K][K+1] = {
        {summ_xi_quad, summ_xi, n, summ_yi},
        {summ_xi_cube, summ_xi_quad, summ_xi, summ_xi_yi},
	{summ_xi_quadro, summ_xi_cube, summ_xi_quad, summ_xi_quad_yi}
    };

    printf("СЛАУ 3x3:\n");
    printf("%f*a0 + %f*a1 + %d*a2 = %f\n%f*a0 + %f*a1 + %f*a2 = %f\n%f*a0 + %f*a1 + %f*a2 = %f\n", summ_xi_quad, summ_xi, n, summ_yi, summ_xi_cube, summ_xi_quad, summ_xi, summ_xi_yi,
	summ_xi_quadro, summ_xi_cube, summ_xi_quad, summ_xi_quad_yi);

    printf("\nРешение СЛАУ 3x3:\n");

    gauss(K, mat2);

    // Вычисление погрешностей для квадратичной зависимости
    double a0, a1, a2;
    a0 = mat2[0][K];
    a1 = mat2[1][K];
    a2 = mat2[2][K];

    for (int i = 0; i < n; i++)
    {
        fab2 += pow(y[i] - (a2 + a1*x[i] + a0*x[i]*x[i]), 2);
    }
    
    printf("Погрешность линейной аппроксимации: %f,\nПогрешность квадратичной аппроксимации: %f.\n", fab, fab2);

    return 0;
}
