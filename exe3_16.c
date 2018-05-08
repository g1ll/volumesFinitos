#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double norma_vetor(double v1[], double v2[], int size);

void jacobi(int n, double **Ap, double x[], double b[], double tol, int ftuv, double tuv);

void gaus_seidel(int n, double **Ap, double x[], double b[], double tol, int ftuv, double tuv);

void sor(int n, double **Ap, double x[], double b[], double tol, double w, int ftuv, double tuv);

void tdma(int n, double **A, double x[], double b[], double tol, int ftuv, double tuv);

int main() {
    unsigned int i, j, //Índices
            n = 10   ,//Número de Volumes
            dtp =1;  //Definir ou não temperatura prescrita no volume a direita
            //(dtp = 1 define uma temperatura prescrita no último volume)

    double l = 3, //Comprimento do domínio computacional (m)
            qw = 10, //FLuxo de calor que entra a direita (W/m²)
            qe = 31, //Fluxo de calor que sai a esquerda (W/m²)
            q = 7, //Fluxo de calor interno (K/m³)
            k = 1, //Condutibilidade do material (K/m²)
            **A, //Matriz de coeficientes do sistema Ax=B
            b[n], //Vetor b do sistema Ax=b
            x[n], //Vetor de soluções com as temperaturas dos volumes
            dx = l / n, //Discretização do domínio computacional
            tol = 1e-60, //Tolerância para os métodos Iterativos
            w = 1, //Fator de relaxação para o método de S.O.R
            tp = 10; //Temperatura prescrita no último volume a direita


    /* aloca as linhas da matriz A */
    A = (double **) calloc(n, sizeof (double *));
    if (A == NULL) {
        printf("** Erro: Memoria Insuficiente **");
        return 0;
    }

    /* aloca as colunas da matriz A */
    for (i = 0; i < n; i++) {
        A[i] = (double*) calloc(n, sizeof (double));
        if (A[i] == NULL) {
            printf("** Erro: Memoria Insuficiente **");
            return 0;
        }
    }

    //Dados do programa
    printf("\nDados do Problema:\n\tdx = %f\n", dx);

    //Montando a matriz de coeficientes
    printf("\n\t Matriz [A b]:\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j && (i == 0 || i == n - 1)) {
                A[i][j] = 1 / dx;
            } else if (i == j) {
                A[i][j] = 2 / dx;
            } else if (j == i + 1 || j == i - 1) {
                A[i][j] = -1 / dx;
            } else {
                A[i][j] = 0;
            }
            printf("\t %-2.2f", A[i][j]);
        }

        //Montando o vetor b
        if (i == 0) {
            //Condiçõa de contorno de flux prescrito a esquerda
            b[i] = q / k * dx + qw;
        } else if (i == n - 1) {
            //Condição de contorno de flux prescrito a direita
            b[i] = q / k * dx - qe;
        } else {
            b[i] = q / k*dx;
        }
        printf("\t %.2f\n", b[i]);
    }

    //INICIALIZNADO O VETOR DE SOLUÇÕES (TEMPERATURAS)
    for (i = 0; i < n; i++) {
        x[i] = 0; //Limpando o vetor de temperaturas
    }

    //RESOLVENDO O SISTEMA POR MÉTODO DE JACOBI:
    jacobi(n, A, x, b, tol, dtp, tp); //tp = Temperatura no ultimo volume  (dtp = 0 = desconhecida)

    //RESOLVENDO O SISTEMA POR GAUS SEIDEL
    for (i = 0; i < n; i++) {
        x[i] = 0;
    }//Limpando o vetor de temperaturas para o próximo método
    gaus_seidel(n, A, x, b, tol, dtp, tp); //tp = Temperatura no ultimo volume  (dtp = 0 = desconhecida)

    //RESOLVENDO O SISTEMA POR S.O.R
    for (i = 0; i < n; i++) {
        x[i] = 0;
    }//Limpando o vetor de temperaturas para o próximo método
    sor(n, A, x, b, tol, w, dtp, tp); //tp = Temperatura no ultimo volume  (dtp = 0 = desconhecida)

    //RESOLVENDO O SISTEMA POR TDMA
    for (i = 0; i < n; i++) {
        x[i] = 0;
    }//Limpando o vetor de temperaturas para o próximo método
    tdma(n, A, x, b, tol, dtp, tp); //tp = Temperatura no ultimo volume  (dtp = 0 = desconhecida)

    free(A);
}

double norma_vetor(double v1[], double v2[], int size) {
    double norma = 0;
    for (int i = 0; i < size; i++) {
        norma += pow(v1[i] - v2[i], 2);
    }
    norma = sqrt(norma);
    //printf("\n%f", norma);
    return norma;
}

void jacobi(int n, double **Ap, double x[], double b[], double tol, int ftuv, double tuv) {
    int i, j, iter = 0;
    double sigma = 0, xa[n];

    printf("\nMÉTODO DE JACOBI:\nMatriz [A b]\n");
    printf("iterando..");
    do {
        for (i = 0; i < n; i++) {
            xa[i] = x[i];
        }

        for (i = 0; i < n; i++) {
            sigma = 0;
            for (j = 0; j < n; j++) {
                if (j != i) {
                    sigma += Ap[i][j] * x[j];
                }
            }
            x[i] = (b[i] - sigma) / Ap[i][i];
        }
        if (ftuv != 0) {
            x[n - 1] = tuv; //forçando a temperatura no último volume
        }
        iter++;
    } while (norma_vetor(xa, x, n) > tol && iter < 100000 * n);

    if (iter < 100000 * n) {
        printf("\n\t%.10f <= %.7f Tolerância\n", norma_vetor(xa, x, n), tol);
        printf("\nNúmero de iterações: %d\n", iter);
        printf("\nSOLUÇÃO x:\n");
        for (i = 0; i < n; i++) {
            printf("\t %.3f", x[i]);
        }
        printf("\n\n");
    } else {
        printf("\nIter: %d\nNÃO CONVERGIU!!!\n", iter);
    }
}

void gaus_seidel(int n, double **Ap, double x[], double b[], double tol, int ftuv, double tuv) {
    int i, j, iter = 0;
    double sigma = 0, xa[n];

    printf("\nMÉTODO DE GAUS-SEIDEL:\nMatriz [A b]\n");
    printf("iterando..");
    do {
        for (i = 0; i < n; i++) {
            xa[i] = x[i];
        }
        for (i = 0; i < n; i++) {
            sigma = 0;
            for (j = 0; j < i; j++) {
                sigma += Ap[i][j] * x[j];
            }
            for (j = i + 1; j < n; j++) {
                sigma += Ap[i][j] * xa[j];
            }
            x[i] = (b[i] - sigma) / Ap[i][i];
        }
        if (ftuv != 0) {
            x[n - 1] = tuv; //forçando a temperatura no último volume
        }
        iter++;
    } while (norma_vetor(xa, x, n) > tol && iter < 100000 * n);

    if (iter < 100000 * n) {
        printf("\n\t%.10f <= %.7f Tolerância\n", norma_vetor(xa, x, n), tol);
        printf("\nNúmero de iterações: %d\n", iter);
        printf("\nSOLUÇÃO x:\n");
        for (i = 0; i < n; i++) {
            printf("\t %.3f", x[i]);
        }
        printf("\n\n");
    } else {
        printf("\nIter: %d\nNÃO CONVERGIU!!!\n", iter);
    }
}

void sor(int n, double **Ap, double x[], double b[], double tol, double w, int ftuv, double tuv) {
    int i, j, iter = 0;
    double sigma = 0, xa[n];

    printf("\nMÉTODO DE S.O.R:\nMatriz [A b]\nFator w: %.3f\n\n", w);
    printf("iterando..");
    do {
        for (i = 0; i < n; i++) {
            xa[i] = x[i];
        }

        for (i = 0; i < n; i++) {
            sigma = 0;
            for (j = 0; j < i; j++) {
                sigma += Ap[i][j] * x[j];
            }

            for (j = i + 1; j < n; j++) {
                sigma += Ap[i][j] * xa[j];
            }
            x[i] = (b[i] - sigma) / Ap[i][i];
            x[i] = w * x[i]+(1 - w) * xa[i];
        }
        if (ftuv != 0) {
            x[n - 1] = tuv; //forçando a temperatura no último volume
        }
        iter++;
    } while (norma_vetor(xa, x, n) > tol && iter < 100000 * n);

    if (iter < 100000 * n) {
        printf("\n\t%.10f <= %.7f Tolerância\n", norma_vetor(xa, x, n), tol);
        printf("\nNúmero de iterações: %d\n", iter);
        printf("\nSOLUÇÃO x:\n");
        for (i = 0; i < n; i++) {
            printf("\t %.3f", x[i]);
        }
        printf("\n\n");
    } else {
        printf("\nIter: %d\nNÃO CONVERGIU!!!\n", iter);
    }
}

void tdma(int n, double **A, double x[], double b[], double tol, int ftuv, double tuv) {
    int i, iter = 0; //indices
    double p[n], q[n], xa[n];

    printf("\nMÉTODO TDMA:\nMatriz [A b]\n\n");
    printf("iterando..");
    do {
        for (i = 0; i < n; i++) {
            xa[i] = x[i];
        }
        p[0] = A[0][1] / A[0][0];
        q[0] = b[0] / A[0][0];
        //printf("\n p0 = %f \tq0 = %f", p[0], q[0]);
        for (i = 1; i < n; i++) {
            if (i < n - 1) {
                p[i] = A[i][i + 1] / (A[i][i] + A[i][i - 1] * p[i - 1]);
                q[i] = (b[i] + A[i][i - 1] * q[i - 1]) / (A[i][i] - A[i][i - 1] * p[i - 1]);
            } else {
                p[i] = 0;
                q[i] = (b[i] + A[i][i - 1] * q[i - 1]) / (A[i][i] - A[i][i - 1] * p[i - 1]);
                x[i] = q[i];
                if (ftuv != 0) {
                    x[n - 1] = tuv; //forçando a temperatura no último volume
                }

            }
            //printf("\n p%d = %f \tq%d = %f", i, p[i], i, q[i]);
        }
        //printf("\n\n Ultimo nó: %f\n", x[9]);
        //printf("\nx[%d]=%f", n - 1, x[n - 1]);
        for (i = n - 2; i >= 0; i--) {
            x[i] = p[i] * x[i + 1] + q[i];
            //printf("\nx[%d]=%f x[%d]=%f", i + 1, x[i + 1], i, x[i]);
        }
        iter++;
    } while (norma_vetor(xa, x, n) > tol && iter < 100000 * n);

    if (iter < 100000 * n) {
        printf("\n\t%.10f <= %.7f Tolerância\n", norma_vetor(xa, x, n), tol);
        printf("\nNúmero de iterações: %d\n", iter);
        printf("\nSOLUÇÃO x:\n");
        for (i = 0; i < n; i++) {
            printf("\t %.3f", x[i]);
        }
        printf("\n\n");
    } else {
        printf("\nIter: %d\nNÃO CONVERGIU!!!\n", iter);
    }
}

