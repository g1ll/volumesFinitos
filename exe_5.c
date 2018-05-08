
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265

double *gaus_seidel(int n, double **Ap, double x[], double b[], double tol);

double norma_vetor(double v1[], double v2[], int size);

int main(int argc, char** argv) {

    int i, j, c = 0, //Índices e contador
            nx = ((argc > 1) ? atoi(argv[1]) : 3),
            ny = ((argc > 2) ? atoi(argv[2]) : nx);
              
    double a = 1, //Comprimento da placa
            b = 1, //Altura da placa
            dx = a / nx, //Delta x
            dy = b / ny, //Delta y
            ro = 1, //Densidade ro
            ga = 1, //Gama = k/Cp (condutividade/calor específico)
            u = ((argc > 3) ? atof(argv[3]) : 0), //Velocidade em x
            v = ((argc > 4) ? atof(argv[4]) : 0), //Velocidade em y
            pex, pey, //Peclet em x e y
            ax, ay, bx, by, //Alfa e Beta em x e y
            Ap[ny][nx], //Matriz de coeficientes Ap
            Ae[ny][nx], //Matriz de coeficientes Ae
            Aw[ny][nx], //Matriz de coeficientes Aw
            An[ny][nx], //Matriz de coeficientes An
            As[ny][nx], //Matriz de coeficientes As
            **A, //Matriz de coeficientes do sistema Ax=b
            B[nx * ny], //Vetor b do sistema Ax=b (fontes)
            *T, //Vetor com Temperaturas
            tol = 1e-80; //Tolerância para os métodos Iterativos

    printf("\n\nDados do Problema:");
    printf("\n\n\t nx: %d \n\t ny: %d", nx, ny);
    printf("\n\n\t dx: %0.5f \n\t dy: %0.5f", dx, dy);
    printf("\n\t u: %0.5f \n\t v: %0.5f", u, v);
    printf("\n\n");
    //CALCULANDO PECLET, ALFA E BETA PARA X E Y (WUDS)
    pex = (u * ro * dx) / ga; //Peclet em x
    pey = (v * ro * dy) / ga; //Peclet em y
    ax = pow(pex, 2) / (10 + 2 * pow(pex, 2)); //Calculando alfa em x
    ay = pow(pey, 2) / (10 + 2 * pow(pey, 2)); //Calculando alfa em y
    bx = (1 + 0.005 * pow(pex, 2)) / 1 + 0.05 * pow(pex, 2); //Cálculo de beta em x       
    by = (1 + 0.005 * pow(pey, 2)) / 1 + 0.05 * pow(pey, 2); //Cálculo de beta em y


    //CALCULANDO COEFICIENTES Ae, Aw, An, As e Ap PARA CADA PONTO DA MALHA
    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) {

            if (i == 0 && j == 0) { //CANTO SUPERIOR ESQUERDO VOLUME 1
                printf("C1");
                Ae[i][j] = -ro * u * (1 / 2 + ax) + bx * ga / dx;
                Aw[i][j] = 0;
                An[i][j] = 0;
                As[i][j] = ro * v * (1 / 2 - ay) + (by * ga) / dy;
                Ap[i][j] = bx * ga / dx + by * ga / dy + 2 * bx * ga / dx + 2 * by * ga / dy;
                B[c] = 0;
            } else if (i == 0 && j == nx - 1) {//CANTO SUP DIR VOLUME NX
                printf("C2");
                Ae[i][j] = 0;
                Aw[i][j] = ro * u * (1 / 2 - ax) + bx * ga / dx;
                An[i][j] = 0;
                As[i][j] = ro * v * (1 / 2 - ay) + by * ga / dy;
                Ap[i][j] = bx * ga / dx + by * ga / dy + 2 * bx * ga / dx + 2 * by * ga / dy;
                B[c] = 0;
            } else if (i == ny - 1 && j == 0) {//CANTO INFERIOR ESQUERDO
                printf("C3");
                Ae[i][j] = -ro * u * (1 / 2 + ax) + bx * ga / dx;
                Aw[i][j] = 0;
                An[i][j] = -ro * v * (1 / 2 + ay) + by * ga / dy;
                As[i][j] = 0;
                Ap[i][j] = bx * ga / dx + by * ga / dy + 2 * bx * ga / dx + 2 * by * ga / dy;
                B[c] = ro * v * sin(PI * (j * dx + dx / 2) / a) + 2 / dy * sin(PI * (j * dx + dx / 2) / a);
            } else if (i == ny - 1 && j == nx - 1) { //CANTO INFERIOR DIREITO
                printf("C4");
                Ae[i][j] = 0;
                Aw[i][j] = ro * u * (1 / 2 - ax) + bx * ga / dx;
                An[i][j] = -ro * v * (1 / 2 + ay) + by * ga / dy;
                As[i][j] = 0;
                Ap[i][j] = bx * ga / dx + by * ga / dy + 2 * bx * ga / dx + 2 * by * ga / dy;
                B[c] = ro * v * sin(PI * (j * dx + dx / 2) / a) + 2 / dy * sin(PI * (j * dx + dx / 2) / a);
            } else if (j == 0 && i > 0 && i < ny - 1) { //LATERAL ESQUERDA
                printf("LE");
                Ae[i][j] = -ro * u * (1 / 2 + ax) + bx * ga / dx;
                Aw[i][j] = 0;
                An[i][j] = -ro * v * (1 / 2 + ay) + by * ga / dy;
                As[i][j] = ro * v * (1 / 2 - ay) + by * ga / dy;
                Ap[i][j] = bx * ga / dx + by * ga / dy + 2 * bx * ga / dx + by * ga / dy;
                B[c] = 0;
            } else if (j == nx - 1 && i > 0 && i < ny - 1) { //LATERAL DIREITA
                printf("LD");
                Ae[i][j] = 0;
                Aw[i][j] = ro * u * (1 / 2 - ax) + bx * ga / dx;
                An[i][j] = -ro * v * (1 / 2 + ay) + by * ga / dy;
                As[i][j] = ro * v * (1 / 2 - ay) + by * ga / dy;
                Ap[i][j] = bx * ga / dx + by * ga / dy + 2 * bx * ga / dx + by * ga / dy;
                B[c] = 0;
            } else if (i == 0 && j > 0 && j < nx - 1) { //TOPO
                printf("TP");
                Ae[i][j] = -ro * u * (1 / 2 + ax) + bx * ga / dx;
                Aw[i][j] = ro * u * (1 / 2 - ax) + bx * ga / dx;
                An[i][j] = 0;
                As[i][j] = ro * v * (1 / 2 - ay) + by * ga / dy;
                Ap[i][j] = bx * ga / dx + by * ga / dy + 2 * bx * ga / dx + by * ga / dy;
                B[c] = 0;
            } else if (i == ny - 1 && j > 0 && j < nx - 1) { //INFERIOR
                printf("IN");
                Ae[i][j] = -ro * u * (1 / 2 + ax) + bx * ga / dx;
                Aw[i][j] = ro * u * (1 / 2 - ax) + bx * ga / dx;
                An[i][j] = -ro * v * (1 / 2 + ay) + by * ga / dy;
                As[i][j] = 0;
                Ap[i][j] = bx * ga / dx + by * ga / dy + 2 * bx * ga / dx + 2 * by * ga / dy;
                B[c] = ro * v * sin(PI * (j * dx + dx / 2) / a) + 2 / dy * sin(PI * (j * dx + dx / 2) / a);
            } else { //VOLUMES DO CENTRO
                printf("CE");
                Ae[i][j] = -ro * u * (1 / 2 + ax) + bx * ga / dx;
                Aw[i][j] = ro * u * (1 / 2 - ax) + bx * ga / dx;
                An[i][j] = -ro * v * (1 / 2 + ay) + by * ga / dy;
                As[i][j] = ro * v * (1 / 2 - ay) + by * ga / dy;
                Ap[i][j] = bx * ga / dx + by * ga / dy + bx * ga / dx + by * ga / dy;
            }
            Ap[i][j] += Ae[i][j] + Aw[i][j] + An[i][j] + As[i][j];
            c++;
            printf(" ");
        }
        printf("\n");
    }

    //IMPRIMINDO OS COEFICIENTES DE CADA PONTO
    printf("\nAe:\n");
    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) {
            printf("%0.3f ", Ae[i][j]);
        }
        printf("\n");
    }
    printf("\n----------\n");

    printf("\nAw:\n");
    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) {
            printf("%0.3f ", Aw[i][j]);
        }
        printf("\n");
    }
    printf("\n----------\n");

    printf("\nAn:\n");
    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) {
            printf("%0.3f ", An[i][j]);
        }
        printf("\n");
    }
    printf("\n----------\n");

    printf("\nAs:\n");
    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) {
            printf("%0.3f ", As[i][j]);
        }
        printf("\n");
    }
    printf("\n----------\n");

    printf("\nAp:\n");
    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) {
            printf("%0.3f ", Ap[i][j]);
        }
        printf("\n");
    }
    printf("\n----------\n");

    //IMPRIMINDO O VETOR B
    printf("\nB\n");
    for (j = 0; j < nx * ny; j++) {
        printf("%0.3f ", B[j]);
    }
    printf("\n--------------\n\n\n");

    //MONTANDO A MATRIZ DE COEFICIENTES

    /* aloca as linhas da matriz A */
    A = (double **) calloc(ny*nx, sizeof (double *));
    if (A == NULL) {
        printf("** Erro: Memoria Insuficiente **");
        return 0;
    }

    /* aloca as colunas da matriz A */
    for (i = 0; i < ny * nx; i++) {
        A[i] = (double*) calloc(ny*nx, sizeof (double));
        if (A[i] == NULL) {
            printf("** Erro: Memoria Insuficiente **");
            return 0;
        }
    }

    //Montando a matriz A com coeficientes
    c = 0;
    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) {
            A[c][c] = Ap[i][j];
            if (c - 1 >= 0) {
                A[c][c - 1] = -Aw[i][j];
            }
            if (c + 1 < nx * ny) {
                A[c][c + 1] = -Ae[i][j];
            }
            if (c - 3 >= 0) {
                A[c][c - 3] = -An[i][j];
            }
            if (c + 3 < nx * ny) {
                A[c][c + 3] = -As[i][j];
            }
            c++;
        }
    }

    //Imprimindo a matriz A b com coeficientes
    printf("MATRIZ DE COEFICIENTES A e vetor B\n\nA=[\n");
    for (i = 0; i < ny * nx; i++) {
        for (j = 0; j < ny * nx; j++) {
            printf(" %.2f", A[i][j]);
        }

        if (i == ny * nx - 1) {
            printf("]\n");
        } else {
            printf(";\n");
        }
    }
    printf("\nb=[ ");
    for (i = 0; i < nx * ny; i++) {
        printf("%.2f ", B[i]);
    }
    printf(" ]'");


    //RESOLVENDO O SISTEMA LINEAR
    printf(" \n\n");
    T = (double *) calloc(ny*nx, sizeof (double *));
    if (T == NULL) {
        printf("** Erro: Memoria Insuficiente **");
        return 0;
    }

    for (i = 0; i < ny * nx; i++) {
        printf(" %.0f", T[i]);
        if ((i + 1) % nx == 0) {
            printf(" \n");
        }
    }

    T = gaus_seidel(ny*nx, A, T, B, tol);

    printf("Dados do Problema:");
    printf("\n\n\t nx: %d \n\t ny: %d", nx, ny);
    printf("\n\n\t dx: %0.5f \n\t dy: %0.5f", dx, dy);
    printf("\n\t u: %0.5f \n\t v: %0.5f", u, v);
    printf("\n\t Pex: %0.5f \n\t Pey: %0.5f", pex, pey);
    printf("\n\t ax: %0.5f \n\t ay: %0.5f", ax, ay);
    printf("\n\t bx: %0.5f \n\t by: %0.5f", bx, by);
    printf("\n\n");
    printf("\n\nResultado:\n\n");
    //IMPRINDO O RESULTADO DO VETOR DE TEMPERATURAS
    for (i = 0; i < ny * nx; i++) {
        printf("  %f ", T[i]);
        if ((i + 1) % nx == 0) {
            printf(" \n");
        }
    }
    printf(" \n\n");
    return (EXIT_SUCCESS);
}

double *gaus_seidel(int n, double **Ap, double x[], double b[], double tol) {
    int i, j, iter = 0;
    double sigma = 0, xa[n];

    printf("\nMÉTODO DE GAUS-SEIDEL:\nMatriz [A b]\n");
    for (i = 0; i < n; i++) {
        xa[i] = x[i] = 0;
    }
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
        iter++;
    } while (norma_vetor(xa, x, n) > tol && iter < 1000000 * n);

    if (iter < 100000 * n) {
        printf("\n\t%.10f <= %.7f Tolerância\n", norma_vetor(xa, x, n), tol);
        printf("\nNúmero de iterações: %d\n", iter);
        printf("\n\n");
    } else {
        printf("\nIter: %d\nNÃO CONVERGIU!!!\n", iter);
    }
    return x;
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