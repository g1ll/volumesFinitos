/**
 * REFERENCES:
 * https://en.wikipedia.org/wiki/C_date_and_time_functions
 **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#define PI 3.14159265

double *gaus_seidel(unsigned int n, double **Ap, double x[], double b[], double tol);

double norma_vetor(double v1[], double v2[], int size);

void make_directory(const char* name);

int main(int argc, char** argv) {

    FILE *fd, *fl; //Gravar dados dos resultados
    time_t current_time;
    struct tm *tm_info;
    char cmd[200], //String comando para gnuplot
            c_time_string[100],
            dirname[100],
            filename[100],
            logfilename[100];

    int i, j, c, //Índices e contador
            dbl = ((argc > 1) ? atoi(argv[1]) : 0),
            dbo = ((argc > 2) ? atoi(argv[2]) : 0),
            nx = ((argc > 3) ? atoi(argv[3]) : 3),
            ny = ((argc > 4) ? atoi(argv[4]) : nx);

    double a = 1, //Comprimento da placa
            b = 1, //Altura da placa
            dx = a / nx, //Delta x
            dy = b / ny, //Delta y
            u = ((argc > 5) ? atof(argv[5]) : 0), //Velocidade em x
            v = ((argc > 6) ? atof(argv[6]) : 0), //Velocidade em y
            q = ((argc > 7) ? atof(argv[7]) : 0), //Geração de calor
            ro = ((argc > 8) ? atof(argv[8]) : 1), //Densidade ro
            ga = ((argc > 9) ? atof(argv[9]) : 1), //Gama = k/Cp (condutividade/calor específico)
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
            Tp[ny][nx], //Campo de temperaturas solução numérica
            Ta[ny][nx], //Campo de temperaturas solução analítica
            tmax = -9999, //Temperatura máxima
            tol = 1e-10; //Tolerância para os métodos Iterativos

    /*DATA ATUAL*/
    current_time = time(NULL);

    if (current_time == ((time_t) - 1)) {
        printf("Erro na data.\n");
        sprintf(dirname, "resultado");
    } else {
        tm_info = localtime(&current_time);
        /* Convert to local time format. */
        //c_time_string = ctime(&current_time);
        strftime(c_time_string, sizeof c_time_string, "%F_%H%M%S", tm_info);
        if (dbo)printf("\n\tDATA: %s.\n", c_time_string);
        if (c_time_string == NULL) {
            printf("Erro na data.\n");
            fprintf(fl, "Erro na data.\n");
            sprintf(dirname, "resultado");
        } else {
            sprintf(dirname, "resultado_%s", c_time_string);
        }
    }

    make_directory(dirname);
    sprintf(logfilename, "%s/exe5_log.txt", dirname);
    fl = fopen(logfilename, "w");
    if (fl == NULL) //if file does not exist, create it
    {
        fl = fopen(logfilename, "wb");
    }
    if (fl == NULL) {
        return 0;
        dbl = 0;
    }

    //CALCULANDO PECLET, ALFA E BETA PARA X E Y (WUDS)
    pex = (u * ro * dx) / ga; //Peclet em x
    pey = (v * ro * dy) / ga; //Peclet em y
    ax = pow(pex, 2) / (10 + 2 * pow(pex, 2)); //Calculando alfa em x
    ay = pow(pey, 2) / (10 + 2 * pow(pey, 2)); //Calculando alfa em y
    bx = (1 + 0.005 * pow(pex, 2)) / 1 + 0.05 * pow(pex, 2); //Cálculo de beta em x       
    by = (1 + 0.005 * pow(pey, 2)) / 1 + 0.05 * pow(pey, 2); //Cálculo de beta em y



    printf("Dados do Problema:");
    printf("\n\n\t nx: %d \n\t ny: %d", nx, ny);
    printf("\n\n\t dx: %0.5f \n\t dy: %0.5f", dx, dy);
    printf("\n\n\t ro: %0.5f \n\t gama: %0.5f", ro, ga);
    printf("\n\t q: %0.5f", q);
    printf("\n\t u: %0.5f \n\t v: %0.5f", u, v);
    printf("\n\t Pex: %0.5f \n\t Pey: %0.5f", pex, pey);
    printf("\n\t ax: %0.5f \n\t ay: %0.5f", ax, ay);
    printf("\n\t bx: %0.5f \n\t by: %0.5f", bx, by);
    printf("\n\n");
    if (dbo)printf("\n\nResultado:\n\n");
    //SAÍDA DO LOG
    fprintf(fl, "Dados do Problema:");
    fprintf(fl, "\n\n\t nx: %d \n\t ny: %d", nx, ny);
    fprintf(fl, "\n\n\t dx: %0.5f \n\t dy: %0.5f", dx, dy);
    fprintf(fl, "\n\n\t ro: %0.5f \n\t gama: %0.5f", ro, ga);
    fprintf(fl, "\n\t q: %0.5f", q);
    fprintf(fl, "\n\t u: %0.5f \n\t v: %0.5f", u, v);
    fprintf(fl, "\n\t Pex: %0.5f \n\t Pey: %0.5f", pex, pey);
    fprintf(fl, "\n\t ax: %0.5f \n\t ay: %0.5f", ax, ay);
    fprintf(fl, "\n\t bx: %0.5f \n\t by: %0.5f", bx, by);
    fprintf(fl, "\n\n");

    //TERMOS FONTE
    for(i=0;i<ny*nx;i++){
        B[i]=q;
    }
    //CANTO SUPERIOR ESQUERDO C1
    i = 0;
    j = 0;
    Ae[i][j] = -(1 / 2 - ax) * ro * u + bx * ga / dx;
    Aw[i][j] = 0;
    An[i][j] = 0;
    As[i][j] = ro * v * (1 / 2 + ay) + (by * ga) / dy;
    Ap[i][j] = (2 * bx * ga) / dx + (2 * by * ga) / dy;
    
    //CANTO SUPERIOR DIREITO C2
    i = 0;
    j = nx - 1;
    Ae[i][j] = 0;
    Aw[i][j] = ro * u * (1 / 2 - ax) + bx * ga / dx;
    An[i][j] = 0;
    As[i][j] = ro * v * (1 / 2 + ay) + by * ga / dy;
    Ap[i][j] = 2 * bx * ga / dx + 2 * by * ga / dy;
    
    //CANTO INFERIOR ESQUERDO C3
    i = ny - 1;
    j = 0;
    Ae[i][j] = -(1 / 2 - ax) * ro * u + bx * ga / dx;
    Aw[i][j] = 0;
    An[i][j] = -(1 / 2 - ay) * ro * v + by * ga / dy;
    As[i][j] = 0;
    Ap[i][j] = 2 * bx * ga / dx + 2 * by * ga / dy;
    B[0] = q + ro * v * sin(PI * (j * dx + dx / 2) / a) + 2 / dy * sin(PI * (j * dx + dx / 2) / a);
    
    //CANTO INFERIOR DIREITO
    i = ny - 1;
    j = nx - 1;
    Ae[i][j] = 0;
    Aw[i][j] = ro * u * (1 / 2 - ax) + bx * ga / dx;
    An[i][j] = -(1 / 2 - ay) * ro * v + by * ga / dy;
    As[i][j] = 0;
    Ap[i][j] = 2 * bx * ga / dx + 2 * by * ga / dy;
    B[j] = q + ro * v * sin(PI * (j * dx + dx / 2) / a) + 2 / dy * sin(PI * (j * dx + dx / 2) / a);
    
    
    for(i=1;i<nx-1;i++){
                B[i] = q + ro * v * sin(PI * (i * dx + dx / 2) / a) + 2 / dy * sin(PI * (i * dx + dx / 2) / a);
    }
    
    //CALCULANDO COEFICIENTES Ae, Aw, An, As e Ap PARA CADA PONTO DA MALHA
    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) {
            if(j>=1 && j<=nx-2 && i>=1 && i<=ny-2){//VOLUMES DO CENTRO
                Ae[i][j] = -(1 / 2 - ax) * ro * u + bx * ga / dx;
                Aw[i][j] = ro * u * (1 / 2 - ax) + bx * ga / dx;
                An[i][j] = -(1 / 2 - ay) * ro * v + by * ga / dy;
                As[i][j] = ro * v * (1 / 2 + ay) + by * ga / dy;
                Ap[i][j] = 0;
            } else if (j == 0 && i > 0 && i < ny - 1) { //LATERAL ESQUERDA

                Ae[i][j] = -(1 / 2 - ax) * ro * u + bx * ga / dx;
                Aw[i][j] = 0;
                An[i][j] = -(1 / 2 - ay) * ro * v + by * ga / dy;
                As[i][j] = ro * v * (1 / 2 + ay) + by * ga / dy;
                Ap[i][j] = 2 * bx * ga / dx;
            } else if (j == nx - 1 && i > 0 && i < ny - 1) { //LATERAL DIREITA

                Ae[i][j] = 0;
                Aw[i][j] = ro * u * (1 / 2 - ax) + bx * ga / dx;
                An[i][j] = -(1 / 2 - ay) * ro * v + by * ga / dy;
                As[i][j] = ro * v * (1 / 2 + ay) + by * ga / dy;
                Ap[i][j] = 2 * bx * ga / dx;
            } else if (i == 0 && j > 0 && j < nx - 1) { //VOLUMES NO TOPO

                Ae[i][j] = -(1 / 2 - ax) * ro * u + bx * ga / dx;
                Aw[i][j] = ro * u * (1 / 2 - ax) + bx * ga / dx;
                An[i][j] = 0;
                As[i][j] = ro * v * (1 / 2 + ay) + by * ga / dy;
                Ap[i][j] = 2 * by * ga / dy;
            } else if (i == ny - 1 && j > 0 && j < nx - 1) { //VOLUMES DA PARTE INFERIOR
                Ae[i][j] = -(1 / 2 - ax) * ro * u + bx * ga / dx;
                Aw[i][j] = ro * u * (1 / 2 - ax) + bx * ga / dx;
                An[i][j] = -(1 / 2 - ay) * ro * v + by * ga / dy;
                As[i][j] = 0;
                Ap[i][j] = 2 * by * ga / dy; 
            } 
            Ap[i][j] += Ae[i][j] + Aw[i][j] + An[i][j] + As[i][j];
            if (dbo)printf(" ");
            if (dbl)fprintf(fl, " ");
        }
        if (dbo)printf("\n");
        if (dbl)fprintf(fl, "\n");
    }

    //IMPRIMINDO OS COEFICIENTES DE CADA PONTO
    if (dbo)printf("\nAe:\n");
    if (dbl)fprintf(fl, "\nAe:\n");
    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) {
            if (dbo)printf("%0.3f ", Ae[i][j]);
            if (dbl)fprintf(fl, "%0.3f ", Ae[i][j]);
        }
        if (dbo)printf("\n");
        if (dbl)fprintf(fl, "\n");
    }
    if (dbo)printf("\n----------\n");
    if (dbl)fprintf(fl, "\n----------\n");

    if (dbo)printf("\nAw:\n");
    if (dbl)fprintf(fl, "\nAw:\n");
    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) {
            if (dbo)printf("%0.3f ", Aw[i][j]);
            if (dbl)fprintf(fl, "%0.3f ", Aw[i][j]);
        }
        if (dbo)printf("\n");
        if (dbl)fprintf(fl, "\n");
    }
    if (dbo)printf("\n----------\n");
    if (dbl)fprintf(fl, "\n----------\n");

    if (dbo)printf("\nAn:\n");
    if (dbl)fprintf(fl, "\nAn:\n");
    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) {
            if (dbo)printf("%0.3f ", An[i][j]);
            if (dbl)fprintf(fl, "%0.3f ", An[i][j]);
        }
        if (dbo)printf("\n");
        if (dbl)fprintf(fl, "\n");
    }
    if (dbo)printf("\n----------\n");
    if (dbl)fprintf(fl, "\n----------\n");

    if (dbo)printf("\nAs:\n");
    if (dbl)fprintf(fl, "\nAs:\n");
    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) {
            if (dbo)printf("%0.3f ", As[i][j]);
            if (dbl)fprintf(fl, "%0.3f ", As[i][j]);
        }
        if (dbo)printf("\n");
        if (dbl)fprintf(fl, "\n");
    }
    if (dbo)printf("\n----------\n");
    if (dbl)fprintf(fl, "\n----------\n");

    if (dbo)printf("\nAp:\n");
    if (dbl)fprintf(fl, "\nAp:\n");

    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) {
            if (dbo)printf("%0.3f ", Ap[i][j]);
            if (dbl)fprintf(fl, "%0.3f ", Ap[i][j]);
        }
        if (dbo)printf("\n");
        if (dbl)fprintf(fl, "\n");
    }
    if (dbo)printf("\n----------\n");
    if (dbl)fprintf(fl, "\n----------\n");

    //IMPRIMINDO O VETOR B
    if (dbo)printf("\nB\n");
    if (dbl)fprintf(fl, "\nB\n");
    for (j = 0; j < nx * ny; j++) {
        if (dbo)printf("%0.3f ", B[j]);
        if (dbl)fprintf(fl, "%0.3f ", B[j]);
    }
    if (dbo)printf("\n--------------\n\n\n");
    if (dbl)fprintf(fl, "\n--------------\n\n\n");

    //MONTANDO A MATRIZ DE COEFICIENTES
    /* aloca as linhas da matriz A */
    A = (double **) calloc(ny*nx, sizeof (double *));
    if (A == NULL) {
        if (dbo)printf("** Erro: Memoria Insuficiente **");
        if (dbl)fprintf(fl, "** Erro: Memoria Insuficiente **");
        return 0;
    }

    /* aloca as colunas da matriz A */
    for (i = 0; i < ny * nx; i++) {
        A[i] = (double*) calloc(ny*nx, sizeof (double));
        if (A[i] == NULL) {
            if (dbo)printf("** Erro: Memoria Insuficiente **");
            if (dbl)fprintf(fl, "** Erro: Memoria Insuficiente **");
            return 0;
        }
    }

    //Montando a matriz A com coeficientes
    for (i = ny - 1, c = 0; i >= 0; i--) {
        for (j = 0; j < nx; j++, c++) {
            A[c][c] = Ap[i][j];
            if (c - 1 >= 0) {
                A[c][c - 1] = -Aw[i][j];
            }
            if (c + 1 < nx * ny) {
                A[c][c + 1] = -Ae[i][j];
            }
            if (c - nx >= 0) {
                A[c][c - nx] = -As[i][j];
            }
            if (c + nx < nx * ny) {
                A[c][c + nx] = -An[i][j];
            }
        }
    }

    //Imprimindo a matriz A b com coeficientes
    if (dbo)printf("MATRIZ DE COEFICIENTES A e vetor B\n\nA=[\n");
    if (dbl)fprintf(fl, "MATRIZ DE COEFICIENTES A e vetor B\n\nA=[\n");
    for (i = 0; i < ny * nx; i++) {
        for (j = 0; j < ny * nx; j++) {
            if (dbo)printf(" %.2f", A[i][j]);
            if (dbl)fprintf(fl, " %.2f", A[i][j]);
        }
        if (i == ny * nx - 1) {
            if (dbo)printf("]\n");
            if (dbl)fprintf(fl, "]\n");
        } else {
            if (dbo)printf(";\n");
            if (dbl)fprintf(fl, ";\n");
        }
    }
    if (dbo)printf("\nb=[ ");
    if (dbl)fprintf(fl, "\nb=[ ");
    for (i = 0; i < nx * ny; i++) {
        if (dbo)printf("%.2f ", B[i]);
        if (dbl)fprintf(fl, "%.2f ", B[i]);
    }
    if (dbo)printf(" ]'");
    if (dbl)fprintf(fl, " ]'");

    //RESOLVENDO O SISTEMA LINEAR
    if (dbo)printf(" \n\n");
    if (dbl)fprintf(fl, " \n\n");
    T = (double *) calloc(ny*nx, sizeof (double *));
    if (T == NULL) {
        if (dbo)printf("** Erro: Memoria Insuficiente **");
        if (dbl)fprintf(fl, "** Erro: Memoria Insuficiente **");
        return 0;
    }

    T = gaus_seidel(ny*nx, A, T, B, tol);

    fprintf(fl, "\n\nResultado:\n\n");
    //IMPRINDO O RESULTADO DO VETOR DE TEMPERATURAS
    for (i = nx * ny - 1, c = 0, j = 0; i >= 0 * nx; i--) {
        if (dbo)printf(" %f", T[i]);
        fprintf(fl, " %f", T[i]);
        Tp[c][j] = T[i];
        if (tmax <= T[i])
            tmax = T[i];
        j++;
        if (j % nx == 0) {
            if (dbo)printf("\n");
            fprintf(fl, "\n");
            c++;
            j = 0;
        }
    }

    printf(" \n\nSALVANDO DADOS:\n");
    fprintf(fl, " \n\nSALVANDO DADOS:\n");

    sprintf(filename, "%s/exe5_data.txt", dirname);
    fd = fopen(filename, "w");
    if (fd == NULL) //if file does not exist, create it
    {
        fd = fopen(filename, "wb");
    }
    if (fd == NULL) {
        printf("Erro salvar %s.\n", filename);
        fprintf(fl, "Erro salvar %s.\n", filename);
    } else {
        for (i = 0; i < ny; i++) {
            for (j = 0; j < nx; j++) {
                fprintf(fd, "%f %f %f\n", j*dx, ((ny - 1) - i) * dy, Tp[i][j]);
            }
            fprintf(fd, "\n");
        }
    }

    printf(" \n\tResultado Numérico salvo em exe5_data.txt\n");
    fprintf(fl, " \n\tResultado Numérico salvo em exe5_data.txt\n");
    fclose(fd);

    tmax = (tmax > 0) ? tmax : 0.5;
    sprintf(cmd, "gnuplot -c exe5_plotdata.gnu %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %s %s", a, b, dx, dy, u, v, tmax, filename, dirname);
    system(cmd);

    //SOLUÇÃO EXATA PARA u = v = 0
    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) {
            if (i == ny - 1) {
                Ta[i][j] = sin((PI * (j * dx + dx / 2)) / a);
            } else {
                Ta[i][j] = (sinh((PI * b) / a) / sinh((PI * ((ny - 1) - i * dy + dy / 2) / a))) * sin((PI * (j * dx + dx / 2)) / a);
            }
        }
    }
    if (dbo)printf("\n\nSOLUÇÃO EXATA PARA u = v = 0\n\n");
    if (dbl)fprintf(fl, "\n\nSOLUÇÃO EXATA PARA u = v = 0\n\n");
    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) {
            if (dbo)printf(" %f", Ta[i][j]);
            if (dbl)fprintf(fl, " %f", Ta[i][j]);
        }
        if (dbo)printf("\n");
        if (dbl)fprintf(fl, "\n");
    }

    if (dbl == 2) {
        if (dbo)printf("Erro salvar %s.\n", logfilename);
    } else {
        if (dbl == 1)fclose(fl);
    }

    if (dbo)printf(" \n\tLog de execução salvo em %s\n", logfilename);

    return (EXIT_SUCCESS);
}

double *gaus_seidel(unsigned int n, double **Ap, double x[], double b[], double tol) {
    int i, j, iter = 0;
    double sigma = 0, xa[n];

    printf("\nMÉTODO DE GAUS-SEIDEL:\nMatriz [A b]\n");
    for (i = 0; i < n; i++) {
        xa[i] = x[i] = 0;
    }
    printf("\n%u\n", n);
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
    } while (norma_vetor(xa, x, n) > tol && iter < 100000 * n);

    if (iter < 100000 * n) {
        printf("\n\t%.10f <= %.7f Tolerância\n", norma_vetor(xa, x, n), tol);
        printf("\nNúmero de iterações: %d\n", iter);
        printf("\n\n");
    } else {
        printf("\nNúmero de iterações: %d\n\t%u", iter, 100000 * n);
        printf("\n\n");
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

void make_directory(const char* name) {
#ifdef __linux__
    mkdir(name, 777);
#else
    _mkdir(name);
#endif
}
