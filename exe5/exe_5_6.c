/**
* 
* REFERENCES:
* https://en.wikipedia.org/wiki/C_date_and_time_functions
*
**/

#define M_PI 3.14159265358979323846
#define __FILENAME__ (strrchr(__FILE__, '\\') ? strrchr(__FILE__, '\\') + 1 : __FILE__)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>


double *gaus_seidel2(int ny, int nx, double **Ap, double **Ae, double **Aw, double **An, double **As, double x[], double b[], double tol);

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
            abe, abw, abn, abs, //Variáveis para guardar a expresão de interpolação
            **Ap, //Matriz de coeficientes Ap
            **Ae, //Matriz de coeficientes Ae
            **Aw, //Matriz de coeficientes Aw
            **An, //Matriz de coeficientes An
            **As, //Matriz de coeficientes As
            B[nx * ny], //Vetor b do sistema Ax=b (fontes)
            *T, //Vetor com Temperaturas
            Tp[ny][nx], //Campo de temperaturas solução numérica
            Ta[ny][nx], //Campo de temperaturas solução analítica
            tmax = -9999, //Temperatura máxima
            tol = 1e-10; //Tolerância para os métodos Iterativos

    /*DATA ATUAL*/
    current_time = time(NULL);
    printf("TESTE.\n");
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
    printf("TESTE2.\n %s",dirname);
   
    #ifdef __linux__
        sprintf(logfilename, "%s/exe5_log.txt", dirname);
    #else
        //sprintf(logfilename, "D:\\Projetos\\volumes_finitos\\exe5\\%s\\exe5_log.txt", dirname);
        sprintf(logfilename, "%s\\exe5_log.txt", dirname);
    #endif

    printf("\n%s\n",logfilename);
    fl = fopen(logfilename, "w");
    if (fl == NULL) //if file does not exist, create it
    {
        fl = fopen(logfilename, "wb");
    }
    
    if (fl == NULL) {
        //return 0;
        dbl = 0;
        printf("NÃO CRIA O ARQUIVO!.\n");
    }

    //CALCULANDO PECLET, ALFA E BETA PARA X E Y (WUDS)
    pex = (u * ro * dx) / ga; //Peclet em x
    pey = (v * ro * dy) / ga; //Peclet em y
    ax = pow(pex, 2) / (10 + 2 * pow(pex, 2)); //Calculando alfa em x
    ay = pow(pey, 2) / (10 + 2 * pow(pey, 2)); //Calculando alfa em y
    bx = (1 + 0.005 * pow(pex, 2)) / 1 + 0.05 * pow(pex, 2); //Cálculo de beta em x       
    by = (1 + 0.005 * pow(pey, 2)) / 1 + 0.05 * pow(pey, 2); //Cálculo de beta em y

    // Expressões de interpolação WUDS
    abe = -(1 / 2 - ax) * ro * u + bx * ga / dx;
    abw = ro * u * (1 / 2 - ax) + bx * ga / dx;
    abn = -(1 / 2 - ay) * ro * v + by * ga / dy;
    abs = ro * v * (1 / 2 + ay) + by * ga / dy;

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

    /* aloca as linhas da matriz Ap */
    Ap = (double **) calloc(ny, sizeof (double *));
    if (Ap == NULL) {
        if (dbo)printf("** Erro: Memoria Insuficiente **");
        if (dbl)fprintf(fl, "** Erro: Memoria Insuficiente **");
        return 0;
    }

    /* aloca as colunas da matriz Ap */
    for (i = 0; i < nx; i++) {
        Ap[i] = (double*) calloc(nx, sizeof (double));
        if (Ap[i] == NULL) {
            if (dbo)printf("** Erro: Memoria Insuficiente **");
            if (dbl)fprintf(fl, "** Erro: Memoria Insuficiente **");
            return 0;
        }
    }

    /* aloca as linhas da matriz Ae */
    Ae = (double **) calloc(ny, sizeof (double *));
    if (Ae == NULL) {
        if (dbo)printf("** Erro: Memoria Insuficiente **");
        if (dbl)fprintf(fl, "** Erro: Memoria Insuficiente **");
        return 0;
    }

    /* aloca as colunas da matriz Ae */
    for (i = 0; i < nx; i++) {
        Ae[i] = (double*) calloc(nx, sizeof (double));
        if (Ae[i] == NULL) {
            if (dbo)printf("** Erro: Memoria Insuficiente **");
            if (dbl)fprintf(fl, "** Erro: Memoria Insuficiente **");
            return 0;
        }
    }

    /* aloca as linhas da matriz Aw */
    Aw = (double **) calloc(ny, sizeof (double *));
    if (Aw == NULL) {
        if (dbo)printf("** Erro: Memoria Insuficiente **");
        if (dbl)fprintf(fl, "** Erro: Memoria Insuficiente **");
        return 0;
    }

    /* aloca as colunas da matriz Aw */
    for (i = 0; i < nx; i++) {
        Aw[i] = (double*) calloc(nx, sizeof (double));
        if (Aw[i] == NULL) {
            if (dbo)printf("** Erro: Memoria Insuficiente **");
            if (dbl)fprintf(fl, "** Erro: Memoria Insuficiente **");
            return 0;
        }
    }

    /* aloca as linhas da matriz An */
    An = (double **) calloc(ny, sizeof (double *));
    if (An == NULL) {
        if (dbo)printf("** Erro: Memoria Insuficiente **");
        if (dbl)fprintf(fl, "** Erro: Memoria Insuficiente **");
        return 0;
    }

    /* aloca as colunas da matriz An */
    for (i = 0; i < nx; i++) {
        An[i] = (double*) calloc(nx, sizeof (double));
        if (An[i] == NULL) {
            if (dbo)printf("** Erro: Memoria Insuficiente **");
            if (dbl)fprintf(fl, "** Erro: Memoria Insuficiente **");
            return 0;
        }
    }

    /* aloca as linhas da matriz As */
    As = (double **) calloc(ny, sizeof (double *));
    if (As == NULL) {
        if (dbo)printf("** Erro: Memoria Insuficiente **");
        if (dbl)fprintf(fl, "** Erro: Memoria Insuficiente **");
        return 0;
    }

    /* aloca as colunas da matriz As */
    for (i = 0; i < nx; i++) {
        As[i] = (double*) calloc(nx, sizeof (double));
        if (As[i] == NULL) {
            if (dbo)printf("** Erro: Memoria Insuficiente **");
            if (dbl)fprintf(fl, "** Erro: Memoria Insuficiente **");
            return 0;
        }
    }


    //TERMOS FONTE
    for (i = 0; i < ny * nx; i++) {
        B[i] = q;
    }
    //CANTO SUPERIOR ESQUERDO C1
    i = 0;
    j = 0;
    Ae[i][j] = abe;
    Aw[i][j] = 0;
    An[i][j] = 0;
    As[i][j] = abs;
    Ap[i][j] = (2 * bx * ga) / dx + (2 * by * ga) / dy;

    //CANTO SUPERIOR DIREITO C2
    i = 0;
    j = nx - 1;
    Ae[i][j] = 0;
    Aw[i][j] = abw;
    An[i][j] = 0;
    As[i][j] = abs;
    Ap[i][j] = 2 * bx * ga / dx + 2 * by * ga / dy;

    //CANTO INFERIOR ESQUERDO C3
    i = ny - 1;
    j = 0;
    Ae[i][j] = abe;
    Aw[i][j] = 0;
    An[i][j] = abn;
    As[i][j] = 0;
    Ap[i][j] = 2 * bx * ga / dx + 2 * by * ga / dy;
    B[0] = q + ro * v * sin(M_PI * (j * dx + dx / 2) / a) + 2 / dy * sin(M_PI * (j * dx + dx / 2) / a);

    //CANTO INFERIOR DIREITO
    i = ny - 1;
    j = nx - 1;
    Ae[i][j] = 0;
    Aw[i][j] = abw;
    An[i][j] = abn;
    As[i][j] = 0;
    Ap[i][j] = 2 * bx * ga / dx + 2 * by * ga / dy;
    B[j] = q + ro * v * sin(M_PI * (j * dx + dx / 2) / a) + 2 / dy * sin(M_PI * (j * dx + dx / 2) / a);


    for (i = 1; i < nx - 1; i++) {
        B[i] = q + ro * v * sin(M_PI * (i * dx + dx / 2) / a) + 2 / dy * sin(M_PI * (i * dx + dx / 2) / a);
    }

    //CALCULANDO COEFICIENTES Ae, Aw, An, As e Ap PARA CADA PONTO DA MALHA
    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) {
            if (j >= 1 && j <= nx - 2 && i >= 1 && i <= ny - 2) {//VOLUMES DO CENTRO
                Ae[i][j] = abe;
                Aw[i][j] = abw;
                An[i][j] = abn;
                As[i][j] = abs;
                Ap[i][j] = 0;
            } else if (j == 0 && i > 0 && i < ny - 1) { //LATERAL ESQUERDA

                Ae[i][j] = abe;
                Aw[i][j] = 0;
                An[i][j] = abn;
                As[i][j] = abs;
                Ap[i][j] = 2 * bx * ga / dx;
            } else if (j == nx - 1 && i > 0 && i < ny - 1) { //LATERAL DIREITA

                Ae[i][j] = 0;
                Aw[i][j] = abw;
                An[i][j] = abn;
                As[i][j] = abs;
                Ap[i][j] = 2 * bx * ga / dx;
            } else if (i == 0 && j > 0 && j < nx - 1) { //VOLUMES NO TOPO

                Ae[i][j] = abe;
                Aw[i][j] = abw;
                An[i][j] = 0;
                As[i][j] = abs;
                Ap[i][j] = 2 * by * ga / dy;
            } else if (i == ny - 1 && j > 0 && j < nx - 1) { //VOLUMES DA PARTE INFERIOR
                Ae[i][j] = abe;
                Aw[i][j] = abw;
                An[i][j] = abn;
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

    //RESOLVENDO O SISTEMA LINEAR
    if (dbo)printf(" \n\n");
    if (dbl)fprintf(fl, " \n\n");
    T = (double *) calloc(ny*nx, sizeof (double *));
    if (T == NULL) {
        if (dbo)printf("** Erro: Memoria Insuficiente **");
        if (dbl)fprintf(fl, "** Erro: Memoria Insuficiente **");
        return 0;
    }

    T = gaus_seidel2(ny, nx, Ap, Ae, Aw, An, As, T, B, tol);

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
   
    #ifdef __linux__
        sprintf(cmd, "gnuplot -c exe5_plotdata.gnu %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %s %s", a, b, dx, dy, u, v, tmax, filename, dirname);
    #else
        sprintf(cmd, "gnuplot -c exe5_plotdata_win.gnu %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %s %s", a, b, dx, dy, u, v, tmax, filename, dirname);
    #endif   
    system(cmd);

    //SOLUÇÃO EXATA PARA u = v = 0
    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) {
            if (i == ny - 1) {
                Ta[i][j] = sin((M_PI * (j * dx + dx / 2)) / a);
            } else {
                Ta[i][j] = (sinh((M_PI * b) / a) / sinh((M_PI * ((ny - 1) - i * dy + dy / 2) / a))) * sin((M_PI * (j * dx + dx / 2)) / a);
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

double *gaus_seidel2(int ny, int nx, double **Ap, double **Ae, double **Aw, double **An, double **As, double x[], double b[], double tol) {
    int i, j, c, iter = 0;
    unsigned int n = nx*ny;
    double sigma = 0, xa[n];

    printf("\nMÉTODO DE GAUS-SEIDEL:\nMatriz [A b]\n");
    for (i = 0; i < n; i++) {
        xa[i] = x[i] = 0;
    }
    printf("\n%d\n", ny);
    printf("\n%d\n", nx);
    printf("\n%u\n", n);
    printf("\n%f\n", Ap[0][0]);
    printf("\n%f\n", Ae[0][0]);
    printf("iterando..");
    do {
        for (i = n - 1; i >= 0; i--) {
            xa[i] = x[i];
        }
        for (i = 0, c = n - 1; i < ny; i++) {
            for (j = 0; j < nx; j++, c--) {
                sigma = 0;
                if (c - nx >= 0) {
                    sigma += -As[i][j] * x[c - nx];
                }
                if (c - 1 >= 0) {
                    sigma += -Aw[i][j] * x[c - 1];
                }
                if (c + 1 < n) {
                    sigma += -Ae[i][j] * xa[c + 1];
                }
                if (c + nx < n) {
                    sigma += -An[i][j] * xa[c + nx];
                }
                x[c] = (b[c] - sigma) / Ap[i][j];
            }
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
    //#ifdef __linux__
    mkdir(name, 777);
    //#else
    //    _mkdir(name);
    //#endif
}
