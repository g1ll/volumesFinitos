#include <stdio.h>
#include <math.h>

int main() {
    int     i, j,               //indices
            p = 5,              //numero de pontos da malha (Limite Máx 54 pontos - Resolver Bug) 
            pd = 2 * p - 1,     //Numero total de pontos contando com os volumes
            n = p - 1,          //Número de Volumes
            linha;              //Indice de linha usado na fatoração LU
    double  l = 1,              //Tamanho do domínio computacional 1D em metros
            dx = l / (pd - 1),  //Discretização em x
            T[pd],              //Temperatura para cada ponto (Solução Exata)
            Td[pd],             //Temperatura para cada ponto (Diferença Finita)
            Tad[pd],            //Temperatura iteracao anterior (jacobi) (DF)
            x,                  //Variável de posição x (Solução Exata)
            gk = 5,             //Termo fonte (geração de calor)
            A[n][n],            //Matriz de coeficientes
            b[n],               //Matriz B 
            Tn[n],              //Vetor com temperaturas nodais
            Tk[n],              //Vetor com temperaturas nodais anteriores (jacobi)
            U[n][n],            //Matriz U (Decompozição LU)
            L[n][n],            //Matriz L (Decompozição LU)
            soma_d,             //Soma da linha abaixo da diagonal (LU)
            multi,              //Multiplo (LU)	
            somajc = 0,         //Soma Jacobi
            tol = 1e-5;         //Tolerância para o método Jacobi

    printf("\n\n\u250F");
    for (i = 0; i < 83; i++) {
        printf("\u2505");
    }
    printf("\u2513\n\u2507");
    for (i = 0; i < 22; i++) {
        printf("\u2591");
    }
    printf("VOLUMES FINITOS - Exercício 3.8 Cap. 3:");
    for (i = 0; i < 22; i++) {
        printf("\u2591");
    }
    printf("\u2507\n\u255A");
    for (i = 0; i < 83; i++) {
        printf("\u2505");
    }
    printf("\u255D");


    printf("\n\n\nDados do Problema:");
    printf("\n\tComprimento l=: %.1f m\n\tMalha: %d pontos", l, p);
    printf("\n\t\u0394x=:%.3f\n\tg/k = :%.2f", dx, gk);
    printf("\nTol: %10f\n", tol);

    T[0] = 0; //Condição de contorno esquerda
    T[pd - 1] = 0; //Condição de contorno direita

    printf("\nCondições de Contorno:\n\tT(0) = %.2f\n\tT(%.1f)= %.2f", T[0], l, T[pd - 1]);
    printf("\n\n\n");
    getchar();
    //-------------------------------------------------------------------------------
    //==========================SOLUCAO EXATA========================================
    for (i = 0; i < 80; i++) {
        printf("\u2508");
    }
    printf("\n\t\u2554");
    for (i = 0; i < 17; i++) {
        printf("\u2550");
    }
    printf("\u2557\n\t\u2551\u222B\u2202x SOLUCAO EXATA\u2551\n\t\u255A");
    for (i = 0; i < 17; i++) {
        printf("\u2550");
    }
    printf("\u255D");

    //solução exata
    for (x = 0, i = 0; i < pd; i++, x += dx) {
        T[i] = -gk / 2 * (x * x) + gk * x / 2; // T = -(g/2k)x² + c_1x + c_2
        //T[i] = i+1;
        printf("\n\tdx= %.4f| T=%.4f\n", x, T[i]);
    }
    printf("\n\n\tResultado:\n\nTemperaturas nos pontos \u26AB (K):");
    for (i = 0; i < pd; i += 2) {
        printf("\t %.4f", T[i]);
    }
    printf("\nTemperaturas nos nós \u26AA (K):");
    for (i = 1; i < pd; i += 2) {
        printf("\t %.4f", T[i]);
    }
    printf("\n\n\nPressione qualquer tecla para continuar!\n\n");
    getchar();
    //--------------------------------------------------------------------------------
    //=====================SOLUÇÃO POR DIFERENCAS FINITAS=============================

    for (i = 0; i < 80; i++) {
        printf("\u2508");
    }
    printf("\n\t\u2554");
    for (i = 0; i < 50; i++) {
        printf("\u2550");
    }
    printf("\u2557\n\t\u2551\u2206\u2245 SOLUÇÃO POR DIFERENCAS FINITAS e Método JACOBI:\u2551\n\t\u255A");
    for (i = 0; i < 50; i++) {
        printf("\u2550");
    }
    printf("\u255D");

    for (i = 0; i < pd; i++) {
        Td[i] = Tad[i] = 0;
    } //Limpando os vetores
    printf("\n\n\tMalha:\t \u0394x %.3f", dx);

    //Método de Jacobi
    x = 0; //contador de iteracoes
    do {
        //Td[0]=0;			//Condição de contorno esquerda
        //Tad[pd-1]=0;		//Condição de contorno direita
        somajc = 0;
        x++;
        for (i = 1; i < pd - 1; i++) {
            Tad[i] = Td[i];
        } //Atualizando T anterior
        for (i = 1; i < pd - 1; i++) {
            Td[i] = (Tad[i - 1] + Tad[i + 1] + gk * (dx * dx)) / 2; //Aproximacao Dif. Finita 2o.
        }
        for (i = 1; i < pd; i++) {
            somajc += Td[i] - Tad[i];
        }
    } while (sqrt(pow(somajc, 2)) > tol && x < 2000); //Critério de parada de jacobi
    if (x >= 2000) {
        printf("\nNÃO CONVERGIU!!!\n");
    }
    printf("\n\tJacobi %.0f iter. diferenca de Td - Tad = %f", x, somajc);
    printf("\n\n\tResultado:\n\nTemperaturas nos pontos \u26AB (K): ");
    for (i = 0; i < pd; i += 2) {
        printf("\t %.4f", Td[i]);
    }
    printf("\nTemperaturas nos nós \u26AA  (K): ");
    for (i = 1; i < pd; i += 2) {
        printf("\t %.4f", Td[i]);
    }
    printf("\n\n\nPressione qualquer tecla para continuar!\n\n");
    getchar();
    //--------------------------------------------------------------------------	
    //=====================SOLUCAO POR VOLUMES FINITOS===========================
    for (i = 0; i < 80; i++) {
        printf("\u2508");
    }
    printf("\n\t\u2554");
    for (i = 0; i < 30; i++) {
        printf("\u2550");
    }
    printf("\u2557\n\t\u2551\u25A6 SOLUCAO POR VOLUMES FINITOS:\u2551\n\t\u255A");
    for (i = 0; i < 30; i++) {
        printf("\u2550");
    }
    printf("\u255D");

    dx = l / n; //Atualizando o dx para Volumes Finitos

    printf("\nSistema Linear Ax = b:\n\n");
    for (i = 0; i < n; i++) {
        b[i] = gk*dx; //Vetor b do sistema Ax = b
        for (j = 0; j < n; j++) { //Montanto a Matriz de coeficientes A
            if (j == i + 1 || j == i - 1) {
                A[i][j] = -1 / dx;
            } else if (j == i) {
                if (i != 0 && i != n - 1) {
                    A[i][j] = 2 / dx;
                } else {
                    A[i][j] = 3 / dx;
                }
            } else {
                A[i][j] = 0;
            }
            printf("\t %-2.2f", A[i][j]);
        }
        printf("\t%-2.2f\n", b[i]);
    }

    //Resolvendo o sistema linear por Jacobi [--ERROS--]
    /*
    printf("\n\nResolvendo o sistema linear por Jacobi :\n");
    for(i=0;i<n;i++){Tn[i]=Tk[i]=1;} //Limpando os vetores
    x=0;//contador de iteracoes
    do{
            x++; 
            for(i=0;i<n;i++){Tk[i]=Tn[i];}
            for(i=0;i<n;i++){
                    for(j=0;j<n;j++){
                            if(j!=i){
                                    Tn[i]+= Tk[j]*A[i][j];
                            }
                    }
                    Tn[i] = (Tn[i]+b[i])/A[i][i];
            }	
    }while(fabs(Tn[n/2]-Tk[n/2]) > 0.0001); //Critério de parada de jacobi
    printf("\n\tJacobi %.0f iter. diferenca de Tn - Tk = %f - %f = %f",x,Tn[n/2],Tk[n/2],fabs(Tn[n/2]-Tk[n/2]));	
    printf("\n\tTemperatura Nodal (K):");
    for(i=0;i<n;i++){ printf("\t %.4f",Tn[i]);}
     */

    //Resolvendo o sistema linear por decomposição LU
    //Calculando as Matrizes U e L
    printf("\n\nResolvendo o sistema por Decomposição LU:\n");

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            U[i][j] = A[i][j];
            L[i][j] = 0;
        }
    }
    printf("\n");
    linha = 0;
    L[0][0] = 1;
    do {
        for (i = linha + 1; i < n; i++) {
            multi = U[i][linha] / U[linha][linha];
            //printf("\n\tmulti: %f",multi);//imprimindo o multiplicador
            for (j = 0; j < n; j++) {
                U[i][j] = U[i][j] - multi * U[linha][j];
            }
            L[i][i] = 1;
            L[i][linha] = multi;
        }

        /*printf("\n\n\n");//imprimindo passo a passo	
        for(i=0;i<n;i++){
                printf("\n");
                for(j=0;j<n;j++){
                        printf("\t %.3f",U[i][j]);
                }
        }*/
        soma_d = 0;
        for (i = 1; i < n; i++) {
            soma_d += U[i][i - 1];
        }
        //printf("\n\n\tsoma: %.3f",soma_d);
        linha++;
    } while (sqrt(pow(soma_d, 2)) > tol);

    printf("\nMatriz U:\n");
    for (i = 0; i < n; i++) {
        printf("\n");
        for (j = 0; j < n; j++) {
            printf("\t %.3f", U[i][j]);
        }
    }

    //Matriz L
    printf("\n\nMatriz L:\n");
    for (i = 0; i < n; i++) {
        printf("\n");
        for (j = 0; j < n; j++) {
            printf("\t %.3f", L[i][j]);
        }
    }
    //Resolvendo o sistema Ly=b
    printf("\n\n(Ly=b) y:\n");
    for (i = 0; i < n; i++) {
        Tk[i] = 0;
    } //Zerando Tk(y)
    Tk[0] = b[0];
    printf("\n\t %.4f", Tk[0]);
    for (i = 1; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i != j) {
                Tk[i] += L[i][j] * Tk[j];
            }
            //printf("\n %d %f %.3f %f",j,L[i][j],Tk[i],Tk[j]);	
        }
        //printf("\n\n\n %.3f\n\n\n",Tk[i]);
        Tk[i] = b[i] - Tk[i];
        printf("\n\t %.4f", Tk[i]);
    }
    //Resolvendo o sistema Ux=y
    printf("\n\n(Ux=y) x:\n");
    for (i = 0; i < n; i++) {
        Tn[i] = 0;
    } //Zerando Tn(x) 
    for (i = n - 1; i >= 0; i--) {
        for (j = 0; j < n; j++) {
            if (i != j) {
                Tn[i] += U[i][j] * Tn[j];
            }
        }
        //printf("\n %d %f %f %f",i,Tn[i],Tk[i], U[i][i]);
        Tn[i] = (Tk[i] - Tn[i]) / U[i][i];
    }
    for (i = 0; i < n; i++) {
        printf("\n\t%d %.4f", i, Tn[i]);
    };
    printf("\n\n\tResultado:\n\nTemperatura nos nós \u26AA  (K):");
    for (i = 0; i < n; i++) {
        printf("\t %.4f", Tn[i]);
    }
    printf("\nTemperaturas  nos pontos \u26AB (K):");
    for (i = 0; i < p; i++) {
        T[i] = 0;
    }
    for (i = 0; i < n - 1; i++) {
        T[i + 1] = (Tn[i] + Tn[i + 1]) / 2;
    }
    for (i = 0; i < p; i++) {
        printf("\t %.4f", T[i]);
    }
    printf("\n\t\u21D2 calculado por diferenças centrais:\n");
    printf("\n\n\n");
}
