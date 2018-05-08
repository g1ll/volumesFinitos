/*
VOLUMES FINITOS
Gill Velleda Gonzales
*/
#include <stdio.h>

int main(){
	
	int t = 6, 				//passo de tempo total
		p = 5, 				//número de pontos
		i,j; 				//indices	

	float r,dt,dx=1,a=1, 	//parâmetros
		  m[t][p]; 			// matriz com resultados

	printf("\nVOLUMES FINITOS - Exercício 3.6 Cap. 3 :\n");

	for(dt=0.25;dt<=0.75;dt+=0.25){		

		printf("\n\nSolução do Problema com dt = %.2f:\n",dt);
		printf("\n T\tp1\tp2\tp3\tp4\tp5\n");
		
		r = a*(dt/(dx*dx));			//calculando parametro r
		
		//printf("\nr = %.3f",r);		

		for(i=0;i<t;i++){			//passo do tempo

			printf("\n %d",i+1);

			for(j=0;j<p;j++){ 		//calculando cada ponto
				if(j == 0 || j == p-1){
					m[i][j] = 1;	//temperatura nas fronteiras direita e esquerda
				}else{
					if(i==0){
						m[i][j] = 0; //temperatura no passo de tempo 0
					}else{						
						m[i][j] = r*(m[i-1][j+1]+m[i-1][j-1]);
					}
				}						

				printf("\t%0.2f",m[i][j]);
			}
		}

		printf("\n\n");
	}

	printf("\n\n\n");
	return 0;
}
