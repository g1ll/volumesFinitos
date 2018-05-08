

#include <stdio.h>
#include <stdlib.h>
#include<math.h>

int main(int argc, char** argv) {

    int i;
    float n, p;

    printf("\n\tPotência:\n\n");
    
    if (argc > 1) {
        for (i = 1; i < argc; i += 2) {
            n = atof(argv[i]);
            if (i + 1 == argc) {
                p = n;
            } else {
                p = atof(argv[i + 1]);
            }
            printf("\n\t %0.f^%0.f = %0.f\n\n", n, p, pow(n, p));
        }
    } else {
        n = p = 2;
        printf("\n\t %0.f^%0.f = %0.f\n\n", n, p, pow(n, p));
    }
    
    //printf("\n\t %d\n\n", argc);
    return (EXIT_SUCCESS);
}

