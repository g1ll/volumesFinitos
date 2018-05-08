/************************************************************************
Main program to call subroutines to plot contours of a function specified
by numerical data on a rectangular grid
contr_lines - generates contour lines and labels with contour heights
contr_shade - generates contour lines with colored shading between contours
Each generates a postscript file.  T.W. Secomb, February 2009.
This code is free to use at your own risk.  Feedback and/or acknowledgement
are appreciated.  Email secomb@u.arizona.edu.
Note: uses nrutil.h and nrutil.cpp as placed in the public domain
by Numerical Recipes at http://www.nrbook.com/a/bookcpdf/c21-1.pdf.
***********************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void contr_lines(FILE *ofp, int m, int n, float scalefac, int nl,float xmin,
		   float xmax, float ymin, float ymax, float *cl, float **zv);
void contr_shade(FILE *ofp, int m, int n, float scalefac, int nl,float xmin,
		   float xmax, float ymin, float ymax, float *cl, float **zv);


int main(int argc, char *argv[])
{
	int i,j,m,n,nl;
	float clmin,clint,xmin,xmax,ymin,ymax,scalefac,x,y;
	float **zv,*cl;
	FILE *ofp;
	nl = 17;
	m = 20;
	n = 20;

	cl = vector(1,nl);
	zv = matrix(1,m,1,n);

	for(i=1; i<=m; i++) for(j=1; j<=n; j++) {
		x = 2.*float(i-1)/float(m-1) - 1.;
		y = 2.*float(j-1)/float(n-1) - 1.;
		zv[i][j] = (cos(4.*x)/(1. + x*x) + sin(4.*y)/(1 + y*y))/2.;
//		zv[i][j] = cos(4.*x)/(1. + x*x)*sin(4.*y)/(1 + y*y);
	}
	xmin = -1.;
	ymin = -1.;
	xmax = 1.;
	ymax = 1.;
	clmin = -0.7;
	clint = 0.1;
	for(i=1; i<=nl; i++) cl[i] = clmin + (i-1)*clint;
	scalefac = 300./__max(xmax-xmin,ymax-ymin);

	ofp = fopen("Contour_lines.ps", "w");
	contr_lines(ofp, m, n, scalefac, nl, xmin, xmax, ymin, ymax, cl, zv);
//These commands are in the main program to allow further additions to page.
	fprintf(ofp, "showpage\n");
	fclose(ofp);

	ofp = fopen("Contour_shade.ps", "w");
	contr_shade(ofp, m, n, scalefac, nl, xmin, xmax, ymin, ymax, cl, zv);
//These commands are in the main program to allow further additions to page.
	fprintf(ofp, "showpage\n");
	fclose(ofp);
	
	return 0;
}