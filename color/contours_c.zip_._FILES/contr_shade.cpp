/*******************************************************************
contr_shade - generates contour lines with colored shading between contours
Variables - 
m n nl:  dimensions of array, no. of contour levels
scalefac: determines size of plot
xmin xmax ymin ymax:  boundaries of box
cl(nl):  array of contour levels
zv(m,n):  array of heights
Output to a postscript file.
Version including color shading of region, TWS February 2009.
Note:  This could be improved by combining polygons with the same color
before writing to the PostScript file.
*********************************************************************/

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void contr_shade(FILE *ofp, int m, int n, float scalefac, int nl,float xmin,
		   float xmax, float ymin, float ymax, float *cl, float **zv)
{
	fprintf(ofp, "%%!PS\n");
	int i,j,k,iwsp,nsq,nsq1,isq,isq1,kmin,kmax;
	int in,in1,in2,in3,in4;
	int *ii,*jj;
	int **corners;
	const int ncmax=10000;
	float dj1,di1,c,zmin,zmax,xz,d,xv12,yv12,xv23,yv23,xv34,yv34,xv41,yv41;
	float cx,cy,cbx,cby,cbbox;
	float *xv,*yv,*red,*blue,*green,*dd;
	float **wsp;

	cx = 50; //origin of contour plot
	cy = 100;//origin of contour plot
	cbbox = 15.; //size of boxes
	cbx = 370; //origin of color bar
	cby = 100;//origin of color bar

	xv = vector(1,m);
	yv = vector(1,n);
	red = vector(0,nl);
	green = vector(0,nl);
	blue = vector(0,nl);
	wsp = matrix(1,ncmax,1,4);
	corners = imatrix(1,4,1,2);
	dd = vector(1,4);
	ii = ivector(1,4);
	jj = ivector(1,4);
	corners[1][1] = 0;
	corners[1][2] = 0;
	corners[2][1] = 1;
	corners[2][2] = 0;
	corners[3][1] = 1;
	corners[3][2] = 1;
	corners[4][1] = 0;
	corners[4][2] = 1;

	for(i=1; i<=m; i++) xv[i] = xmin + (i - 1)*(xmax - xmin)/(m - 1);
	for(j=1; j<=n; j++) yv[j] = ymin + (j - 1)*(ymax - ymin)/(n - 1);
	fprintf(ofp, "/mx {%g sub %g mul %g add} def\n",xmin,scalefac,cx);
	fprintf(ofp, "/my {%g sub %g mul %g add} def\n",ymin,scalefac,cy);
	fprintf(ofp, "/m {moveto} def\n");
	fprintf(ofp, "/l {lineto} def\n");
	fprintf(ofp, "/c {closepath setrgbcolor fill} def\n");
	fprintf(ofp, "0.5 setlinewidth\n");
	fprintf(ofp, "/Times-Roman findfont\n");
	fprintf(ofp, "8 scalefont\n");
	fprintf(ofp, "setfont\n");

//Set up colors using Matlab 'jet' scheme
	for(k=0; k<=nl; k++){
		xz = float(k)/float(nl);
		blue[k] = __min(__max(1.5-4*abs(xz-0.25), 0.), 1.);
		green[k]= __min(__max(1.5-4*abs(xz-0.5), 0.), 1.);
		red[k]  = __min(__max(1.5-4*abs(xz-0.75), 0.), 1.);
	}
//Analyze each rectangle separately.  Find max and min z levels.
	iwsp = 0;
	for(i=1; i<m; i++) for(j=1; j<n; j++){
		for(in=1; in<=4; in++){
			ii[in] = i + corners[in][1];
			jj[in] = j + corners[in][2];
			if(in == 1){
				zmin = zv[ii[in]][jj[in]];
				zmax = zv[ii[in]][jj[in]];
			}
			else{
				zmin = __min(zv[ii[in]][jj[in]],zmin);
				zmax = __max(zv[ii[in]][jj[in]],zmax);
			}
		}
		kmin = 0;	//below the lowest contour
		kmax = nl;	//above the highest contour
		for(k=1; k<=nl; k++) if(zmin > cl[k]) kmin = k;
		for(k=nl; k>0; k--)
			if(zmax < cl[k]) kmax = k-1;
//Color each rectangle.  If more colors are needed, overwrite the lower colors
		for(k=kmin; k<=kmax; k++){
			in1 = 0;
			if(k == 0) in1 = 1;
			else for(in=4; in>=1; in--) if(zv[ii[in]][jj[in]] > cl[k]) in1 = in;	//find a corner this color or higher
			if(in1 == 0) printf("*** Error in contrsub: search failed\n");
			in2 = in1%4 + 1;
			in3 = in2%4 + 1;
			in4 = in3%4 + 1;
			for(in=1; in<=4; in++) if(k == 0) dd[in] = 1.; else dd[in] = zv[ii[in]][jj[in]] - cl[k];
			if(dd[in1] != dd[in2]){
				xv12 = (dd[in1]*xv[ii[in2]]-dd[in2]*xv[ii[in1]])/(dd[in1]-dd[in2]);
				yv12 = (dd[in1]*yv[jj[in2]]-dd[in2]*yv[jj[in1]])/(dd[in1]-dd[in2]);
			}
			if(dd[in2] != dd[in3]){
				xv23 = (dd[in2]*xv[ii[in3]]-dd[in3]*xv[ii[in2]])/(dd[in2]-dd[in3]);
				yv23 = (dd[in2]*yv[jj[in3]]-dd[in3]*yv[jj[in2]])/(dd[in2]-dd[in3]);
			}
			if(dd[in3] != dd[in4]){
				xv34 = (dd[in3]*xv[ii[in4]]-dd[in4]*xv[ii[in3]])/(dd[in3]-dd[in4]);
				yv34 = (dd[in3]*yv[jj[in4]]-dd[in4]*yv[jj[in3]])/(dd[in3]-dd[in4]);
			}
			if(dd[in4] != dd[in1]){
				xv41 = (dd[in4]*xv[ii[in1]]-dd[in1]*xv[ii[in4]])/(dd[in4]-dd[in1]);
				yv41 = (dd[in4]*yv[jj[in1]]-dd[in1]*yv[jj[in4]])/(dd[in4]-dd[in1]);
			}
			fprintf(ofp, "newpath %g mx %g my m\n",xv[ii[in1]],yv[jj[in1]]);
			if(dd[in2] > 0){													//corners 1,2 are this color
				fprintf(ofp, "%g mx %g my l\n",xv[ii[in2]],yv[jj[in2]]);
				if(dd[in3] > 0){												//corners 1,2,3 are this color
					fprintf(ofp, "%g mx %g my l\n",xv[ii[in3]],yv[jj[in3]]);
					if(dd[in4] > 0)
						fprintf(ofp, "%g mx %g my l\n",xv[ii[in4]],yv[jj[in4]]);//corners 1,2,3,4 are this color
					else{														//corners 1,2,3,not 4 are this color
						fprintf(ofp, "%g mx %g my l\n",xv34,yv34);
						fprintf(ofp, "%g mx %g my l\n",xv41,yv41);
						iwsp += 1;
						wsp[iwsp][1] = xv34;
						wsp[iwsp][2] = yv34;
						wsp[iwsp][3] = xv41;
						wsp[iwsp][4] = yv41;
					}
				}
				else{															//corners 1,2,not 3 are this color
					fprintf(ofp, "%g mx %g my l\n",xv23,yv23);
					iwsp += 1;
					wsp[iwsp][1] = xv23;
					wsp[iwsp][2] = yv23;
					if(dd[in4] > 0){											//corners 1,2,not 3,4 are this color
						fprintf(ofp, "%g mx %g my l\n",xv34,yv34);
						wsp[iwsp][3] = xv34;
						wsp[iwsp][4] = yv34;
						fprintf(ofp, "%g mx %g my l\n",xv[ii[in4]],yv[jj[in4]]);
					}
					else{
						fprintf(ofp, "%g mx %g my l\n",xv41,yv41);				//corners 1,2,not 3,not 4 are this color
						wsp[iwsp][3] = xv41;
						wsp[iwsp][4] = yv41;
					}
				}
			}
			else{																//corners 1,not 2 are this color
				fprintf(ofp, "%g mx %g my l\n",xv12,yv12);
				iwsp += 1;
				wsp[iwsp][1] = xv12;
				wsp[iwsp][2] = yv12;
				if(dd[in3] > 0){												//corners 1,not 2,3 are this color
					fprintf(ofp, "%g mx %g my l\n",xv23,yv23);
					wsp[iwsp][3] = xv23;
					wsp[iwsp][4] = yv23;
					fprintf(ofp, "%g mx %g my l\n",xv[ii[in3]],yv[jj[in3]]);
					if(dd[in4] > 0)
						fprintf(ofp, "%g mx %g my l\n",xv[ii[in4]],yv[jj[in4]]);//corners 1,not 2,3,4 are this color
					else{														//corners 1,not 2,3,not 4 are this color
						fprintf(ofp, "%g mx %g my l\n",xv34,yv34);
						fprintf(ofp, "%g mx %g my l\n",xv41,yv41);
						iwsp += 1;
						wsp[iwsp][1] = xv34;
						wsp[iwsp][2] = yv34;
						wsp[iwsp][3] = xv41;
						wsp[iwsp][4] = yv41;
					}
				}
				else{															//corners 1,not 2,not 3 are this color
					if(dd[in4] > 0){											//corners 1,not 2,not 3,4 are this color
						fprintf(ofp, "%g mx %g my l\n",xv34,yv34);
						wsp[iwsp][3] = xv34;
						wsp[iwsp][4] = yv34;
						fprintf(ofp, "%g mx %g my l\n",xv[ii[in4]],yv[jj[in4]]);
					}
					else{
						fprintf(ofp, "%g mx %g my l\n",xv41,yv41);				//corners 1,not 2,not 3,not 4 are this color
						wsp[iwsp][3] = xv41;
						wsp[iwsp][4] = yv41;
					}
				}
			}
			if(iwsp > ncmax-4) printf("*** Error: ncmax too small in contr\n");						
			fprintf(ofp, "%6.3f %6.3f %6.3f c\n",red[k],green[k],blue[k]);
		}
	}
//Now outline contours
	fprintf(ofp, "0 0 0 setrgbcolor\n");//black
	for(in=1; in<=iwsp; in++) fprintf(ofp, "newpath %g mx %g my m %g mx %g my l stroke\n",
			wsp[in][1],wsp[in][2],wsp[in][3],wsp[in][4]);
//Draw a box
	fprintf(ofp, "0 0 0 setrgbcolor\n");//black
	fprintf(ofp, "newpath\n");
	fprintf(ofp, "%g mx %g my m %g mx %g my l %g mx %g my l %g mx %g my l\n",
		xmin,ymin,xmax,ymin,xmax,ymax,xmin,ymax);
	fprintf(ofp, "closepath stroke\n");
//Create a color bar
	for(k=0; k<=nl; k++){
		fprintf(ofp, "newpath %g %g m %g %g l %g %g l %g %g l closepath\n",
			cbx,cby+k*cbbox,cbx+cbbox,cby+k*cbbox,cbx+cbbox,cby+(k+1)*cbbox,cbx,cby+(k+1)*cbbox);
		fprintf(ofp, "%f %f %f setrgbcolor fill\n",red[k],green[k],blue[k]);
		if(k>0) fprintf(ofp, "%g %f m 0 0 0 setrgbcolor (%g) show\n",cbx+cbbox*1.1,cby+cbbox*(k-0.1),cl[k]);
	}
	fprintf(ofp, "newpath %g %g m %g %g l %g %g l %g %g l closepath stroke\n",
		cbx,cby,cbx+cbbox,cby,cbx+cbbox,cby+cbbox*(nl+1),cbx,cby+cbbox*(nl+1));
}