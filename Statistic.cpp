




#include <stdio.h>
#include <stdlib.h> 
#include <mpi.h>
#include <omp.h>
#include <math.h>


#include "Resolution.h"
#include "Pre_selection.h"

extern int X_np;

void Statistic
(
// ============================================================================ //
int myid,
int step,
int iteration,
int statistic_step,

double deltaT,

double (*U1_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*U1q)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2q)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3q)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4q)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5q)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*J_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*xidx_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidy_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidz_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*etdx_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdy_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdz_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ztdx_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdy_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdz_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*MR1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ML1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*EpY)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
// ============================================================================ //
)

{

	
#include "ijk.h"
#include "Viscous_terms.h"

#include "MPI_prm.h"
#include "Mpi_division.h"

	char LESdata[100];
	
	double rho,U,V,W,VV,P,C,T,h,H;
	double r,u,v,w;
	double Ek0,Ek,Ekq,ek,t_non;

	Ek = 0.0;


// ============================================================================================================= //

if ((step%statistic_step == 0)) {

	//// ============================================ ////
			 istart = 3;		     	              ////	
	//// ============================================ ////
			iend = gend[myid];				    	  ////
	//// ============================================ ////
		
			for (i = istart; i <= iend; i++) {
				for (j = 2; j <= nz; j++) {
					for (k = 2; k <= nz; k++) {

						rho = U1_[i][j][k];
						U = U2_[i][j][k]/rho;
						V = U3_[i][j][k]/rho;
						W = U4_[i][j][k]/rho;     
						
						ek = 0.5*rho*(U*U+V*V+W*W);

						Ek = Ek + 1.0/1.1842*ek;

						r = U1q[i][j][k];
						u = U2q[i][j][k]/r;
						v = U3q[i][j][k]/r;
						w = U4q[i][j][k]/r; 

						
						ek = 0.5*r*(u*u+v*v+w*w);

						Ekq = Ekq + 1.0/1.1842*ek;

					}
				}
			}

			Ek = -(Ek-Ekq);


	MPI_Comm comm;
	comm=MPI_COMM_WORLD;
	icount = 1;
	idest = 0;
	MPI_Reduce ((void*)&Ek, (void*)&Ek0, icount, MPI_DOUBLE, MPI_SUM, idest, comm);	

	if (myid == 0) {
	// =============================================================================================================== //

		FILE *fptr;
		sprintf(LESdata,"Ek.dat");
		fptr = fopen(LESdata,"a");

		fprintf(fptr,"%f\t%f\n",step*deltaT,Ek0);

		fclose(fptr);


	// =============================================================================================================== //
	}

}


}
