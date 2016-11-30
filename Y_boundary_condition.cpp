




#include <stdlib.h> 
#include <mpi.h>
#include <omp.h>
#include <omp.h>
#include "Resolution.h"

extern int X_np;

void Y_boundary_condition
(
// ============================================================================ //
int myid,

double (*U1_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
// ============================================================================ //
)

{

#include "ijk.h"
#include "prm.h"

#include "MPI_prm.h"
#include "Mpi_division.h"
	
	double rho,U,V,W,VV,P,C,T,h,H;
	double temp;

//// ============================================ ////
		if (myid ==0) istart = 2;		          ////	
		else istart = 0;            			  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid]+1;  ////
		else iend = gend[myid]+3;				  ////
//// ============================================ ////

#pragma omp parallel for private(k,rho,U,V,W,VV,P,temp,T)
	for (i = istart; i <= iend; i++) {
		for (k = 2; k <= nz; k++) {
			
			
			U1_[i][1][k] = U1_[i][ny][k];
			U2_[i][1][k] = U2_[i][ny][k];
			U3_[i][1][k] = U3_[i][ny][k];
			U4_[i][1][k] = U4_[i][ny][k];
			U5_[i][1][k] = U5_[i][ny][k];


			U1_[i][nyy][k] = U1_[i][2][k];
			U2_[i][nyy][k] = U2_[i][2][k];
			U3_[i][nyy][k] = U3_[i][2][k];
			U4_[i][nyy][k] = U4_[i][2][k];
			U5_[i][nyy][k] = U5_[i][2][k];

			
		}
	}
#pragma omp barrier

}
