



/**** computational parameters ****/

#define e 0.0001 //0.005

#define K 1.4
#define R 287 


double f = 0.0;  

double Q0 = 0.0;
double Qold = Q0;

double L = deltaXI*X_out/2.0/3.1415928;

double high = 2.0*3.1415*L;

double psi = 0.001;

double gamma = 1.0e-8;

/**** computational parameters-end ****/

double rho0 = 1.1842;
double U0 = 34.6064;
double V0 = 0.0;
double W0 = 0.0;
double P0 = 101300.0;
double T0 = P0/(R*T0);