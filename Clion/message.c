/*      program CENTRAL
c***********************************************************************
c
c                           1D CFD code for
c
c                             students of
c                     CFD in turbomachinery
c
c
c     solution of the 1-D unsteady Euler equations for nozzle flows
c     with a time-iterative central method
c
c     numerical algorithm:
c     - MacCormack algorithm
c	  - Lax-Wendroff algorithm
c	  - explicit central method with simple and "sophisticated" 4th order
c	    artifical dissipation
c
c	  Remarks to central method: CFL < 0.1
c	  values of eps, k2, k4 Parameter are important
c	  good results for eps=50, k4=0.04,k2=1
c	  (if k2 is too small, Dissip4 is not switched off)
c
c	  index 0 is included!!
c	  imax points, so (imax-1) cells
c
c***********************************************************************/

/*------------------ SUBROUTINES -------------------------*/
void init();
void init_shocktube();
void grid();
void input();
void boundary();
void boundary_q();
void calc_uq();
void calc_uqq();
void calc_f();
void dissip_simple();
void dissip_complex();
void calc_f_star_LW();
void calc_f_star_central();
void calc_f_star_roe();
int conv(int itr);
void output();
void timestep();
void mmult();

/*------------------ DECLARATIONS -------------------------*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef max
#define max(a, b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef fabs
#define fabs(a) (((a) < (0.0)) ? (-(a)) : (a))
#endif

/* max. dimension of the arrays*/
#define Mat_dim 1001

/*------------------------- INPUT VARIABLES----------------------------*/
int iread;			/*0:std, 1:read from input file */
int max_iter;		/*max number of iterations*/
int nprint;			/*print at xxth iteration*/
double convergence; /*convergence limit*/
int imax;			/*number of mesh points*/
double x_min;		/*[m]*/
double x_max;		/*[m]*/
double y_min;		/*[m] for laval nozzle*/
double y_max;		/*[m] for laval nozzle*/
double R;			/*R*/
double gamma;		/*kappa*/
double p_tot;		/*total inlet pressure [bar]*/
double p_exit;		/*static exit pressure [bar]*/
double T_tot;		/*total inlet temperature [K]*/
int sub_exit;		/*sub/supersonic conditions*/
int method;			/*=1: MCC, =2: LW, =3:central*/
int simple_dissip;	/*1: simple dissipation, else sophisticated dissipation*/
double cfl;			/*courant number   <0 => transient*/
double eps_s;		// Parameter of simple dissipation
double k4;			// 4th order parameter for complex dissipation
double k2;			// 2nd order parameter for complex dissipation
/*--------------------------   END INPUT VARIABLES  ----------------------*/

/*------------------------- STATIC VARIABLES----------------------------*/
/*Files*/
FILE *logfile;
/*control*/
double dt;
double dx;
double time;
/*geometry*/
double area[Mat_dim] = {0.0};
double da_dx[Mat_dim] = {0.0};
/*vektor*/
double u[Mat_dim][3] = {{0.0}};
double u_q[Mat_dim][3] = {{0.0}};
double u_qq[Mat_dim][3] = {{0.0}};
double dissip[Mat_dim][3] = {{0.0}};
double delta_u[Mat_dim][3] = {{0.0}};
double f[Mat_dim][3] = {{0.0}};
double f_star[Mat_dim][3] = {{0.0}};
double source[Mat_dim][3] = {{0.0}};
/*boundary*/
double rho_tot;
double h_tot;
/*x_coordinaten*/
double x[Mat_dim] = {0.0};
/*residuum*/
double resid0;
double resid1;
double resid2;
/*------------------------- END STATIC VARIABLES----------------------------*/

void main()
{
	int itr, end = 0, i, k;

	input();
	grid();
	init();
	//init_shocktube();

	//----------------loop start----------------------------------------
	for (itr = 1; itr <= max_iter; itr++)
	{
		timestep();

		switch (method)
		{
		case 1: // MCC
			calc_uq();
			boundary_q();
			calc_uqq();
			for (i = 1; i < imax - 1; i++)
			{
				for (k = 0; k < 3; k++)
				{
					delta_u[i][k] = 0.5 * (u_q[i][k] + u_qq[i][k]) - u[i][k];
					u[i][k] = u[i][k] + delta_u[i][k];
				}
			}
			break;

		case 2: // LW
			calc_f();
			calc_f_star_LW();
			for (i = 1; i < imax - 1; i++)
			{
				for (k = 0; k < 3; k++)
				{
					// delta_u[i][k] = ...  //  Calculation of delta_u using the conservative star fluxes
					u[i][k] = u[i][k] + delta_u[i][k];
				}
			}
			break;

		case 3: // central

			calc_f();

			if (simple_dissip == 1)
				dissip_simple();
			else
				dissip_complex();

			calc_f_star_central();

			for (i = 1; i < imax - 1; i++)
			{
				for (k = 0; k < 3; k++)
				{
					//  Calculation of delta_u using the conservative star fluxes
					delta_u[i][k] = -dt / dx * (f_star[i][k] - f_star[i - 1][k]) + dt * source[i][k];
					u[i][k] = u[i][k] + delta_u[i][k];
				}
			}
			break;

		case 4: // roe

			calc_f();
			calc_f_star_roe();

			for (i = 1; i < imax - 1; i++)
			{
				for (k = 0; k < 3; k++)
				{
					//  Calculation of delta_u using the conservative star fluxes
					delta_u[i][k] = -dt / dx * (f_star[i][k] - f_star[i - 1][k]) + dt * source[i][k];
					u[i][k] = u[i][k] + delta_u[i][k];
				}
			}
			break;

		default:
			printf("\n No numerical method was selected (1-3)!!!\n");
			exit(0);
		}

		boundary();

		if ((itr % nprint == 0) || (itr == max_iter))
			end = conv(itr);

		if (end == 1)
		{
			printf("\n Convergence limit of %lf achieved! \n", convergence);
			printf("\n Number of iterations: %d \n", itr);
			break;
		}
	}

	//---------------------loop end-----------------------------------------

	output();
	fclose(logfile);
	return;
}

/*--------------------------data input-----------------------------------*/
void input()
{
	FILE *input;
	char dummy[81];

	input = fopen("nozzle_3.dat", "rt");
	fgets(dummy, 80, input);
	fgets(dummy, 80, input);
	fgets(dummy, 80, input);
	fgets(dummy, 80, input);
	fscanf(input, "%d%d%d%lf\n", &iread, &max_iter, &nprint, &convergence);
	fgets(dummy, 80, input);
	fscanf(input, "%d%lf%lf%lf%lf\n", &imax, &x_min, &x_max, &y_min, &y_max);
	fgets(dummy, 80, input);
	fscanf(input, "%lf%lf\n", &R, &gamma);
	fgets(dummy, 80, input);
	fscanf(input, "%lf%lf%lf%d\n", &p_tot, &T_tot, &p_exit, &sub_exit);
	fgets(dummy, 80, input);
	fscanf(input, "%d%d\n", &method, &simple_dissip);
	fgets(dummy, 80, input);
	fscanf(input, "%lf%lf%lf%lf\n", &cfl, &eps_s, &k4, &k2);

	if (imax > Mat_dim)
	{
		printf("** Specified number of points is higher than allocated vector size ! **\n");
	}

	p_tot = p_tot * 100000.0;
	p_exit = p_exit * 100000.0;
	h_tot = gamma / (gamma - 1) * R * T_tot;
}

/*--------------------------grid definition -----------------------------------*/
void grid()
{
	int i;
	// Berechnen von dx, x[i], area[i], da_dx[i]
	// Calculation of dx, x[i], area[i], da_dx[i]8
	dx = (x_max - x_min) / (imax - 1.0);

	for (i = 0; i < imax; i++)
	{
		x[i] = x_min + dx * i;
		area[i] = y_min + (y_max - y_min) * x[i] * x[i] / (x_max * x_max);
		da_dx[i] = 2.0 * (y_max - y_min) * x[i] / (x_max * x_max);
	}
}

/*--------------------------initializing-----------------------------------*/
void init()
{
	FILE *old_data;
	int i;

	// Berechnen von rho_tot, am Eintritt fuer die Randbedingungen
	// Calculaton of rho_tot at the inlet for the boundary condition algorithm
	rho_tot = p_tot / (R * T_tot);

	if (iread == 0)
	{
		// Initialisieren des Stroemungsfeldes (Zustandsvektor U) mit den Ruhezustandswerten
		// Initialisation of the flow field (state vector U) with the stagnation values (=total values)
		for (i = 0; i < imax; i++)
		{
			u[i][0] = rho_tot;
			u[i][1] = 0.0;
			u[i][2] = p_tot / (gamma - 1);
		}
	}
	else
	{
		old_data = fopen("nozzle.out", "r");

		for (i = 0; i < imax; i++)
			fscanf(old_data, "%lf%lf%lf\n", &u[i][0], &u[i][1], &u[i][2]);

		fclose(old_data);
	}

	logfile = fopen("nozzle.log", "w");
	fprintf(logfile, "%s\n", " Iter cont_resid imp_resid energy_resid");
}

void init_shocktube()
{
	FILE *old_data;
	int i;

	// Berechnen von rho_tot, am Eintritt fuer die Randbedingungen
	// Calculaton of rho_tot at the inlet for the boundary condition algorithm
	rho_tot = p_tot / (R * T_tot);

	if (iread == 0)
	{
		// Initialisieren des Stroemungsfeldes (Zustandsvektor U) mit den Ruhezustandswerten
		// Initialisation of the flow field (state vector U) with the stagnation values (=total values)
		for (i = 0; i < imax / 2; i++)
		{
			u[i][0] = rho_tot;
			u[i][1] = 0.0;
			u[i][2] = p_tot / (gamma - 1);
		}
		for (i = imax / 2; i < imax; i++)
		{
			u[i][0] = rho_tot / 8.0;
			u[i][1] = 0.0;
			u[i][2] = p_tot / 10.0 / (gamma - 1);
		}
	}
	else
	{
		old_data = fopen("nozzle.out", "r");

		for (i = 0; i < imax; i++)
			fscanf(old_data, "%lf%lf%lf\n", &u[i][0], &u[i][1], &u[i][2]);

		fclose(old_data);
	}

	logfile = fopen("nozzle.log", "w");
	fprintf(logfile, "%s\n", " Iter cont_resid imp_resid energy_resid");
}

//--------------------------timestep calculation--------------------------------------
void timestep()
{
	int i;
	double eigenmax = 0.0;
	double vel, p, c, eigen;

	for (i = 1; i < imax; i++)
	{
		double rho = u[i][0];
		vel = u[i][1] / rho;
		p = (u[i][2] - rho * vel * vel / 2.0) * (gamma - 1.0);
		c = pow(gamma * p / rho, 0.5);
		eigen = max(fabs(vel + c), fabs(vel - c));
		if (eigen > eigenmax)
		{
			eigenmax = eigen;
		}
	}

	dt = cfl * dx / eigenmax;

	/*
Bestimmen des maximalen Eigenwertes eigenmax fuer das gesamte Stroemungsfeld
		 eigen = max(fabs(vel+c),fabs(vel-c))
Calculation of the maximum eigenvalue eigenmax for the total flow field
		 eigen = max(fabs(vel+c),fabs(vel-c))

Bestimmen von dt als Funktion von cfl und max. Eigenwert
Find dt as function of cfl and maximium eigenvalue

*/
	time = time + dt;
}

//-----------------------flux and source vector------------------------------------------
void calc_f()
{
	double rho, vel, p;
	int i;

	for (i = 0; i < imax; i++)
	{
		rho = u[i][0];
		vel = u[i][1] / rho;
		p = (u[i][2] - rho * vel * vel / 2.0) * (gamma - 1.0);

		f[i][0] = rho * vel;
		f[i][1] = rho * vel * vel + p;
		f[i][2] = (u[i][2] + p) * vel;
		source[i][0] = -da_dx[i] / area[i] * f[i][0];
		source[i][1] = -da_dx[i] / area[i] * rho * vel * vel;
		source[i][2] = -da_dx[i] / area[i] * f[i][2];
	}

	// Berechnung des des Flussvektors F und des Source-Vektors in allen Punkten
	// Calculaton of flux vector F and source vector S in all grid points
}

//------------------------simple dissipation vector---------------------------------------
void dissip_simple()
{
	int i, k;

	// dissipation vector at i+1/2
	for (i = 1; i < imax - 2; i++)
	{
		for (k = 0; k < 3; k++)
		{
			dissip[i][k] = eps_s * dx * (u[i + 2][k] - 3.0 * u[i + 1][k] + 3.0 * u[i][k] - u[i - 1][k]);
		}
	}

	// dissipation vector at i=0 and i=imax-2
	for (k = 0; k < 3; k++)
	{
		dissip[0][k] = eps_s * dx * (u[2][k] - 2.0 * u[1][k] + u[0][k]);
		dissip[imax - 2][k] = eps_s * dx * (-u[imax - 1][k] + 2.0 * u[imax - 2][k] - u[imax - 3][k]);
	}
}

//-----------------------complex dissipation vector-----------------------------------------
void dissip_complex()
{
	// int i,k;
	// double rho,vel,c,p,pm,pp;
	// double eigen, eigenp;
	// double eps2,eps4;
	// double sensor[Mat_dim];

	// dissipation vector at i+1/2

	// dissipation vector for i=0 and i=imax-2
}

//---------------------------f_star central------------------------------------------
void calc_f_star_central()
{
	int i, k;

	// calculation of f_star at i+1/2 for central method
	for (i = 0; i < imax - 1; i++)
	{
		for (k = 0; k < 3; k++)
		{
			f_star[i][k] = 0.5 * (f[i + 1][k] + f[i][k]) + dissip[i][k];
		}
	}
}

void calc_f_star_roe()
{
	int i, k;
	double R_avg, rho_avg, u_avg, H_avg, c_avg;
	double H_i, H_ip1, p_i, p_ip1;  // p ... plus
	double lambda_avg_1, lambda_avg_2, lambda_avg_3;
	double c_i, c_ip1;
	double lambda_i_1, lambda_i_2, lambda_i_3;
	double lambda_ip1_1, lambda_ip1_2, lambda_ip1_3;
	double epsilon_1, epsilon_2, epsilon_3;
	double lambda_mod_1, lambda_mod_2, lambda_mod_3;
	double R_eig[3][3] = {{0.0}};	   // Right eigenvector
	double L_eig[3][3] = {{0.0}};	   // Left eigenvector
	double Lambda_tmp[3][3] = {{0.0}}; // Eigenvalue matrix
	double _R_Lambda_tmp[3][3] = {{0.0}};
	double _tmp_eig[3][3] = {{0.0}};
	double _tmp_u[3][1] = {{0.0}};
	double _tmp[3][1] = {{0.0}};

	for (i = 0; i < imax - 1; i++)
	{
		R_avg = sqrt(u[i + 1][0] / u[i][0]);										   // Eq. 11-75
		rho_avg = R_avg * u[i][0];													   // Eq. 11-76
		u_avg = (R_avg * u[i + 1][1] / u[i + 1][0] + u[i][1] / u[i][0]) / (R_avg + 1); // Eq. 11-77

		p_i = (u[i][2] - u[i][1] * u[i][1] / u[i][0] / 2.0) * (gamma - 1.0);
		p_ip1 = (u[i + 1][2] - u[i + 1][1] * u[i + 1][1] / u[i + 1][0] / 2.0) * (gamma - 1.0);
		H_i = gamma / (gamma - 1) * p_i / u[i][0] + u[i][1] * u[i][1] / u[i][0] / u[i][0];
		H_ip1 = gamma / (gamma - 1) * p_ip1 / u[i + 1][0] + u[i + 1][1] * u[i + 1][1] / u[i + 1][0] / u[i + 1][0];
		H_avg = (R_avg * H_ip1 + H_i) / (R_avg + 1.0); // Eq. 11-78

		c_avg = sqrt((gamma - 1.0) * (H_avg - u_avg * u_avg / 2.0)); // Eq. 11-80

		// Calc Eigenvalues according to Eq. 11-81
		lambda_avg_1 = u_avg;
		lambda_avg_2 = u_avg + c_avg;
		lambda_avg_3 = u_avg - c_avg;

		c_i = pow(gamma * p_i / u[i][0], 0.5);
		lambda_i_1 = u[i][1] / u[i][0];
		lambda_i_2 = lambda_i_1 + c_i;
		lambda_i_3 = lambda_i_1 - c_i;

		c_ip1 = pow(gamma * p_ip1 / u[i+1][0], 0.5);
		lambda_ip1_1 = u[i+1][1] / u[i+1][0];
		lambda_ip1_2 = lambda_ip1_1 + c_ip1;
		lambda_ip1_3 = lambda_ip1_1 - c_ip1;

		epsilon_1 = max(max(lambda_avg_1 - lambda_i_1, lambda_ip1_1 - lambda_avg_1), 0);
		epsilon_2 = max(max(lambda_avg_2 - lambda_i_2, lambda_ip1_2 - lambda_avg_2), 0);
		epsilon_3 = max(max(lambda_avg_3 - lambda_i_3, lambda_ip1_3 - lambda_avg_3), 0);

		lambda_mod_1 = max(fabs(lambda_avg_1), epsilon_1);
		lambda_mod_2 = max(fabs(lambda_avg_2), epsilon_2);
		lambda_mod_3 = max(fabs(lambda_avg_3), epsilon_3);

		L_eig[0][0] = 1.0 - (gamma - 1.0) / 2.0 * u_avg * u_avg / c_avg / c_avg;
		L_eig[0][1] = (gamma - 1.0) * u_avg / c_avg / c_avg;
		L_eig[0][2] = -(gamma - 1.0) / c_avg / c_avg;
		L_eig[1][0] = ((gamma - 1.0) / 2.0 * u_avg * u_avg - u_avg * c_avg) / rho_avg / c_avg;
		L_eig[1][1] = (c_avg - (gamma - 1.0) * u_avg) / rho_avg / c_avg;
		L_eig[1][2] = (gamma - 1.0) / rho_avg / c_avg;
		L_eig[2][0] = ((gamma - 1.0) / 2.0 * u_avg * u_avg + u_avg * c_avg) / rho_avg / c_avg;
		L_eig[2][1] = -(c_avg + (gamma - 1.0) * u_avg) / rho_avg / c_avg;
		L_eig[2][2] = L_eig[1][2];

		R_eig[0][0] = 1.0;
		R_eig[0][1] = rho_avg / 2.0 / c_avg;
		R_eig[0][2] = R_eig[0][1];
		R_eig[1][0] = lambda_avg_1;
		R_eig[1][1] = R_eig[0][1] * lambda_avg_2;
		R_eig[1][2] = R_eig[0][1] * lambda_avg_3;
		R_eig[2][0] = u_avg * u_avg / 2.0;
		R_eig[2][1] = R_eig[0][1] * (c_avg * c_avg / (gamma - 1.0) + R_eig[2][0] + u_avg * c_avg);
		R_eig[2][2] = R_eig[0][1] * (c_avg * c_avg / (gamma - 1.0) + R_eig[2][0] - u_avg * c_avg);

		Lambda_tmp[0][0] = (lambda_avg_1 - lambda_mod_1) / 2.0;
		Lambda_tmp[1][1] = (lambda_avg_2 - lambda_mod_2) / 2.0;
		Lambda_tmp[2][2] = (lambda_avg_3 - lambda_mod_3) / 2.0;

		for (k = 0; k < 3; k++)
		{
			_tmp_u[k][0] = u[i + 1][k] - u[i][k];
		}

		mmult(3, 3, 3, R_eig, Lambda_tmp, _R_Lambda_tmp);
		mmult(3, 3, 3, _R_Lambda_tmp, L_eig, _tmp_eig);
		mmult(3, 3, 1, _tmp_eig, _tmp_u, _tmp);

		for (k = 0; k < 3; k++)
		{
			f_star[i][k] = f[i][k] + _tmp[k][0];
		}
	}
}

void mmult(int n, int m, int k, double A[n][m], double B[m][k], double C[n][k])
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < k; j++)
		{
			C[i][j] = 0;
			for (int p = 0; p < m; p++)
			{
				C[i][j] += A[i][p] * B[p][j];
			}
		}
	}
}

//--------------------------f_star Lax-Wendroff-----------------------------------------
void calc_f_star_LW()
{
	// int i,k;

	// calculation of f_star at i+1/2 for Lax-Wendroff method
}

//--------------------------MCC U_q vector--------------------------------------------
void calc_uq()
{
	// double rho,vel,p;
	// int i;

	/*
	Berechnung des Flussvektors F und des Source-Vektors in allen Punkten fuer U
	Bestimmung von U_q (forward)

	Calculaton of flux vector F and source vector S in all grid points for U
	Calculate U_q (forward)
	*/
}

//----------------------------MCC_U_qq vector------------------------------------------
void calc_uqq()
{
	// double rho,vel,p;
	// int i;

	/*
	Berechnung des des Flussvektors F und des Source-Vektors in allen Punkten fuer U_q
	Bestimmung von U_qq (backward)

	Calculaton of flux vector F and source vector S in all grid points for U_q
	Calculate U_qq (backward)
	*/
}

//-------------------------------U_q boundary conditions-----------------------------------
void boundary_q()
{
	// double p,vel;
	double rho;

	//	Bestimmen der Randwerte fuer i=0 und i=imax-1 fuer U_q-Vektor
	//	Calculation of boundary values for i=0 and i=imax-1 for U_q vector

	/*inlet i=0*/

	if (rho > rho_tot)
		rho = rho_tot;

	/*outlet i=imax-1*/
}

//----------------------------U boundary conditions------------------------------------------
void boundary()
{
	double p, vel, tmp;
	double rho;

	//	Bestimmen der Randwerte fuer i=0 und i=imax-1 fuer U-Vektor
	//	Calculation of boundary values for i=0 and i=imax-1 for U vector

	/*inlet i=0*/
	rho = u[1][0] - (u[2][0] - u[1][0]);
	if (rho > rho_tot)
	{
		rho = rho_tot;
	}

	//p = p_tot * pow((R * T_tot * rho / p_tot), gamma); // oder?: p_tot*pow(((gamma-1)/gamma*h_tot*rho/p_tot),gamma);
	//vel = sqrt(2.0 * gamma / (gamma - 1.0) * R * T_tot * (1.0 - pow((p / p_tot), (gamma - 1.0) / gamma)));

    p = p_tot*pow(((gamma-1)/gamma*h_tot*rho/p_tot),gamma);
	tmp = pow((p/p_tot),((gamma-1)/gamma));
	vel = pow(2*h_tot*(1-tmp),0.5);

	u[0][0] = rho;
	u[0][1] = rho * vel;
	u[0][2] = p / (gamma - 1.0) + rho * pow(vel, 2.0) / 2.0;

	/*outlet i=imax-1*/
	if (sub_exit == 1)
	{
	    u[imax - 1][0] = u[imax - 2][0] + (u[imax - 2][0] - u[imax - 3][0]);
        u[imax - 1][1] = u[imax - 2][1] + (u[imax - 2][1] - u[imax - 3][1]);
        vel = u[imax - 1][1] / u[imax - 1][0];
		u[imax - 1][2] = p_exit / (gamma - 1.0) + u[imax - 1][0] * pow(vel, 2.0) / 2.0;
	}
	else
	{
		//u[imax - 1][2] = 2 * u[imax - 2][2] - u[imax - 3][2];
		u[imax - 1][0] = u[imax - 2][0];
		u[imax - 1][1] = u[imax - 2][1];
		u[imax - 1][2] = u[imax - 2][2];
	}
}

//----------------------------------------------------------------------
int conv(int itr)
{
	int i, k, end;
	double resid[3] = {0.0};

	if (itr == 1)
		return 0;

	printf("calc timestep = %d\t", itr);

	for (i = 1; i < imax - 1; i++)
	{
		for (k = 0; k < 3; k++)
		{
			resid[k] = resid[k] + fabs(delta_u[i][k]);
		}
	}

	if ((itr == nprint) || (itr == 2))
	{
		resid0 = resid[0];
		resid1 = resid[1];
		resid2 = resid[2];
	}

	fprintf(logfile, "%d %lf %lf %lf\n", itr, resid[0] / resid0, resid[1] / resid1, resid[2] / resid2);

	if (resid[2] / resid2 < convergence)
		end = 1;

	printf("resid = %lf\n", resid[2] / resid2);

	return end;
}

//----------------------------------------------------------------------
void output()
{
	int i;
	double temp, cont0, rho, vel;
	double p[Mat_dim], mach[Mat_dim], p_tot_is[Mat_dim], cont[Mat_dim];
	FILE *result, *machzahl, *pressure, *continuity, *total_pressure;

	result = fopen("nozzle.out", "wt");
	machzahl = fopen("mach.out", "wt");
	pressure = fopen("press.out", "wt");
	continuity = fopen("cont.out", "wt");
	total_pressure = fopen("ptot.out", "wt");

	for (i = 0; i < imax; i++)
		fprintf(result, "%lf\t%lf\t%lf\n", u[i][0], u[i][1], u[i][2]);
	fprintf(result, "\n %lf \n", time);

	for (i = 0; i < imax; i++)
	{
		rho = u[i][0];
		vel = u[i][1] / rho;
		p[i] = (gamma - 1) * (u[i][2] - rho * vel * vel / 2);
		mach[i] = vel / pow((gamma * p[i] / rho), 0.5);
		temp = 1 + (gamma - 1) / 2 * mach[i] * mach[i];
		p_tot_is[i] = p[i] * pow(temp, (gamma / (gamma - 1)));
		if (i == 0)
		{
			cont0 = rho * vel * area[i];
			if (cont0 <= 0.)
				cont0 = 1.e-5;
		}
		cont[i] = rho * vel * area[i] / cont0;
	}

	for (i = 0; i < imax; i++)
		fprintf(pressure, "%lf\t%lf\n", x[i], p[i] / 1.e5);
	for (i = 0; i < imax; i++)
		fprintf(machzahl, "%lf\t%lf\n", x[i], mach[i]);
	for (i = 0; i < imax; i++)
		fprintf(continuity, "%lf\t%lf\n", x[i], cont[i]);
	for (i = 0; i < imax; i++)
		fprintf(total_pressure, "%lf\t%lf\n", x[i], p_tot_is[i] / p_tot);

	fclose(pressure);
	fclose(machzahl);
	fclose(continuity);
	fclose(total_pressure);
	fclose(result);
}
