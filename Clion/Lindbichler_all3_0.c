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
int conv(int itr);
void output();
void timestep();

/*------------------ DECLARATIONS -------------------------*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef fabs
#define fabs(a) (((a) < (0.0)) ? (-(a)) : (a))
#endif

/* max. dimension of the arrays*/
#define Mat_dim 1001

/*------------------------- INPUT VARIABLES----------------------------*/
int 	iread;  			/*0:std, 1:read from input file */
int 	max_iter;	 		/*max number of iterations*/
int 	nprint; 			/*print at xxth iteration*/
double	convergence;		/*convergence limit*/
int 	imax;				/*number of mesh points*/
double	x_min;	 			/*[m]*/
double	x_max;	 			/*[m]*/
double	y_min;				/*[m] for laval nozzle*/
double	y_max;				/*[m] for laval nozzle*/
double	R;					/*R*/
double	gamma;   			/*kappa*/
double	p_tot; 				/*total inlet pressure [bar]*/
double	p_exit;				/*static exit pressure [bar]*/
double	T_tot;				/*total inlet temperature [K]*/
int		sub_exit;			/*sub/supersonic conditions*/
int		method;				/*=1: MCC, =2: LW, =3:central*/
int		simple_dissip;		/*1: simple dissipation, else sophisticated dissipation*/
double	cfl;   				/*courant number   <0 => transient*/
double	eps_s;				// Parameter of simple dissipation
double	k4;					// 4th order parameter for complex dissipation
double	k2;					// 2nd order parameter for complex dissipation
/*--------------------------   END INPUT VARIABLES  ----------------------*/

/*------------------------- STATIC VARIABLES----------------------------*/
/*Files*/
FILE *logfile;
/*control*/
double dt;
double dx ;
double time ;
/*geometry*/
double area[Mat_dim] = {0.0};
double da_dx[Mat_dim] = {0.0};
/*vektor*/
double u[Mat_dim][3] = {0.0};
double u_q[Mat_dim][3] = {0.0};
double u_qq[Mat_dim][3] = {0.0};
double dissip[Mat_dim][3] = {0.0};
double delta_u[Mat_dim][3] = {0.0};
double f[Mat_dim][3] = {0.0};
double f_star[Mat_dim][3] = {0.0};
double source[Mat_dim][3] = {0.0};
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
	int itr, end=0,i,k;

	input();

	grid();

	init();

	//----------------loop start----------------------------------------
	for (itr=1; itr<=max_iter; itr++)
	{
        timestep();

        switch (method)
		{
			case 1: //MCC

				calc_uq();

				boundary_q();

				calc_uqq();

				for (i=1; i<imax-1; i++)
				{
					for (k=0;k<3;k++)
					{
						delta_u[i][k] = 0.5*(u_q[i][k] + u_qq[i][k])-u[i][k];
						u[i][k] = u[i][k] + delta_u[i][k];
					}
				}

				break;

			case 2: //LW

				calc_f();

      			calc_f_star_LW();

				for (i=1; i<imax-1; i++)
				{
					for (k=0;k<3;k++)
					{
						//delta_u[i][k] = ...  //  Calculation of delta_u using the conservative star fluxes
						u[i][k] = u[i][k] + delta_u[i][k];
					}
				}

				break;

			case 3: //central

				calc_f();

				if (simple_dissip == 1)
					dissip_simple();
				else
					dissip_complex();

      			calc_f_star_central();

                for (i=1; i<imax-1; i++)
                {
                    for (k=0;k<3;k++)
                    {
                        delta_u[i][k] = -dt/dx*(f_star[i][k]-f_star[i-1][k]) + dt*source[i][k];
                        u[i][k] = u[i][k] + delta_u[i][k];
                    }
                }

				break;

			default:

				printf("\n No numerical method was selected (1-3)!!!\n");
				exit(0);

		}



        boundary();

        if((itr%nprint == 0) || (itr == max_iter))		end = conv(itr);

        if (end == 1)
		{
			printf("\n Convergence limit of %lf achieved! \n",convergence);
			printf("\n Number of iterations: %d \n",itr);
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
	char  dummy[81];

	input=fopen("nozzle_3.dat","rt");
	fgets(dummy,80,input);
	fgets(dummy,80,input);
	fgets(dummy,80,input);
	fgets(dummy,80,input);
	fscanf(input,"%d%d%d%lf\n",&iread,&max_iter,&nprint,&convergence);
	fgets(dummy,80,input);
	fscanf(input,"%d%lf%lf%lf%lf\n",&imax,&x_min,&x_max,&y_min,&y_max);
	fgets(dummy,80,input);
	fscanf(input,"%lf%lf\n",&R,&gamma);
	fgets(dummy,80,input);
	fscanf(input,"%lf%lf%lf%d\n",&p_tot,&T_tot,&p_exit,&sub_exit);
	fgets(dummy,80,input);
	fscanf(input,"%d%d\n",&method, &simple_dissip);
	fgets(dummy,80,input);
	fscanf(input,"%lf%lf%lf%lf\n",&cfl,&eps_s, &k4, &k2);

	if (imax > Mat_dim)
	{
		printf("** Specified number of points is higher than allocated vector size ! **\n");
	}

	p_tot=p_tot*100000.0;
	p_exit=p_exit*100000.0;
	h_tot = gamma/(gamma-1)*R*T_tot;
}


/*--------------------------grid definition -----------------------------------*/
void grid()
{
    int i;
    // Berechnen von dx, x[i], area[i], da_dx[i]
    // Calculation of dx, x[i], area[i], da_dx[i]
    dx = (x_max - x_min) / (imax - 1);
    for (i=0; i<imax; i++){
        x[i] = x_min + dx * i;
        area[i] = y_min + (y_max - y_min) * pow(x[i],2) / pow(x_max,2);
        da_dx[i] = (y_max-y_min)*x[i]*2./pow(x_max,2);
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
		// Initialisation of the flow field (state vector U) with the stagnation values (=total values)

		for(i=0; i<=(imax-1); i++) {

		    u[i][0] = rho_tot;
		    u[i][1] = 0.0;
		    u[i][2] = p_tot / (gamma - 1);

		}

	}
	else
	{
		old_data = fopen("nozzle.out","r");

		for(i=0; i<imax; i++)
			fscanf(old_data,"%lf%lf%lf\n",&u[i][0],&u[i][1],&u[i][2]);

		fclose(old_data);
	}


    logfile = fopen("nozzle.log","w");
	fprintf (logfile, "%s\n", " Iter cont_resid imp_resid energy_resid");
}


//--------------------------timestep calculation--------------------------------------
void timestep()
{
	int i;
	double eigenmax,vel,p,c,eigen,T;


/*Calculation of the maximum eigenvalue eigenmax for the total flow field
         eigen = max(fabs(vel+c),fabs(vel-c))

Find dt as function of cfl and maximium eigenvalue

*/

    eigenmax = 0;

    for(i=0; i<=(imax-1); i++) {

        vel = u[i][1] / u[i][0];

        T = (u[i][2] / u[i][0] - pow(vel, 2) / 2) * (gamma - 1) / R;
        c = pow(gamma * R * T, 0.5);

        eigen = max(fabs(vel+c),fabs(vel-c));

        eigenmax = max(eigen, eigenmax);
    }

    dt = cfl * dx / eigenmax;

	time = time + dt;
}

//-----------------------flux and source vector------------------------------------------
void calc_f()
{
    double rho,vel,p, temp;
	int i;

	//Berechnung des des Flussvektors F und des Source-Vektors in allen Punkten
	//Calculaton of flux vector F and source vector S in all grid points

}

//------------------------simple dissipation vector---------------------------------------
void dissip_simple()
{
	int i,k;

	// dissipation vector at i+1/2


	// dissipation vector at i=0 and i=imax-2


}


//-----------------------complex dissipation vector-----------------------------------------
void dissip_complex()
{
	int i,k;
	double rho,vel,c,p,pm,pp;
	double eigen, eigenp;
	double eps2,eps4;
	double sensor[Mat_dim];

	// dissipation vector at i+1/2


	// dissipation vector for i=0 and i=imax-2


}


//---------------------------f_star central------------------------------------------
void calc_f_star_central()
{
   	int i,k;

	// calculation of f_star at i+1/2 for central method


}

//--------------------------f_star Lax-Wendroff-----------------------------------------
void calc_f_star_LW()
{
   	int i,k;

	// calculation of f_star at i+1/2 for Lax-Wendroff method

}

//--------------------------MCC U_q vector--------------------------------------------
void calc_uq()
{
    double rho,vel,p;
	int i;

	/*

	Calculaton of flux vector F and source vector S in all grid points for U
	Calculate U_q (forward)
	*/

	for(i=1; i<=(imax-2); i++) {

        rho = u[i][0];
        vel = u[i][1] / rho;

        p = (gamma - 1) * (u[i][2] - rho * pow(vel, 2) / 2);

        f[i][0] = rho * vel;
        f[i][1] = rho * pow(vel, 2) + p;
        f[i][2] = vel * (u[i][2] + p);



        source[i][0] = da_dx[i] / area[i] * f[i][0];
        source[i][1] = da_dx[i] / area[i] * rho * pow(vel, 2);
        source[i][2] = da_dx[i] / area[i] * f[i][2];

        u_q[i][0] = u[i][0] - dt / dx * (f[i+1][0]-f[i][0]) - source[i][0];
        u_q[i][1] = u[i][1] - dt / dx * (f[i+1][1]-f[i][1]) - source[i][1];
        u_q[i][2] = u[i][2] - dt / dx * (f[i+1][2]-f[i][2]) - source[i][2];
    }


}

//----------------------------MCC_U_qq vector------------------------------------------
void calc_uqq()
{
    double rho,vel,p;
	int i;

	/*
	Berechnung des des Flussvektors F und des Source-Vektors in allen Punkten fuer U_q
	Bestimmung von U_qq (backward)

	Calculaton of flux vector F and source vector S in all grid points for U_q
	Calculate U_qq (backward)
	*/

	for(i=1; i<=(imax-2); i++) {

        rho = u_q[i][0];
        vel = u_q[i][1] / rho;

        p = (gamma - 1) * (u[i][2] - rho * pow(vel, 2) / 2);

        f[i][0] = rho * vel;
        f[i][1] = rho * pow(vel, 2) + p;
        f[i][2] = vel * (u[i][2] + p);

        source[i][0] = da_dx[i] / area[i] * f[i][0];
        source[i][1] = da_dx[i] / area[i] * rho * pow(vel, 2);
        source[i][2] = da_dx[i] / area[i] * f[i][2];

        u_qq[i][0] = u[i][0] - dt / dx * (f[i][0]-f[i-1][0]) - source[i][0];
        u_qq[i][1] = u[i][1] - dt / dx * (f[i][1]-f[i-1][1]) - source[i][1];
        u_qq[i][2] = u[i][2] - dt / dx * (f[i][2]-f[i-1][2]) - source[i][2];
    }

}

//-------------------------------U_q boundary conditions-----------------------------------
void boundary_q()
{
    double rho,p,vel,e, rho_vel;

//	Calculation of boundary values for i=0 and i=imax-1 for U_q vector


	/*inlet i=0*/

	rho = 2 * u_q[1][0] - u_q[2][0];

	if (rho > rho_tot){

        rho = rho_tot;
	}

	p = p_tot * pow((R * T_tot * rho / p_tot), gamma);

    vel = pow(2 * gamma * R * T_tot / (gamma - 1) * (1 - pow(p / p_tot, ((gamma - 1)/gamma))),0.5);

    e = p / (gamma - 1) * pow(vel, 2) * rho / 2;

    u_q[0][0] = rho;
    u_q[0][1] = rho * vel;
    u_q[0][2] = e;

	/*outlet i=imax-1*/

	rho = 2 * u_q[imax-2][0] - u_q[imax-3][0];
    rho_vel = 2 * u_q[imax-2][1] - u_q[imax-3][1];

    if (sub_exit == 1){

       vel = rho_vel / rho;

       e = p_exit / (gamma-1) + rho * pow(vel,2) / 2;
    }
    else
    {
        e = 2 * u_q[imax-2][2] - u_q[imax-3][2];
    }

    u_q[imax-1][0] = rho;
    u_q[imax-1][1] = rho_vel;
    u_q[imax-1][2] = e;
}

//----------------------------U boundary conditions------------------------------------------
void boundary()
{
    double rho,p,vel,temp,e,rho_vel;

	//	Bestimmen der Randwerte fuer i=0 und i=imax-1 fuer U-Vektor
	//	Calculation of boundary values for i=0 and i=imax-1 for U vector


	/*inlet i=0*/

	rho = 2 * u[1][0] - u[2][0];

	if (rho > rho_tot){

        rho = rho_tot;
	}

	p = p_tot * pow((R * T_tot * rho / p_tot), gamma);

    vel = pow(2 * gamma * R * T_tot / (gamma - 1) * (1 - pow(p / p_tot, ((gamma - 1)/gamma))),0.5);

    e = p / (gamma - 1) + pow(vel, 2) * rho / 2;

    u[0][0] = rho;
    u[0][1] = rho * vel;
    u[0][2] = e;

	/*outlet i=imax-1*/

	rho = 2 * u[imax-2][0] - u[imax-3][0];
    rho_vel = 2 * u[imax-2][1] - u[imax-3][1];

    if (sub_exit == 1){

       vel = rho_vel / rho;

       e = p_exit / (gamma-1) + rho * pow(vel,2) / 2;
    }
    else
    {
        e = 2 * u[imax-2][2] - u[imax-3][2];
    }

    u[imax-1][0] = rho;
    u[imax-1][1] = rho_vel;
    u[imax-1][2] = e;
}

//----------------------------------------------------------------------
int conv(int itr)
{
	int i,k,end;
	double resid[3] = {0.0};

	if (itr == 1) return 0;

	printf("calc timestep = %d\t", itr);


	for(i=1; i<imax-1; i++)
	{
        for (k=0;k<3;k++)
		{
			resid[k] = resid[k] + fabs(delta_u[i][k]);
		}
	}

	if ((itr == nprint) || (itr==2))
	{
		resid0 = resid[0];
		resid1 = resid[1];
		resid2 = resid[2];
	}

	fprintf(logfile,"%d %lf %lf %lf\n",itr, resid[0]/resid0, resid[1]/resid1, resid[2]/resid2);

	if (resid[2]/resid2 < convergence) end = 1;

	printf("resid = %lf\n", resid[2]/resid2);

    return end;
}



//----------------------------------------------------------------------
void output()
{
    int i;
	double temp,cont0, rho,vel;
	double p[Mat_dim], mach[Mat_dim],p_tot_is[Mat_dim], cont[Mat_dim];
	FILE *result,*machzahl,*pressure,*continuity,*total_pressure;

    result = fopen("nozzle.out","wt");
	machzahl = fopen("mach.out","wt");
	pressure = fopen("press.out","wt");
	continuity = fopen("cont.out","wt");
	total_pressure = fopen("ptot.out","wt");

	for (i=0;i<imax;i++) fprintf(result,"%lf\t%lf\t%lf\n",u[i][0],u[i][1],u[i][2]);
	fprintf(result,"\n %lf \n", time);

	for (i=0; i<imax; i++)
	{
        rho = u[i][0];
        vel = u[i][1]/rho;
        p[i] = (gamma-1)*(u[i][2]-rho*vel*vel/2);
        mach[i] = vel/pow((gamma*p[i]/rho),0.5);
        temp = 1+(gamma-1)/2*mach[i]*mach[i];
        p_tot_is[i] = p[i]*pow(temp,(gamma/(gamma-1)));
        if (i == 0)
		{
            cont0 = rho*vel*area[i];
            if(cont0<=0.) cont0 = 1.e-5;
		}
        cont[i] = rho*vel*area[i]/cont0;
	}


    for (i=0;i<imax;i++) fprintf(pressure,"%lf\t%lf\n", x[i], p[i]/1.e5);
    for (i=0;i<imax;i++) fprintf(machzahl,"%lf\t%lf\n", x[i], mach[i]);
    for (i=0;i<imax;i++) fprintf(continuity,"%lf\t%lf\n", x[i], cont[i]);
    for (i=0;i<imax;i++) fprintf(total_pressure,"%lf\t%lf\n", x[i], p_tot_is[i]/p_tot);

    fclose (pressure);
    fclose (machzahl);
    fclose (continuity);
    fclose (total_pressure);
    fclose (result);
}

