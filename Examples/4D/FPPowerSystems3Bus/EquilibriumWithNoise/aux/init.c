#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

static void ComputeElectricalPower(double theta1,double theta2,
                                   double E1,double E2,double E3,
                                   double *G,double *B,
                                   double *Pe1,double *Pe2,double *Pe3)
{
  double Eph[3][2];
  Eph[0][0] = E1*cos(theta1);   Eph[0][1] = E1*sin(theta1);
  Eph[1][0] = E2*cos(theta2);   Eph[1][1] = E2*sin(theta2);
  Eph[2][0] = E3;               Eph[2][1] = 0.0;

  double Y[3][3][2];
  int i,j;
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      Y[i][j][0] = G[i*3+j];
      Y[i][j][1] = B[i*3+j];
    }
  }

  double YEph[3][2];
  YEph[0][0] =   Y[0][0][0]*Eph[0][0] - Y[0][0][1]*Eph[0][1]
               + Y[0][1][0]*Eph[1][0] - Y[0][1][1]*Eph[1][1]
               + Y[0][2][0]*Eph[2][0] - Y[0][2][1]*Eph[2][1];
  YEph[0][1] =   Y[0][0][0]*Eph[0][1] + Y[0][0][1]*Eph[0][0]
               + Y[0][1][0]*Eph[1][1] + Y[0][1][1]*Eph[1][0]
               + Y[0][2][0]*Eph[2][1] + Y[0][2][1]*Eph[2][0];
  YEph[1][0] =   Y[1][0][0]*Eph[0][0] - Y[1][0][1]*Eph[0][1]
               + Y[1][1][0]*Eph[1][0] - Y[1][1][1]*Eph[1][1]
               + Y[1][2][0]*Eph[2][0] - Y[1][2][1]*Eph[2][1];
  YEph[1][1] =   Y[1][0][0]*Eph[0][1] + Y[1][0][1]*Eph[0][0]
               + Y[1][1][0]*Eph[1][1] + Y[1][1][1]*Eph[1][0]
               + Y[1][2][0]*Eph[2][1] + Y[1][2][1]*Eph[2][0];
  YEph[2][0] =   Y[2][0][0]*Eph[0][0] - Y[2][0][1]*Eph[0][1]
               + Y[2][1][0]*Eph[1][0] - Y[2][1][1]*Eph[1][1]
               + Y[2][2][0]*Eph[2][0] - Y[2][2][1]*Eph[2][1];
  YEph[2][1] =   Y[2][0][0]*Eph[0][1] + Y[2][0][1]*Eph[0][0]
               + Y[2][1][0]*Eph[1][1] + Y[2][1][1]*Eph[1][0]
               + Y[2][2][0]*Eph[2][1] + Y[2][2][1]*Eph[2][0];

  YEph[0][1] = - YEph[0][1];
  YEph[1][1] = - YEph[1][1];
  YEph[2][1] = - YEph[2][1];

  *Pe1 = Eph[0][0]*YEph[0][0] - Eph[0][1]*YEph[0][1];
  *Pe2 = Eph[1][0]*YEph[1][0] - Eph[1][1]*YEph[1][1];
  *Pe3 = Eph[2][0]*YEph[2][0] - Eph[2][1]*YEph[2][1];
}


int main()
{
  const double PI = 4.0*atan(1.0);
  int ferr,ierr=0, n, i, N[4]={1,1,1,1},ndims,offset;
  char ip_file_type[50];
  FILE *out, *in;

  /* reading grid inputs */
  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) {
    printf("Error: Input file \"solver.inp\" not found.\n");
    return(0);
  } else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in,"%s",word);
        if (!strcmp(word, "ndims")) fscanf(in,"%d",&ndims);
        else if (!strcmp(word, "size")) {
          fscanf(in,"%d",&N[0]);
          fscanf(in,"%d",&N[1]);
          fscanf(in,"%d",&N[2]);
          fscanf(in,"%d",&N[3]);
        } else if (!strcmp(word, "ip_file_type")) fscanf(in,"%s",ip_file_type);
      }
    } else printf("Error: Illegal format in solver.inp. Crash and burn!\n");
  }
  fclose(in);
  if (ndims != 4) {
    printf("Error: ndims is not 4 in solver.inp. this code is to generate 4D initial conditions\n");
    return(0);
  }
	printf("Grid:\t\t\t%d X %d X %d X %d\n",N[0],N[1],N[2],N[3]);
  /* check if ip_file_type is binary */
  if (strcmp(ip_file_type,"binary") && strcmp(ip_file_type,"bin")) {
    printf("Error: ip_file_type needs to be bin or binary.\n");
    return(0);
  }

  double pi = 4.0*atan(1.0);
  /* default values of model parameters */
  double Pm1_avg      = 0.8;
  double Pm2_avg      = 1.6;
  double Pmref_avg    = 0.79330781761651;
  double H1           = 3.20;
  double H2           = 6.40;
  double Href         = 13.60;
  double E1           = 1.01556070860155;
  double E2           = 1.0491099265981;
  double Eref         = 1.05623172878954;
  double omegaB       = 2*pi*60.0;
  
  double sigma[2][2];
  sigma[0][0]  = 0.0125;
  sigma[0][1]  = 0.0;
  sigma[1][0]  = 0.0;
  sigma[1][1]  = 0.0125;

  double lambda[2][2];
  lambda[0][0] = 10.0/omegaB;
  lambda[0][1] = 0.0;
  lambda[1][0] = 0.0;
  lambda[1][1] = 10.0/omegaB;

  double gamma        = 0.25;

  double G[9], B[9];
  
  G[0*3+0] = 0.276805493111691;
  G[0*3+1] = 0.213024867595501;
  G[0*3+2] = 0.209205876527443;
  G[1*3+0] = 0.213024867595501;
  G[1*3+1] = 0.419642083051144;
  G[1*3+2] = 0.286592141665043;
  G[2*3+0] = 0.209205876527443;
  G[2*3+1] = 0.286592141665044;
  G[2*3+2] = 0.844559256324453;
  
  B[0*3+0] = -2.36794416971567;
  B[0*3+1] =  1.08817493992579;
  B[0*3+2] =  1.22601259339234;
  B[1*3+0] =  1.08817493992579;
  B[1*3+1] = -2.72352378723346;
  B[1*3+2] =  1.51348094527252;
  B[2*3+0] =  1.22601259339234;
  B[2*3+1] =  1.51348094527252;
  B[2*3+2] = -2.98729895217208;

  /* domain bounds */
  double xmin[4], xmax[4];
  xmin[0] = -1.0;
  xmin[1] = -1.0;
  xmin[2] = -0.1;
  xmin[3] = -0.1;
  xmax[0] = 1.0;
  xmax[1] = 1.0;
  xmax[2] = 0.1;
  xmax[3] = 0.1;

  /* location of initial Dirac - default values */
  double x0[4];
  x0[0]   = 0.159629725683297;
  x0[1]   = 0.281854642778853;
  x0[2]   = 0.0;
  x0[3]   = 0.0;

  /* reading physical model specific inputs - all processes */
  printf("Reading physical model inputs from file \"physics.inp\".\n");
  in = fopen("physics.inp","r");
  if (!in) {
    fprintf(stderr,"Error: File \"physics.inp\" not found.\n");
    return(1);
  } else {
    char word[500];
    ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
    if (!strcmp(word, "begin")){
	    while (strcmp(word, "end")){
		    ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
        if      (!strcmp(word,"Pm1_avg"   ))  {ferr=fscanf(in,"%lf",&Pm1_avg   ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"Pm2_avg"   ))  {ferr=fscanf(in,"%lf",&Pm2_avg   ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"Pmref_avg" ))  {ferr=fscanf(in,"%lf",&Pmref_avg ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"H1"        ))  {ferr=fscanf(in,"%lf",&H1        ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"H2"        ))  {ferr=fscanf(in,"%lf",&H2        ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"Href"      ))  {ferr=fscanf(in,"%lf",&Href      ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"E1"        ))  {ferr=fscanf(in,"%lf",&E1        ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"E2"        ))  {ferr=fscanf(in,"%lf",&E2        ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"Eref"      ))  {ferr=fscanf(in,"%lf",&Eref      ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"omegaB"    ))  {ferr=fscanf(in,"%lf",&omegaB    ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"gamma"     ))  {ferr=fscanf(in,"%lf",&gamma     ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"sigma"  ))  {
          ferr=fscanf(in,"%lf",&sigma[0][0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&sigma[0][1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&sigma[1][0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&sigma[1][1]) ;if(ferr!=1)return(1);
        } else if (!strcmp(word,"lambda"))  {
          ferr=fscanf(in,"%lf",&lambda[0][0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&lambda[0][1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&lambda[1][0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&lambda[1][1]) ;if(ferr!=1)return(1);
        } else if (!strcmp(word,"G"))  {
          ferr=fscanf(in,"%lf",&G[0*3+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&G[0*3+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&G[0*3+2]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&G[1*3+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&G[1*3+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&G[1*3+2]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&G[2*3+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&G[2*3+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&G[2*3+2]) ;if(ferr!=1)return(1);
        } else if (!strcmp(word,"B"))  {
          ferr=fscanf(in,"%lf",&B[0*3+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&B[0*3+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&B[0*3+2]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&B[1*3+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&B[1*3+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&B[1*3+2]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&B[2*3+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&B[2*3+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&B[2*3+2]) ;if(ferr!=1)return(1);
        } else if (!strcmp(word,"xmin"))  {
          ferr=fscanf(in,"%lf",&xmin[0]); if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&xmin[1]); if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&xmin[2]); if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&xmin[3]); if(ferr!=1)return(1);
        } else if (!strcmp(word,"xmax"))  {
          ferr=fscanf(in,"%lf",&xmax[0]); if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&xmax[1]); if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&xmax[2]); if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&xmax[3]); if(ferr!=1)return(1);
        } else if (!strcmp(word,"x0"))  {
          ferr=fscanf(in,"%lf",&x0[0]); if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&x0[1]); if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&x0[2]); if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&x0[3]); if(ferr!=1)return(1);
        }
      }
	  } else {
    	fprintf(stderr,"Error: Illegal format in file \"physics.inp\".\n");
      return(1);
	  }
    fclose(in);
  }

  /* compute drift velocity at specified equilibrium */
  double theta1 = x0[0];
  double theta2 = x0[1];
  double Omega1 = x0[2];
  double Omega2 = x0[3];
  printf("Specified equilibrium location: %1.16E, %1.16E, %1.16E, %1.16E\n",
         theta1, theta2, Omega1, Omega2);

  double Pe1, Pe2, Peref;
  ComputeElectricalPower(theta1,theta2,E1,E2,Eref,G,B,&Pe1,&Pe2,&Peref);

  double F1 = Pm1_avg / (2*H1) - Pmref_avg / (2*Href);
  double F2 = Pm2_avg / (2*H2) - Pmref_avg / (2*Href);
  double S1 = Pe1 / (2*H1) - Peref / (2*Href);
  double S2 = Pe2 / (2*H2) - Peref / (2*Href);

  double drift[4];
  drift[0] = omegaB * Omega1;
  drift[1] = omegaB * Omega2;
  drift[2] = F1 - gamma*Omega1 - S1;
  drift[3] = F2 - gamma*Omega2 - S2;
  printf("Computed drift velocity at equilibrium: %1.16E, %1.16E, %1.16E, %1.16E.\n",
         drift[0],drift[1],drift[2],drift[3]);

  /* allocate and generate the grid */
  printf("Generating grid.\n");
	double dx[4];
  dx[0] = (xmax[0]-xmin[0]) / ((double)(N[0]));
  dx[1] = (xmax[1]-xmin[1]) / ((double)(N[1]));
  dx[2] = (xmax[2]-xmin[2]) / ((double)(N[2]));
  dx[3] = (xmax[3]-xmin[3]) / ((double)(N[3]));
	double *x;
  int sizex = N[0] + N[1] + N[2] + N[3];
	x = (double*) calloc (sizex, sizeof(double));
  offset = 0;
  for (n = 0; n < 4; n++) {
    for (i = 0; i < N[n]; i++) x[i+offset] = xmin[n] + i*dx[n];
    offset += N[n];
  }

  /* find the grid point nearest to this initial Dirac position */
  printf("Finding grid point nearest to initial Dirac location.\n");
  int i0[4] = {-1,-1,-1,-1};
  offset = 0;
  for (n=0; n<4; n++) {
    double mindist = (x0[n]-x[offset])*(x0[n]-x[offset]);
    for (i=0; i<N[n]; i++) {
      double dist = (x0[n]-x[i+offset])*(x0[n]-x[i+offset]);
      if (dist <= mindist) {
        mindist = dist;
        i0[n] = i;
      }
    }
    offset += N[n];
  }
  printf("Grid point nearest to initial Dirac location: %d,%d,%d,%d.\n",
          i0[0],i0[1],i0[2],i0[3]);

  /* adjust location of initial Dirac to coincide with this grid point */
  offset = 0;
  for (n=0; n<4; n++) { x0[n] = x[i0[n]+offset]; offset += N[n]; }

  /* adjust Pm1_avg and Pm2_avg so that drift is zero at this adjusted Dirac location */
  printf("Adjusting Pm1_avg and Pm2_avg so that equilibrium point is at this grid location.\n");
  theta1 = x0[0];
  theta2 = x0[1];
  Omega1 = x0[2];
  Omega2 = x0[3];

  ComputeElectricalPower(theta1,theta2,E1,E2,Eref,G,B,&Pe1,&Pe2,&Peref);
  Pm1_avg = (2*H1) * (Pmref_avg / (2*Href) + Pe1 / (2*H1) - Peref / (2*Href));
  Pm2_avg = (2*H2) * (Pmref_avg / (2*Href) + Pe2 / (2*H2) - Peref / (2*Href));
  F1 = Pm1_avg / (2*H1) - Pmref_avg / (2*Href);
  F2 = Pm2_avg / (2*H2) - Pmref_avg / (2*Href);
  S1 = Pe1 / (2*H1) - Peref / (2*Href);
  S2 = Pe2 / (2*H2) - Peref / (2*Href);
  drift[0] = omegaB * Omega1;
  drift[1] = omegaB * Omega2;
  drift[2] = F1 - gamma*Omega1 - S1;
  drift[3] = F2 - gamma*Omega2 - S2;
  printf("Computed drift velocity at equilibrium: %1.16E, %1.16E, %1.16E, %1.16E.\n",
         drift[0],drift[1],drift[2],drift[3]);

  /* rewrite physics.inp with these adjusted values */
  printf("Rewriting physics.inp.\n");
  out = fopen("physics.inp","w");
  fprintf(out,"begin\n");
  fprintf(out,"\tPm1_avg            %1.16E\n",Pm1_avg);
  fprintf(out,"\tPm2_avg            %1.16E\n",Pm2_avg);
  fprintf(out,"\tPmref_avg          %1.16E\n",Pmref_avg);
  fprintf(out,"\tH1                 %1.16E\n",H1);       
  fprintf(out,"\tH2                 %1.16E\n",H2);       
  fprintf(out,"\tHref               %1.16E\n",Href);       
  fprintf(out,"\tE1                 %1.16E\n",E1);       
  fprintf(out,"\tE2                 %1.16E\n",E2);       
  fprintf(out,"\tEref               %1.16E\n",Eref);       
  fprintf(out,"\tomegaB             %1.16E\n",omegaB);
  fprintf(out,"\tgamma              %1.16E\n",gamma);
  fprintf(out,"\n");
  fprintf(out,"\tG\n");
  fprintf(out,"\t                   %1.16E %1.16E %1.16E\n",G[0*3+0],G[0*3+1],G[0*3+2]);
  fprintf(out,"\t                   %1.16E %1.16E %1.16E\n",G[1*3+0],G[1*3+1],G[1*3+2]);
  fprintf(out,"\t                   %1.16E %1.16E %1.16E\n",G[2*3+0],G[2*3+1],G[2*3+2]);
  fprintf(out,"\n");
  fprintf(out,"\tB\n");
  fprintf(out,"\t                   %1.16E %1.16E %1.16E\n",B[0*3+0],B[0*3+1],B[0*3+2]);
  fprintf(out,"\t                   %1.16E %1.16E %1.16E\n",B[1*3+0],B[1*3+1],B[1*3+2]);
  fprintf(out,"\t                   %1.16E %1.16E %1.16E\n",B[2*3+0],B[2*3+1],B[2*3+2]);
  fprintf(out,"\n");
  fprintf(out,"\txmin               %1.16E %1.16E %1.16E %1.16E\n",xmin[0],xmin[1],xmin[2],xmin[3]);
  fprintf(out,"\txmax               %1.16E %1.16E %1.16E %1.16E\n",xmax[0],xmax[1],xmax[2],xmax[3]);
  fprintf(out,"\tx0                 %1.16E %1.16E %1.16E %1.16E\n",x0  [0],x0  [1],x0  [2],x0  [3]);
  fprintf(out,"\n");
  fprintf(out,"\tlambda\n");
  fprintf(out,"\t                   %1.16E %1.16E\n",lambda[0][0], lambda[0][1]);
  fprintf(out,"\t                   %1.16E %1.16E\n",lambda[1][0], lambda[1][1]);
  fprintf(out,"\tsigma\n");
  fprintf(out,"\t                   %1.16E %1.16E\n",sigma[0][0], sigma[0][1]);
  fprintf(out,"\t                   %1.16E %1.16E\n",sigma[1][0], sigma[1][1]);
  fprintf(out,"end\n");
  fclose(out);

  /* Generate the solution and place the Dirac */
  printf("Generating solution.\n");
  int sizeu = N[0] * N[1] * N[2] * N[3];
	double *u;
	u = (double*) calloc (sizeu, sizeof(double));
  for (i = 0; i < sizeu; i++) u[i] = 0.0;
  /* Placing Dirac */
  printf("Placing Dirac at grid cell (%d, %d, %d, %d).\n",i0[0],i0[1],i0[2],i0[3]);
  double dv = dx[0] * dx[1] * dx[2] * dx[3];
  int p = i0[0] + N[0] * (i0[1] + N[1] * (i0[2] + N[2] * (i0[3])));
  u[p] = 1.0/dv;

  /* writing to file */
  printf("Writing to file.\n");
	out = fopen("initial.inp","wb");
  fwrite(x,sizeof(double),sizex,out);
  fwrite(u,sizeof(double),sizeu,out);
	fclose(out);

	free(x);
	free(u);

	return(0);
}
