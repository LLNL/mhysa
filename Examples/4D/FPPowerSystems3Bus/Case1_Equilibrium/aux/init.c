#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


int LUDecomp(double *A, double *rhs, int N)
{
	int i, j, k;
  for (i = 0; i < N; i++) {
    if (A[i*N+i] == 0){
      fprintf(stderr,"Error: Zero encountered on main diagonal!\n");
			return(1);
		}
		for (j = i+1; j < N; j++){
			double factor = A[j*N+i] / A[i*N+i];
			A[j*N+i] = 0;
			for (k = i+1; k < N; k++) A[j*N+k] -= factor * A[i*N+k];
			rhs[j] -= factor * rhs[i];
		}
	}
	for (i = N-1; i >=0; i--){
		double sum = 0;
		for (j = i+1; j < N; j++) sum += A[i*N+j] * rhs[j];
		rhs[i] = (rhs[i] - sum) / A[i*N+i];
	}
	return(0);
}

int MatInverse(double *A, double *B, int N)
{
	int i, j;
	double *r  = (double*) calloc(N,sizeof(double));
	double *AA = (double*) calloc(N*N,sizeof(double));
	for (i = 0; i < N; i++){
    for (j=0; j<N*N; j++) AA[j] = A[j];
		for (j = 0; j < N; j++){
			if (j == i)	r[j] = 1.0;
			else		    r[j] = 0.0;
		}
		int ierr = LUDecomp(AA,r,N); if (ierr) return(ierr);
		for (j = 0; j < N; j++)	B[j*N+i] = r[j];
	}
  free(r);
  free(AA);
	return(0);
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
  if (strcmp(ip_file_type,"binary") && strcmp(ip_file_type,"bin")) {
    printf("Error: Invalid ip_file_type: Binary only!\n");
    return(1);
  }

  /* check if ip_file_type is binary */
  if (strcmp(ip_file_type,"binary") && strcmp(ip_file_type,"bin")) {
    printf("Error: ip_file_type needs to be bin or binary.\n");
    return(0);
  }

  /* declarations */
  int M = 6;
  double PM1, H1, D1, E1, Xd1;
  double PM2, H2, D2, E2, Xd2;
  double omegaB,alpha,beta,G[M*M/4],B[M*M/4],lambda[2][2],sigma[2][2];

  /* default values of model parameters */
  PM1     = 3.109260511864138;
  PM2     = 1.0;
  H1      = 1000.640;
  H2      = 3.64;
  omegaB  = 2*PI*60;
  D1      = 450;
  D2      = 1.0;
  E1      = 1.044225672060674;
  E2      = 1.034543707656856;
  Xd1     = 0.02;
  Xd2     = 0.2;
  alpha   = 4.386890097147679;
  beta    = -1.096722524286920;

  lambda[0][0] = 0.1;
  lambda[0][1] = 0.0;
  lambda[1][0] = 0.0;
  lambda[1][1] = 0.1;

  sigma[0][0] = 1.0;
  sigma[0][1] = 0.0;
  sigma[1][0] = 0.0;
  sigma[1][1] = 1.0;

  G[0*M/2+0] =  7.631257631257632;
  G[0*M/2+1] = -3.815628815628816;
  G[0*M/2+2] = -3.815628815628816;
  G[1*M/2+0] = -3.815628815628816;
  G[1*M/2+1] =  6.839334669523348;
  G[1*M/2+2] = -3.023705853894533;
  G[2*M/2+0] = -3.815628815628816;
  G[2*M/2+1] = -3.023705853894533;
  G[2*M/2+2] =  6.839334669523348;

  B[0*M/2+0] = -38.053788156288157;
  B[0*M/2+1] =  19.078144078144078;
  B[0*M/2+2] =  19.078144078144078;
  B[1*M/2+0] =  19.078144078144078;
  B[1*M/2+1] = -34.081673347616743;
  B[1*M/2+2] =  15.118529269472663;
  B[2*M/2+0] =  19.078144078144078;
  B[2*M/2+1] =  15.118529269472663;
  B[2*M/2+2] = -34.081673347616743;

  /* domain bounds */
  double xmin[4], xmax[4];
  xmin[0] = 0.0;
  xmin[1] = 0.0;
  xmin[2] = 0.95;
  xmin[3] = 0.95;
  xmax[0] = 1.0;
  xmax[1] = 1.0;
  xmax[2] = 1.05;
  xmax[3] = 1.05;

  /* location of initial Dirac - default values */
  double x0[4];
  x0[0]   = 0.058417047872067;
  x0[1]   = 0.169798193956055;
  x0[2]   = 1.0;
  x0[3]   = 1.0;

  /* reading physical model specific inputs - all processes */
  printf("Reading physical model inputs from file \"physics.inp\".\n");
  in = fopen("physics.inp","r");
  if (!in) printf("Warning: File \"physics.inp\" not found. Using default values.\n");
  else {
    char word[100];
    ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
    if (!strcmp(word, "begin")){
	    while (strcmp(word, "end")){
		    ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
        if      (!strcmp(word,"PM1"    ))  {ferr=fscanf(in,"%lf",&PM1   ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"H1"     ))  {ferr=fscanf(in,"%lf",&H1    ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"D1"     ))  {ferr=fscanf(in,"%lf",&D1    ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"E1"     ))  {ferr=fscanf(in,"%lf",&E1    ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"Xd1"    ))  {ferr=fscanf(in,"%lf",&Xd1   ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"PM2"    ))  {ferr=fscanf(in,"%lf",&PM2   ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"H2"     ))  {ferr=fscanf(in,"%lf",&H2    ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"D2"     ))  {ferr=fscanf(in,"%lf",&D2    ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"E2"     ))  {ferr=fscanf(in,"%lf",&E2    ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"Xd2"    ))  {ferr=fscanf(in,"%lf",&Xd2   ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"omegaB" ))  {ferr=fscanf(in,"%lf",&omegaB) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"alpha"  ))  {ferr=fscanf(in,"%lf",&alpha ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"beta"   ))  {ferr=fscanf(in,"%lf",&beta  ) ;if(ferr!=1)return(1);}
        else if (!strcmp(word,"G"))  {
          ferr=fscanf(in,"%lf",&G[0*M/2+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&G[0*M/2+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&G[0*M/2+2]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&G[1*M/2+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&G[1*M/2+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&G[1*M/2+2]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&G[2*M/2+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&G[2*M/2+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&G[2*M/2+2]) ;if(ferr!=1)return(1);
        } else if (!strcmp(word,"B"))  {
          ferr=fscanf(in,"%lf",&B[0*M/2+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&B[0*M/2+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&B[0*M/2+2]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&B[1*M/2+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&B[1*M/2+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&B[1*M/2+2]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&B[2*M/2+0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&B[2*M/2+1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&B[2*M/2+2]) ;if(ferr!=1)return(1);
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
        } else if (!strcmp(word,"sigma"  ))  {
          ferr=fscanf(in,"%lf",&sigma[0][0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&sigma[0][1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&sigma[1][0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&sigma[1][1]) ;if(ferr!=1)return(1);
        } else if (!strcmp(word,"lambda"))  {
          ferr=fscanf(in,"%lf",&lambda[0][0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&lambda[0][1]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&lambda[1][0]) ;if(ferr!=1)return(1);
          ferr=fscanf(in,"%lf",&lambda[1][1]) ;if(ferr!=1)return(1);
        }
      }
	  } else {
    	fprintf(stderr,"Error: Illegal format in file \"physics.inp\".\n");
      return(1);
	  }
    fclose(in);
  }

  /* Some initial calculations of the physical parameters */
  double A[M*M],Ainv[M*M];

  A[0*M+0] =  G[0*M/2+0];
  A[0*M+1] = -B[0*M/2+0] + 1.0/Xd1;
  A[1*M+0] =  B[0*M/2+0] - 1.0/Xd1;
  A[1*M+1] =  G[0*M/2+0];

  A[0*M+2] =  G[0*M/2+1];
  A[0*M+3] = -B[0*M/2+1];
  A[1*M+2] =  B[0*M/2+1];
  A[1*M+3] =  G[0*M/2+1];

  A[0*M+4] =  G[0*M/2+2];
  A[0*M+5] = -B[0*M/2+2];
  A[1*M+4] =  B[0*M/2+2];
  A[1*M+5] =  G[0*M/2+2];

  A[2*M+0] =  G[1*M/2+0];
  A[2*M+1] = -B[1*M/2+0];
  A[3*M+0] =  B[1*M/2+0];
  A[3*M+1] =  G[1*M/2+0];

  A[2*M+2] =  G[1*M/2+1];
  A[2*M+3] = -B[1*M/2+1] + 1.0/Xd2;
  A[3*M+2] =  B[1*M/2+1] - 1.0/Xd2;
  A[3*M+3] =  G[1*M/2+1];

  A[2*M+4] =  G[1*M/2+2];
  A[2*M+5] = -B[1*M/2+2];
  A[3*M+4] =  B[1*M/2+2];
  A[3*M+5] =  G[1*M/2+2];

  A[4*M+0] =  G[2*M/2+0];
  A[4*M+1] = -B[2*M/2+0];
  A[5*M+0] =  B[2*M/2+0];
  A[5*M+1] =  G[2*M/2+0];

  A[4*M+2] =  G[2*M/2+1];
  A[4*M+3] = -B[2*M/2+1];
  A[5*M+2] =  B[2*M/2+1];
  A[5*M+3] =  G[2*M/2+1];

  A[4*M+4] =  G[2*M/2+2] + alpha;
  A[4*M+5] = -B[2*M/2+2] - beta;
  A[5*M+4] =  B[2*M/2+2] + beta;
  A[5*M+5] =  G[2*M/2+2] + alpha;
  
  int inverr = MatInverse(A,Ainv,M);
  if (inverr) {
    fprintf(stderr,"Error in FPPowerSystem3BusInitialize(): Unable to invert matrix!\n");
    return(0);
  }

#if 0
  /* Verifying Ainv is correct */
  double eye[M*M];
  int j,k;
  for (j=0; j<M; j++) {
    for (k=0; k<M; k++) {
      eye[j*M+k] = 0.0;
      for (i=0; i<M; i++) eye[j*M+k] += (A[j*M+i] * Ainv[i*M+k]);
    }
  }
  printf("\n");
  printf("Verifying Ainv is correct. A*Ainv = \n");
  for (j=0; j<M; j++) {
    for (k=0; k<M; k++) printf("%+1.16E  ",eye[j*M+k]);
    printf("\n");
  }
  printf("\n");
#endif

  /* calculate location of equilibrium */
  double theta1, theta2, omega1, omega2;
  theta1 = x0[0];
  theta2 = x0[1];
  omega1 = x0[2];
  omega2 = x0[3];
  printf("Initial Dirac location: %1.16E, %1.16E, %1.16E, %1.16E\n",
          theta1, theta2, omega1, omega2);

  double VR1, VI1, VR2, VI2;
  VR1 = Ainv[0*M+0]*E1*sin(theta1)/Xd1 - Ainv[0*M+1]*E1*cos(theta1)/Xd1 
      + Ainv[0*M+2]*E2*sin(theta2)/Xd2 - Ainv[0*M+3]*E2*cos(theta2)/Xd2;
  VI1 = Ainv[1*M+0]*E1*sin(theta1)/Xd1 - Ainv[1*M+1]*E1*cos(theta1)/Xd1 
      + Ainv[1*M+2]*E2*sin(theta2)/Xd2 - Ainv[1*M+3]*E2*cos(theta2)/Xd2;
  VR2 = Ainv[2*M+0]*E1*sin(theta1)/Xd1 - Ainv[2*M+1]*E1*cos(theta1)/Xd1 
      + Ainv[2*M+2]*E2*sin(theta2)/Xd2 - Ainv[2*M+3]*E2*cos(theta2)/Xd2;
  VI2 = Ainv[3*M+0]*E1*sin(theta1)/Xd1 - Ainv[3*M+1]*E1*cos(theta1)/Xd1 
      + Ainv[3*M+2]*E2*sin(theta2)/Xd2 - Ainv[3*M+3]*E2*cos(theta2)/Xd2;


  /* verify drift is zero for these values */
  double zero1 = PM1 - VR1*E1*sin(theta1)/Xd1 + VI1*E1*cos(theta1)/Xd1;
  double zero2 = PM2 - VR2*E2*sin(theta2)/Xd2 + VI2*E2*cos(theta2)/Xd2;
  printf("Computed drift velocity at equilibrium: %1.16E, %1.16E.\n",zero1,zero2);

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

  printf("Adjusting PM1 and PM2 so that equilibrium point is at this grid location.\n");
  theta1 = x0[0];
  theta2 = x0[1];
  omega1 = x0[2];
  omega2 = x0[3];

  /* recompute the voltages */
  VR1 = Ainv[0*M+0]*E1*sin(theta1)/Xd1 - Ainv[0*M+1]*E1*cos(theta1)/Xd1 
      + Ainv[0*M+2]*E2*sin(theta2)/Xd2 - Ainv[0*M+3]*E2*cos(theta2)/Xd2;
  VI1 = Ainv[1*M+0]*E1*sin(theta1)/Xd1 - Ainv[1*M+1]*E1*cos(theta1)/Xd1 
      + Ainv[1*M+2]*E2*sin(theta2)/Xd2 - Ainv[1*M+3]*E2*cos(theta2)/Xd2;
  VR2 = Ainv[2*M+0]*E1*sin(theta1)/Xd1 - Ainv[2*M+1]*E1*cos(theta1)/Xd1 
      + Ainv[2*M+2]*E2*sin(theta2)/Xd2 - Ainv[2*M+3]*E2*cos(theta2)/Xd2;
  VI2 = Ainv[3*M+0]*E1*sin(theta1)/Xd1 - Ainv[3*M+1]*E1*cos(theta1)/Xd1 
      + Ainv[3*M+2]*E2*sin(theta2)/Xd2 - Ainv[3*M+3]*E2*cos(theta2)/Xd2;

  /* adjust PM1 and PM2 so that drift is zero at this adjusted Dirac location */
  PM1 = VR1*E1*sin(theta1)/Xd1 - VI1*E1*cos(theta1)/Xd1;
  PM2 = VR2*E2*sin(theta2)/Xd2 - VI2*E2*cos(theta2)/Xd2;
  printf("Adjusted Dirac location: %1.16E, %1.16E, %1.16E, %1.16E\n",
          theta1, theta2, omega1, omega2);
  /* verify drift is zero for these values */
  zero1 = PM1 - VR1*E1*sin(theta1)/Xd1 + VI1*E1*cos(theta1)/Xd1;
  zero2 = PM2 - VR2*E2*sin(theta2)/Xd2 + VI2*E2*cos(theta2)/Xd2;
  printf("Recomputed drift velocity at equilibrium: %1.16E, %1.16E.\n",zero1,zero2);

  /* rewrite physics.inp with these adjusted values */
  printf("Rewriting physics.inp.\n");
  out = fopen("physics.inp","w");
  fprintf(out,"begin\n");
  fprintf(out,"\tPM1                %1.16E\n",PM1);
  fprintf(out,"\tH1                 %1.16E\n",H1);
  fprintf(out,"\tD1                 %1.16E\n",D1);
  fprintf(out,"\tE1                 %1.16E\n",E1);
  fprintf(out,"\tXd1                %1.16E\n",Xd1);
  fprintf(out,"\n");
  fprintf(out,"\tPM2                %1.16E\n",PM2);
  fprintf(out,"\tH2                 %1.16E\n",H2);
  fprintf(out,"\tD2                 %1.16E\n",D2);
  fprintf(out,"\tE2                 %1.16E\n",E2);
  fprintf(out,"\tXd2                %1.16E\n",Xd2);
  fprintf(out,"\n");
  fprintf(out,"\tomegaB             %1.16E\n",omegaB);
  fprintf(out,"\talpha              %1.16E\n",alpha );
  fprintf(out,"\tbeta               %1.16E\n",beta  );
  fprintf(out,"\n");
  fprintf(out,"\tG\n");
  fprintf(out,"\t                   %1.16E %1.16E %1.16E\n",G[0*M/2+0],G[0*M/2+1],G[0*M/2+2]);
  fprintf(out,"\t                   %1.16E %1.16E %1.16E\n",G[1*M/2+0],G[1*M/2+1],G[1*M/2+2]);
  fprintf(out,"\t                   %1.16E %1.16E %1.16E\n",G[2*M/2+0],G[2*M/2+1],G[2*M/2+2]);
  fprintf(out,"\n");
  fprintf(out,"\tB\n");
  fprintf(out,"\t                   %1.16E %1.16E %1.16E\n",B[0*M/2+0],B[0*M/2+1],B[0*M/2+2]);
  fprintf(out,"\t                   %1.16E %1.16E %1.16E\n",B[1*M/2+0],B[1*M/2+1],B[1*M/2+2]);
  fprintf(out,"\t                   %1.16E %1.16E %1.16E\n",B[2*M/2+0],B[2*M/2+1],B[2*M/2+2]);
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
