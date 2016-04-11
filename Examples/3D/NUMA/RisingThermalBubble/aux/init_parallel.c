/*

  This code generates the initial solution for the 
  3D rising thermal bubble case for the non-hydrostatic
  unified model of the atmosphere (NUMA).

  The initial solution is written for the local (MPI) 
  sub-domains, and the number of files written 
  depend on the number of I/O ranks. Although this
  code is serial, it will compute the decomposition
  of the global domain and generate the initial
  solution accordingly.

*/

#define _MAX_STRING_SIZE_ 50

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define _ArrayIncrementIndex_(N,imax,i,done) \
  { \
    int arraycounter = 0; \
    while (arraycounter < (N)) { \
      if (i[arraycounter] == imax[arraycounter]-1) { \
        i[arraycounter] = 0; \
        arraycounter++; \
      } else { \
        i[arraycounter]++; \
        break; \
      } \
    } \
    if (arraycounter == (N)) done = 1; \
    else          done = 0; \
  }

#define _ArrayIndex1D_(N,imax,i,ghost,index)  \
  { \
    index = i[N-1]+(ghost); \
    int arraycounter; \
    for (arraycounter = (N)-2; arraycounter > -1; arraycounter--) { \
      index = ((index*(imax[arraycounter]+2*(ghost))) + (i[arraycounter]+(ghost))); \
    } \
  }

double raiseto(double x, double a)
{
  return(exp(a*log(x)));
}

void GetStringFromInteger(int a,char *A,int width)
{
  int i;
  for (i=0; i<width; i++) {
    char digit = (char) (a%10 + '0'); 
    a /= 10;
    A[width-1-i] = digit;
  }
  return;
}

void MPIGetFilename(char *root,int rank,char *filename)
{
  char  tail[_MAX_STRING_SIZE_] = "";

  GetStringFromInteger(rank,tail,4);
  strcpy(filename,"");
  strcat(filename,root);
  strcat(filename,"." );
  strcat(filename,tail);

  return;
}

int MPIRanknD(int ndims,int rank,int* iproc,int *ip)
{
  int i,term    = 1;
  for (i=0; i<ndims; i++) term *= iproc[i];
  for (i=ndims-1; i>=0; i--) {
    term /= iproc[i];
    ip[i] = rank/term;
    rank -= ip[i]*term;
  }
  return(0);
}

int MPIPartition1D(int nglobal,int nproc,int rank)
{
  int nlocal;
  if (nglobal%nproc == 0) nlocal = nglobal/nproc;
  else {
    if (rank == nproc-1)  nlocal = nglobal/nproc + nglobal%nproc;
    else                  nlocal = nglobal/nproc;
  }
  return(nlocal);
}

int MPILocalDomainLimits(int ndims,int p,int *iproc,int *dim_global,int *is, int *ie) 
{
  int i;
  int ip[ndims];
  MPIRanknD(ndims,p,iproc,ip);

  for (i=0; i<ndims; i++) {
    int imax_local, isize, root = 0;
    imax_local = MPIPartition1D(dim_global[i],iproc[i],root );
    isize      = MPIPartition1D(dim_global[i],iproc[i],ip[i]);
    if (is)  is[i] = ip[i]*imax_local;
    if (ie)  ie[i] = ip[i]*imax_local + isize;
  }
  return(0);
}

int main()
{
  FILE    *in, *out;
  int     NI, NJ, NK, ndims, nvars, size, bytes, N_IORanks, i, j, k;
  int     *dim_global,*dim_local,*iproc;
  char    ip_file_type[_MAX_STRING_SIZE_], input_mode[_MAX_STRING_SIZE_], fnameout[_MAX_STRING_SIZE_];
  strcpy  (ip_file_type,"ascii");
  strcpy  (fnameout,"initial_par.inp");

  printf("Reading file \"solver.inp\"...\n");
  in = fopen("solver.inp","r");
  if (!in) {
    fprintf(stderr,"Error: File \"solver.inp\" not found.\n");
    return(0);
  } else {
	  char word[_MAX_STRING_SIZE_];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")){
	    while (strcmp(word, "end")){
		    fscanf(in,"%s",word);
        if (!strcmp(word, "ndims")) {
          fscanf(in,"%d",&ndims);
          dim_global = (int*) calloc (ndims,sizeof(int));
          dim_local  = (int*) calloc (ndims,sizeof(int));
          iproc      = (int*) calloc (ndims,sizeof(int));
        }	else if (!strcmp(word, "nvars")) {
          fscanf(in,"%d",&nvars);
        } else if (!strcmp(word, "size")) {
          int i;
          if (!dim_global) {
            fprintf(stderr,"Error in ReadInputs(): dim_global not allocated.\n");
            fprintf(stderr,"Please specify ndims before dimensions.\n"         );
            return;
          } else {
            for (i=0; i<ndims; i++) fscanf(in,"%d",&dim_global[i]);
          }
        } else if (!strcmp(word, "iproc")) {
          int i;
          if (!iproc) {
            fprintf(stderr,"Error in ReadInputs(): iproc not allocated.\n");
            fprintf(stderr,"Please specify ndims before iproc.\n"         );
            return;
          } else {
            for (i=0; i<ndims; i++) fscanf(in,"%d",&iproc[i]);
          }
        } else if (!strcmp(word, "ip_file_type" )) {
          fscanf(in,"%s",ip_file_type);
        } else if (!strcmp(word, "input_mode")) {
          fscanf(in,"%s",input_mode);
          if (strcmp(input_mode,"serial")) fscanf(in,"%d",&N_IORanks);
        }
      }
    } else {
  	  fprintf(stderr,"Error: Illegal format in file \"solver.inp\".\n");
      return;
    }
    fclose(in);

    /* Print to screen the inputs read */
	  printf("\tNo. of dimensions                          : %d\n",ndims);
	  printf("\tNo. of variables                           : %d\n",nvars);
	  printf("\tDomain size                                : ");
    for (i=0; i<ndims; i++) printf ("%d ",dim_global[i]);
    printf("\n");
	  printf("\tProcesses along each dimension             : ");
    for (i=0; i<ndims; i++) printf ("%d ",iproc[i]);
    printf("\n");
    printf("\tInitial solution file type                 : %s\n",ip_file_type);
    printf("\tInitial solution read mode                 : %s\n",input_mode  );
    printf("\tNumber of IO ranks                         : %d\n",N_IORanks   );
  }

  if (ndims != 3) {
    printf("ndims is not 3 in solver.inp. this code is to generate 3D exact solution\n");
    return(0);
  }
  if (!strcmp(ip_file_type,"ascii")) {
    printf("Error: ip_file_type *must* be specified and set to \"binary\" in solver.inp.\n");
    return(0);
  }
  if (strcmp(input_mode,"parallel")) {
    printf("Error: input_mode is not \"parallel\".\n");
    return;
  }

  /* read in zero-altitude pressure and temperature */
  double P0 = 100000.0;
  double T0 = 300.0;
  double g     = 9.8;
  double R     = 287.058;
  double gamma = 1.4;
  printf("Reading file \"physics.inp\"...\n");
  in = fopen("physics.inp","r");
  if (!in) printf ("Error: Input file \"physics.inp\" not found. Default values will be used.\n");
  else {
    char word[500];
    fscanf(in,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(in,"%s",word);
        if      (!strcmp(word, "Pref" )) fscanf(in,"%lf",&P0);
        else if (!strcmp(word, "Tref" )) fscanf(in,"%lf",&T0);
        else if (!strcmp(word, "gamma")) fscanf(in,"%lf",&gamma);
        else if (!strcmp(word, "g"    )) fscanf(in,"%lf",&g);
        else if (!strcmp(word, "R"    )) fscanf(in,"%lf",&R);
      }
    } else printf("Error: Illegal format in physics.inp. Crash and burn!\n");
  }
  fclose(in);
  double inv_gamma_m1 = 1.0 / (gamma-1.0);
  double Cp = gamma * inv_gamma_m1 * R;

  /* Define the domain */
  double xmin, xmax, ymin, ymax, zmin, zmax;
  xmin    =  0.0;
  xmax    =  1000.0;
  ymin    =  0.0;
  ymax    =  1000.0;
  zmin    =  0.0;
  zmax    =  1000.0 ;

  NI = dim_global[0];
  NJ = dim_global[1];
  NK = dim_global[2];
  double Lx = xmax - xmin;
  double Ly = ymax - ymin;
  double Lz = zmax - zmin;
	double dx = Lx / ((double)NI-1);
	double dy = Ly / ((double)NJ-1);
	double dz = Lz / ((double)NK-1);

  /* Initial perturbation center */
  double xc = 500;
  double yc = 500;
  double zc = 260;

  /* initial perturbation parameters */
  double pi = 4.0*atan(1.0);
  double rc = 250.0;

  printf("Generating grid.\n");
  double *Xg = (double*) calloc (NI+NJ+NK, sizeof(double));
  double *x = Xg, *y = Xg+NI, *z = Xg+NI+NJ;
	for (i = 0; i < NI; i++){
  	for (j = 0; j < NJ; j++){
  	  for (k = 0; k < NK; k++){
  	  	x[i] = xmin + i*dx;
	    	y[j] = ymin + j*dy;
	    	z[k] = zmin + k*dz;
      }
	  }
	}

  int nproc = 1;
  for (i=0; i<ndims; i++) nproc *= iproc[i];
  if (nproc%N_IORanks != 0) N_IORanks = 1;
  printf("Splitting data into %d processes. Will generate %d files (one for each file IO rank).\n",
         nproc,N_IORanks);

  int proc,IORank;
  int GroupSize = nproc / N_IORanks;
  for (IORank = 0; IORank < N_IORanks; IORank++) {
    printf("Generating and writing local solutions for IORank %d.\n",IORank);
    char out_filename[_MAX_STRING_SIZE_];
    MPIGetFilename(fnameout,IORank,out_filename);

    int Start = IORank      * GroupSize;
    int End   = (IORank+1)  * GroupSize;

    out = fopen(out_filename,"wb");
    for (proc=Start; proc < End; proc++) {

      int ip[ndims],is[ndims],ie[ndims];
      double *Xl, *Ul;
      MPIRanknD(ndims,proc,iproc,ip);
      MPILocalDomainLimits(ndims,proc,iproc,dim_global,is,ie);
      for (i=0; i<ndims; i++) dim_local[i] = ie[i]-is[i];

      size = 0; for (i=0; i<ndims; i++) size += dim_local[i];
      Xl = (double*) calloc (size, sizeof(double));
      int offsetl=0, offsetg=0;
      for (i=0; i<ndims; i++) {
        int p; for (p=0; p<dim_local[i]; p++) Xl[p+offsetl] = Xg[p+is[i]+offsetg];
        offsetl += dim_local[i];
        offsetg += dim_global[i];
      }
      x = Xl;
      y = Xl + dim_local[0];
      z = Xl + dim_local[0] + dim_local[1];

      size = nvars; for (i=0; i<ndims; i++) size *= dim_local[i];
      Ul = (double*) calloc (size, sizeof(double));
      int done = 0; int index[ndims]; for(i=0; i<ndims; i++) index[i]=0;
      while (!done) {
        int p; _ArrayIndex1D_(ndims,dim_local ,index,0,p);
        i = index[0]; j = index[1]; k = index[2];

        /* temperature peturbation */
        double r = sqrt((x[i]-xc)*(x[i]-xc)+(y[j]-yc)*(y[j]-yc)+(z[k]-zc)*(z[k]-zc));
        double dtheta = (r>rc ? 0.0 : (0.5*(1.0+cos(pi*r/rc))) );

        double theta      = T0 + dtheta; /* temperature potential */
        double Pexner     = 1.0 - (g/(Cp*T0))*z[k]; /* Exner pressure */
        double rho        = (P0/(R*theta)) * raiseto(Pexner,inv_gamma_m1); /* density */

        double theta_ref  = T0; /* reference temperature potential  */
        double Pexner_ref = 1.0 - (g/(Cp*T0))*z[k]; /* reference Exner pressure */
        double rho_ref    = (P0/(R*theta_ref)) * raiseto(Pexner_ref,inv_gamma_m1); /* reference density */

        Ul[5*p+0] = rho - rho_ref;
        Ul[5*p+1] = 0.0;
        Ul[5*p+2] = 0.0;
        Ul[5*p+3] = 0.0;
        Ul[5*p+4] = (rho*theta) - (rho_ref*theta_ref);
        _ArrayIncrementIndex_(ndims,dim_local,index,done);
      }

      size = 0; for (i=0; i<ndims; i++) size += dim_local[i];
      bytes = fwrite(Xl,sizeof(double),size,out);
      if (bytes != size) printf("Error: Unable to write grid data to file %s.\n",fnameout);
      size = nvars; for (i=0; i<ndims; i++) size *= dim_local[i];
      bytes = fwrite(Ul,sizeof(double),size,out);
      if (bytes != size) printf("Error: Unable to write solution data to file %s.\n",fnameout);

      free(Xl);
      free(Ul);
    }
    fclose(out);

  }

  free(dim_global);
  free(dim_local);
  free(iproc);
  free(Xg);
  return(0);
}
