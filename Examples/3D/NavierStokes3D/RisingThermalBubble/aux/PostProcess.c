/*
  Rising Thermal Bubble:-
  The code takes a binary solution file (that contains the 
  conserved variable (rho, rho*u, rho*v, rho*w, e) as its input 
  and calculates the primitive atmospheric flow variables:
  rho, u, v, w, P, theta, pi, rho0, P0, theta0, pi0
  and writes them to a tecplot or text file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct _parameters_{
  double grav_x, grav_y, grav_z, R, gamma, P_ref, rho_ref;
  int HB;
} Parameters;

#define _ArraySetValue_(x,size,value)                                                                               \
  {                                                                                                                 \
    int arraycounter;                                                                                               \
    for (arraycounter = 0; arraycounter < (size); arraycounter++)  x[arraycounter] = (value);                       \
  }

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

void IncrementFilename(char *f)
{
  if (f[7] == '9') {
    f[7] = '0';
    if (f[6] == '9') {
      f[6] = '0';
      if (f[5] == '9') {
        f[5] = '0';
        if (f[4] == '9') {
          f[4] = '0';
          if (f[3] == '9') {
            f[3] = '0';
            fprintf(stderr,"Warning: file increment hit max limit. Resetting to zero.\n");
          } else {
            f[3]++;
          }
        } else {
          f[4]++;
        }
      } else {
        f[5]++;
      }
    } else {
      f[6]++;
    }
  } else {
    f[7]++;
  }
}

double raiseto(double x, double a)
{
  return(exp(a*log(x)));
}

int WriteTecplot3D(int ndims,int nvars,int *dim,double *x,double *u,char *f,int *index)
{
  if (ndims !=3) {
    fprintf(stderr,"Error in WriteTecplot3D(): This functions is hardcoded for 3-dimensional ");
    fprintf(stderr,"problems only. Instead, ndims=%d.\n",ndims);
    return(1);
  }
  int i;
  int imax = dim[0];
  int jmax = dim[1];
  int kmax = dim[2];

  printf("Writing tecplot solution file %s.\n",f);
  FILE *out;
  out = fopen(f,"w");
  if (!out) {
    fprintf(stderr,"Error: could not open %s for writing.\n",f);
    return(1);
  }

  /* writing tecplot data file headers */
  fprintf(out,"VARIABLES=\"I\",\"J\",\"K\",\"X\",\"Y\",\"Z\",");
  char varname[3] = "00";
  for (i = 0; i < nvars; i++) {
    fprintf(out,"\"%s\",",varname);
    if (varname[1] == '9') { varname[0]++; varname[1] = '0'; }
    else                     varname[1]++;
  }
  fprintf(out,"\n");
  fprintf(out,"ZONE I=%d,J=%d,K=%d,F=POINT\n",imax,jmax,kmax);

  /* writing the data */
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int i, p;
    _ArrayIndex1D_(ndims,dim,index,0,p);
    for (i=0; i<ndims; i++) fprintf(out,"%4d ",index[i]);
    for (i=0; i<ndims; i++) {
      int j,offset = 0; for (j=0; j<i; j++) offset += dim[j];
      fprintf(out,"%+E ",x[offset+index[i]]);
    }
    for (i=0; i<nvars; i++) fprintf(out,"%+E ",u[nvars*p+i]);
    fprintf(out,"\n");
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  fclose(out);
  return(0);
}

int WriteText(int ndims,int nvars,int *dim,double *x,double *u,char *f,int *index)
{
  FILE *out;
  out = fopen(f,"w");
  if (!out) {
    fprintf(stderr,"Error: could not open %s for writing.\n",f);
    return(1);
  }

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int i, p;
    _ArrayIndex1D_(ndims,dim,index,0,p);
    for (i=0; i<ndims; i++) fprintf(out,"%4d ",index[i]);
    for (i=0; i<ndims; i++) {
      int j,offset = 0; for (j=0; j<i; j++) offset += dim[j];
      fprintf(out,"%+E ",x[offset+index[i]]);
    }
    for (i=0; i<nvars; i++) fprintf(out,"%+E ",u[nvars*p+i]);
    fprintf(out,"\n");
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  fclose(out);
  return(0);
}

int PostProcess(char *fname, char *oname, void *p, int flag)
{
  Parameters *params = (Parameters*) p;
  FILE *in; in = fopen(fname,"rb");

  if (!in) return(-1);

  printf("Reading file %s.\n",fname);
  int ndims, nvars;
  double *U,*x;

  /* read the file headers */
  fread(&ndims,sizeof(int),1,in);
  fread(&nvars,sizeof(int),1,in);

  /* some checks */
  if (ndims != 3) {
    printf("Error: ndims in %s not equal to 3!\n",fname);
    return(1);
  }
  if (nvars != 5) {
    printf("Error: nvars in %s not equal to 5!\n",fname);
    return(1);
  }

  /* read dimensions */
  int dims[ndims];
  fread(dims,sizeof(int),ndims,in);
  printf("Dimensions: %d x %d x %d\n",dims[0],dims[1],dims[2]);
  printf("Nvars     : %d\n",nvars);

  /* allocate grid and solution arrays */
  x = (double*) calloc (dims[0]+dims[1]+dims[2]       ,sizeof(double));
  U = (double*) calloc (dims[0]*dims[1]*dims[2]*nvars ,sizeof(double));

  /* read grid and solution */
  fread(x,sizeof(double),dims[0]+dims[1]+dims[2]      ,in);
  fread(U,sizeof(double),dims[0]*dims[1]*dims[2]*nvars,in);
  /* done reading */
  fclose(in);

  int imax = dims[0];
  int jmax = dims[1];
  int kmax = dims[2];

  /* allocate primitive variable array (rho, u, v, w, P, theta, rho0, P0, pi0, theta0) */
  int evars = 5;
  double *Q = (double*) calloc ((nvars+evars)*imax*jmax*kmax,sizeof(double));
    
  /* calculate primitive variables */
  int i, j, k;
  double *X           = x;
  double *Y           = x+imax;
  double *Z           = x+imax+jmax;
  double grav_y       = params->grav_y;
  double R            = params->R;
  double gamma        = params->gamma;
  double P_ref        = params->P_ref;
  double rho_ref      = params->rho_ref;
  double T_ref        = P_ref / (R*rho_ref);
  double inv_gamma_m1 = 1.0 / (gamma-1.0);
  double Cp           = gamma * inv_gamma_m1 * R;

  for (i=0; i<imax; i++) {
    for (j=0; j<jmax; j++) {
      for (k=0; k<kmax; k++) {
        int p = i + imax*j + imax*jmax*k;

        double rho0, theta0, Pexner, P0;
        theta0  = T_ref;
        Pexner  = 1.0 - (grav_y*Y[j])/(Cp*T_ref);
        rho0    = (P_ref/(R*theta0)) * raiseto(Pexner, inv_gamma_m1);
        P0      = P_ref   * raiseto(Pexner, gamma*inv_gamma_m1);

        double rho, uvel, vvel, wvel, E, P, theta;
        rho   = U[nvars*p+0];
        uvel  = U[nvars*p+1] / rho;
        vvel  = U[nvars*p+2] / rho;
        wvel  = U[nvars*p+3] / rho;
        E     = U[nvars*p+4];
        P     = (gamma-1.0) * (E - 0.5*rho*(uvel*uvel+vvel*vvel+wvel*wvel));
        theta = (E-0.5*rho*(uvel*uvel+vvel*vvel+wvel*wvel))/(Pexner*rho) * ((gamma-1.0)/R);

        Q[(nvars+evars)*p+0] = rho;
        Q[(nvars+evars)*p+1] = uvel;
        Q[(nvars+evars)*p+2] = vvel;
        Q[(nvars+evars)*p+3] = wvel;
        Q[(nvars+evars)*p+4] = P;
        Q[(nvars+evars)*p+5] = theta;
        Q[(nvars+evars)*p+6] = rho0;
        Q[(nvars+evars)*p+7] = P0;
        Q[(nvars+evars)*p+8] = Pexner;
        Q[(nvars+evars)*p+9] = theta0;
      }
    }
  }

  int dim[ndims], index[ndims];
  dim[0] = imax;
  dim[1] = jmax;
  dim[2] = kmax;

  /* write Tecplot/Text file */
  if (flag) WriteTecplot3D(ndims,nvars+evars,&dim[0],x,Q,oname,&index[0]);
  else      WriteText     (ndims,nvars+evars,&dim[0],x,Q,oname,&index[0]);

  /* clean up */
  free(U);
  free(Q);
  free(x);
}

int main()
{
  FILE *out1, *out2, *in, *inputs;
  char filename[50], op_file_format[50], tecfile[50], overwrite[50];

  int flag;
  printf("Write tecplot file (1) or plain text file (0): ");
  scanf("%d",&flag);

  if ((flag != 1) && (flag != 0)) {
    printf("Error: Invalid input. Should be 1 or 0.\n");
    return(0);
  }

  printf("Reading solver.inp.\n");
  inputs = fopen("solver.inp","r");
  if (!inputs) {
    fprintf(stderr,"Error: File \"solver.inp\" not found.\n");
    return(1);
  } else {
	  char word[100];
    fscanf(inputs,"%s",word);
    if (!strcmp(word, "begin")){
	    while (strcmp(word, "end")){
		    fscanf(inputs,"%s",word);
   			if      (!strcmp(word, "op_file_format"   ))  fscanf(inputs,"%s" ,op_file_format);
   			else if (!strcmp(word, "op_overwrite"     ))  fscanf(inputs,"%s" ,overwrite      );
      }
    }
    fclose(inputs);
  }
  if (strcmp(op_file_format,"binary") && strcmp(op_file_format,"bin")) {
    printf("Error: solution output needs to be in binary files.\n");
    return(0);
  }

  Parameters params;
  /* default values */
  params.grav_x   = 0.0;
  params.grav_y   = 9.8;
  params.grav_z   = 0.0;
  params.R        = 287.058;
  params.gamma    = 1.4;
  params.P_ref    = 100000.0;
  params.rho_ref  = 100000.0 / (params.R * 300.0);
  params.HB       = 0;
  /* read these parameters from file */
  printf("Reading physics.inp.\n");
  inputs = fopen("physics.inp","r");
  if (!inputs) {
    fprintf(stderr,"Error: File \"physics.inp\" not found.\n");
    return(1);
  } else {
	  char word[100];
    fscanf(inputs,"%s",word);
    if (!strcmp(word, "begin")){
	    while (strcmp(word, "end")){
		    fscanf(inputs,"%s",word);
   			if      (!strcmp(word, "gamma"))    fscanf(inputs,"%lf",&params.gamma);
   			else if (!strcmp(word, "gravity")) {
          fscanf(inputs,"%lf",&params.grav_x);
          fscanf(inputs,"%lf",&params.grav_y);
          fscanf(inputs,"%lf",&params.grav_z);
        } else if (!strcmp(word,"p_ref"))   fscanf(inputs,"%lf",&params.P_ref);
        else if   (!strcmp(word,"rho_ref")) fscanf(inputs,"%lf",&params.rho_ref);
        else if   (!strcmp(word,"HB"))      fscanf(inputs,"%d",&params.HB);
      }
    }
    fclose(inputs);
  }

  if (!strcmp(overwrite,"no")) {
    strcpy(filename,"op_00000.bin");
    while(1) {
      /* set filename */
      strcpy(tecfile,filename);
      tecfile[9]  = 'd';
      tecfile[10] = 'a';
      tecfile[11] = 't';
      int err = PostProcess(filename, tecfile, &params, flag);
      if (err == -1) {
        printf("No more files found. Exiting.\n");
        break;
      }
      IncrementFilename(filename);
    }
  } else if (!strcmp(overwrite,"yes")) {
    strcpy(filename,"op.bin");
    /* set filename */
    strcpy(tecfile,filename);
    tecfile[3] = 'd';
    tecfile[4] = 'a';
    tecfile[5] = 't';
    int err = PostProcess(filename, tecfile, &params, flag);
    if (err == -1) {
      printf("Error: op.bin not found.\n");
      return(0);
    }
  }

  return(0);
}
