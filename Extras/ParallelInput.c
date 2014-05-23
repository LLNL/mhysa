/* This code generates the files for reading the initial 
 * (and exact, if available) solutions in parallel, from
 * the "initial.inp" and "exact.inp" files containing the 
 * initial and exact solutions of the complete domain
 *
 * The new files are "initial_par.inp.xxxx" and 
 * "exact_par.inp.xxxx" (where xxxx corresponds to
 * the IO rank that the file is meant for).
 *
 * Binary files are supported.
*/

#define _MAX_STRING_SIZE_ 50

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

#define _ArrayIndex1DWO_(N,imax,i,offset,ghost,index) \
  { \
    index = i[N-1]+(ghost)+ offset[N-1];\
    int arraycounter; \
    for (arraycounter = (N)-2; arraycounter > -1; arraycounter--) { \
      index = ((index*(imax[arraycounter]+2*(ghost))) + (i[arraycounter]+(ghost)+offset[arraycounter]));\
    } \
  }

void SplitDomain(char*, char*);

int main()
{
  char filename_in [_MAX_STRING_SIZE_];
  char filename_out[_MAX_STRING_SIZE_];

  strcpy(filename_in ,"initial.inp"    );
  strcpy(filename_out,"initial_par.inp");
  SplitDomain(filename_in,filename_out);

  strcpy(filename_in ,"exact.inp"    );
  strcpy(filename_out,"exact_par.inp");
  SplitDomain(filename_in,filename_out);

  return(0);
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

void SplitDomain(char *fnamein, char *fnameout)
{
  FILE  *in;
  int   ndims, nvars, size, bytes,i,N_IORanks;
  int   *dim_global,*dim_local,*iproc;
  char  ip_file_type[_MAX_STRING_SIZE_], input_mode[_MAX_STRING_SIZE_];
  double *Xg, *Ug;

  in = fopen("solver.inp","r");
  if (!in) {
    fprintf(stderr,"Error: File \"solver.inp\" not found.\n");
    return;
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

  /* checks */
  if (strcmp(ip_file_type,"bin") && strcmp(ip_file_type,"binary")) {
    printf("Error: this script is for binary files only.\n");
    return;
  }
  if (strcmp(input_mode,"parallel")) {
    printf("Error: input_mode is not \"parallel\".\n");
    return;
  }

  /* Read the solution for the entire domain */
  in = fopen(fnamein,"rb");
  if (!in) {
    printf("File %s not found.\n",fnamein);
    free(dim_local);
    free(dim_global);
    free(iproc);
    return;
  }
  size = 0;
  for (i=0; i<ndims; i++) size += dim_global[i];
  Xg = (double*) calloc (size, sizeof(double));
  bytes = fread(Xg,sizeof(double),size,in);
  if (bytes != size) {
    printf("Error: unable to read grid.\n");
    return;
  }
  size = nvars;
  for (i=0; i<ndims; i++) size *= dim_global[i];
  Ug = (double*) calloc (size, sizeof(double));
  bytes = fread(Ug,sizeof(double),size,in);
  if (bytes != size) {
    printf("Error: unable to read solution.\n");
    return;
  }
  fclose(in);

  int nproc = 1;
  for (i=0; i<ndims; i++) nproc *= iproc[i];
  if (nproc%N_IORanks != 0) N_IORanks = 1;
  printf("Splitting data into %d processes. Will generate %d files (one for each file IO rank).\n",nproc,N_IORanks);

  int proc,IORank;
  FILE *out;
  int GroupSize = nproc / N_IORanks;
  for (IORank = 0; IORank < N_IORanks; IORank++) {
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

      size = nvars; for (i=0; i<ndims; i++) size *= dim_local[i];
      Ul = (double*) calloc (size, sizeof(double));
      int done = 0; int index[ndims]; for(i=0; i<ndims; i++) index[i]=0;
      while (!done) {
        int p1; _ArrayIndex1DWO_(ndims,dim_global,index,is,0,p1);
        int p2; _ArrayIndex1D_  (ndims,dim_local ,index,   0,p2);
        int v; for (v=0; v<nvars; v++) Ul[nvars*p2+v] = Ug[nvars*p1+v];
        _ArrayIncrementIndex_(ndims,dim_local,index,done);
      }

      size = 0; for (i=0; i<ndims; i++) size += dim_local[i];
      bytes = fwrite(Xl,sizeof(double),size,out);
      if (bytes != size) printf("Error: Unable to write data to file %s.\n",fnameout);
      size = nvars; for (i=0; i<ndims; i++) size *= dim_local[i];
      bytes = fwrite(Ul,sizeof(double),size,out);
      if (bytes != size) printf("Error: Unable to write data to file %s.\n",fnameout);

      free(Xl);
      free(Ul);
    }
    fclose(out);

  }
  
  free(Xg);
  free(Ug);
  free(dim_local);
  free(dim_global);
  free(iproc);
  return;
}
