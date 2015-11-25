/* 
  This code stitches the solution written out in parallel mode into
  a single file.

  The input files are the solution files written out in parallel mode:
  op.bin.xxxx

  For unsteady output, the code will write out binary files containing 
  global solution at each output time step (one file for each time step)
  op_xxxxx.bin

  For steady output, the code will write out one binary file containing
  the global solution:
  op.bin

*/

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

#define _MAX_STRING_SIZE_ 50

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
  FILE    *in;
  int     ndims, nvars, size, bytes,i,N_IORanks, proc, GroupSize;
  int     *dim_global,*dim_local,*iproc,IORank; proc;
  char    output_mode[_MAX_STRING_SIZE_], op_overwrite[_MAX_STRING_SIZE_];
  double  *Xg, *Ug;

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
            return(0);
          } else for (i=0; i<ndims; i++) fscanf(in,"%d",&dim_global[i]);
        } else if (!strcmp(word, "iproc")) {
          int i;
          if (!iproc) {
            fprintf(stderr,"Error in ReadInputs(): iproc not allocated.\n");
            fprintf(stderr,"Please specify ndims before iproc.\n"         );
            return(0);
          } else for (i=0; i<ndims; i++) fscanf(in,"%d",&iproc[i]);
        } else if (!strcmp(word, "output_mode")) {
          fscanf(in,"%s",output_mode);
          if (strcmp(output_mode,"serial")) fscanf(in,"%d",&N_IORanks);
        } else if (!strcmp(word, "op_overwrite")) {
          fscanf(in,"%s",op_overwrite);
        }
      }
    } else {
  	  fprintf(stderr,"Error: Illegal format in file \"solver.inp\".\n");
      return(0);
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
    printf("\tSolution output  mode                      : %s\n",output_mode  );
    printf("\tNumber of IO ranks                         : %d\n",N_IORanks   );
  }

  if (strcmp(output_mode,"parallel")) {
    printf("Error: output_mode is not \"parallel\". Why are you using this code?\n");
    return(0);
  }

  FILE *inps[N_IORanks];

  /* open the files for each IO group */
  printf("Opening op.bin.xxxx files.\n");
  for (IORank = 0; IORank < N_IORanks; IORank++) {
    char filename[_MAX_STRING_SIZE_];
    MPIGetFilename("op.bin",IORank,filename);
    inps[IORank] = fopen(filename,"rb");
    if (!inps[IORank]) {
      printf("Error: Could not open %s for reading.\n",filename);
      return(0);
    }    
  }

  int nproc = 1;
  for (i=0; i<ndims; i++) nproc *= iproc[i];
  if (nproc%N_IORanks != 0) {
    printf("Error: nproc is not a multiple of N_IORanks. HyPar could not have run in this mode. Something wrong/fishy!\n");
    return(0);
  }

  if (!strcmp(op_overwrite,"no")) {
    /* for unsteady solution output, the output file for each IO group
       contains the solutions for all output time steps */
    char out_filename[_MAX_STRING_SIZE_];
    strcpy(out_filename,"op_00000.bin");
    int flag = 1;
    while (flag) {
      /* allocate global grid and solution arrays */
      size = 0;
      for (i=0; i<ndims; i++) size += dim_global[i];
      Xg = (double*) calloc (size, sizeof(double));
      size = nvars;
      for (i=0; i<ndims; i++) size *= dim_global[i];
      Ug = (double*) calloc (size, sizeof(double));

      for (IORank=0; IORank < N_IORanks; IORank++) {
        /* for each IO group, calculate its range of processes */
        GroupSize = nproc / N_IORanks;
        int Start = IORank      * GroupSize;
        int End   = (IORank+1)  * GroupSize;
        /* for each process in this IO group, read the solution and
           put it in the global arrays */
        for (proc=Start; proc < End; proc++) {
          /* calculate local dimensions */
          int ip[ndims],is[ndims],ie[ndims];
          MPIRanknD(ndims,proc,iproc,ip);
          MPILocalDomainLimits(ndims,proc,iproc,dim_global,is,ie);
          for (i=0; i<ndims; i++) dim_local[i] = ie[i]-is[i];

          /* allocate local arrays for grid and solution */
          double *Xl, *Ul;
          size = 0; for (i=0; i<ndims; i++) size += dim_local[i];
          Xl = (double*) calloc (size, sizeof(double));
          size = nvars; for (i=0; i<ndims; i++) size *= dim_local[i];
          Ul = (double*) calloc (size, sizeof(double));

          /* read the local solution */
          size = 0; for (i=0; i<ndims; i++) size += dim_local[i];
          bytes = fread(Xl,sizeof(double),size,inps[IORank]);
          if (bytes != size) printf("Error: Unable to read grid.\n");
          size = nvars; for (i=0; i<ndims; i++) size *= dim_local[i];
          bytes = fread(Ul,sizeof(double),size,inps[IORank]);
          if (bytes != size) printf("Error: Unable to read solution.\n");

          /* copy to global grid and solution */
          int offsetl=0, offsetg=0;
          for (i=0; i<ndims; i++) {
            int p; for (p=0; p<dim_local[i]; p++) Xg[p+is[i]+offsetg] = Xl[p+offsetl];
            offsetl += dim_local[i];
            offsetg += dim_global[i];
          }
          int done = 0; int index[ndims]; for(i=0; i<ndims; i++) index[i]=0;
          while (!done) {
            int p1; _ArrayIndex1DWO_(ndims,dim_global,index,is,0,p1);
            int p2; _ArrayIndex1D_  (ndims,dim_local ,index,   0,p2);
            int v; for (v=0; v<nvars; v++) Ug[nvars*p1+v] = Ul[nvars*p2+v];
            _ArrayIncrementIndex_(ndims,dim_local,index,done);
          }

          /* free local grid and solution arrays */
          free(Xl);
          free(Ul);
        }
      }
     
      /* write global grid and solution to file */
      printf("Writing solution file %s.\n",out_filename);
      FILE *out;
      out = fopen(out_filename,"wb");
      fwrite(&ndims,sizeof(int),1,out);
      fwrite(&nvars,sizeof(int),1,out);
      fwrite(dim_global,sizeof(int),ndims,out);
      size = 0;
      for (i=0; i<ndims; i++) size += dim_global[i];
      bytes = fwrite(Xg,sizeof(double),size,out);
      if (bytes != size) {
        printf("Error: unable to write grid to %s.\n",out_filename);
        return(0);
      }
      size = nvars;
      for (i=0; i<ndims; i++) size *= dim_global[i];
      bytes = fwrite(Ug,sizeof(double),size,out);
      if (bytes != size) {
        printf("Error: unable to write solution to %s.\n",out_filename);
        return(0);
      }
      fclose(out);

      /* free global grid and solution arrays and increment filename */
      free(Xg);
      free(Ug);
      IncrementFilename(out_filename);

      /* check if the local solution files have reached their end */
      for (IORank=0; IORank<N_IORanks; IORank++) {
        int check = fgetc(inps[IORank]);
        if (check == EOF) flag = 0;
        else fseek(inps[IORank],-sizeof(char),SEEK_CUR);
      }
    }
    
  } else {
    /* for steady solution */
    char out_filename[_MAX_STRING_SIZE_];
    strcpy(out_filename,"op.bin");
    /* allocate global grid and solution arrays */
    size = 0;
    for (i=0; i<ndims; i++) size += dim_global[i];
    Xg = (double*) calloc (size, sizeof(double));
    size = nvars;
    for (i=0; i<ndims; i++) size *= dim_global[i];
    Ug = (double*) calloc (size, sizeof(double));

    for (IORank=0; IORank < N_IORanks; IORank++) {
      /* for each IO group, calculate its range of processes */
      GroupSize = nproc / N_IORanks;
      int Start = IORank      * GroupSize;
      int End   = (IORank+1)  * GroupSize;
      /* for each process in this IO group, read the solution and
         put it in the global arrays */
      for (proc=Start; proc < End; proc++) {
        /* calculate local dimensions */
        int ip[ndims],is[ndims],ie[ndims];
        MPIRanknD(ndims,proc,iproc,ip);
        MPILocalDomainLimits(ndims,proc,iproc,dim_global,is,ie);
        for (i=0; i<ndims; i++) dim_local[i] = ie[i]-is[i];

        /* allocate local arrays for grid and solution */
        double *Xl, *Ul;
        size = 0; for (i=0; i<ndims; i++) size += dim_local[i];
        Xl = (double*) calloc (size, sizeof(double));
        size = nvars; for (i=0; i<ndims; i++) size *= dim_local[i];
        Ul = (double*) calloc (size, sizeof(double));

        /* read the local solution */
        size = 0; for (i=0; i<ndims; i++) size += dim_local[i];
        bytes = fread(Xl,sizeof(double),size,inps[IORank]);
        if (bytes != size) printf("Error: Unable to read grid.\n");
        size = nvars; for (i=0; i<ndims; i++) size *= dim_local[i];
        bytes = fread(Ul,sizeof(double),size,inps[IORank]);
        if (bytes != size) printf("Error: Unable to read solution.\n");

        /* copy to global grid and solution */
        int offsetl=0, offsetg=0;
        for (i=0; i<ndims; i++) {
          int p; for (p=0; p<dim_local[i]; p++) Xg[p+is[i]+offsetg] = Xl[p+offsetl];
          offsetl += dim_local[i];
          offsetg += dim_global[i];
        }
        int done = 0; int index[ndims]; for(i=0; i<ndims; i++) index[i]=0;
        while (!done) {
          int p1; _ArrayIndex1DWO_(ndims,dim_global,index,is,0,p1);
          int p2; _ArrayIndex1D_  (ndims,dim_local ,index,   0,p2);
          int v; for (v=0; v<nvars; v++) Ug[nvars*p1+v] = Ul[nvars*p2+v];
          _ArrayIncrementIndex_(ndims,dim_local,index,done);
        }

        /* free local grid and solution arrays */
        free(Xl);
        free(Ul);
      }
    }
     
    /* write global grid and solution to file */
    printf("Writing solution file %s.\n",out_filename);
    FILE *out;
    out = fopen(out_filename,"wb");
    fwrite(&ndims,sizeof(int),1,out);
    fwrite(&nvars,sizeof(int),1,out);
    fwrite(dim_global,sizeof(int),ndims,out);
    size = 0;
    for (i=0; i<ndims; i++) size += dim_global[i];
    bytes = fwrite(Xg,sizeof(double),size,out);
    if (bytes != size) {
      printf("Error: unable to write grid to %s.\n",out_filename);
      return(0);
    }
    size = nvars;
    for (i=0; i<ndims; i++) size *= dim_global[i];
    bytes = fwrite(Ug,sizeof(double),size,out);
    if (bytes != size) {
      printf("Error: unable to write solution to %s.\n",out_filename);
      return(0);
    }
    fclose(out);

    /* free global grid and solution arrays and increment filename */
    free(Xg);
    free(Ug);
    
  }

  /* close the files for each IO group */
  for (IORank = 0; IORank < N_IORanks; IORank++) fclose(inps[IORank]);

  free(dim_local);
  free(dim_global);
  free(iproc);
  return(0);
}
