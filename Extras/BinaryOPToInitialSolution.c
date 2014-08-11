/*
  This code converts a binary solution file to file 
  solution.inp in the same format as a binary initial.inp.
  To use is as an initial solution (for example, in a 
  restart run), rename it to initial.inp and place it in
  the run directory.
  If input mode is parallel, use ParallelInput.c
  thereafter.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main()
{
  FILE *out, *in, *inputs;
  char filename[50], op_file_format[50];

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
      }
    }
    fclose(inputs);
  }
  if (strcmp(op_file_format,"binary") && strcmp(op_file_format,"bin")) {
    printf("Error: solution output needs to be in binary files.\n");
    return(0);
  }

  printf("Enter filename (op_xxxxx.bin): ");
  scanf("%s",filename);
  
  in = fopen(filename,"rb");

  if (!in) {
    printf("Error: File %s not found.\n",filename);
    return(1);
  } else {
    printf("Reading file %s.\n",filename);

    /* read the file headers */
    int ndims, nvars, d;
    fread(&ndims,sizeof(int),1,in);
    fread(&nvars,sizeof(int),1,in);
    int dims[ndims];
    /* read dimensions */
    fread(dims,sizeof(int),ndims,in);

    /* allocate grid and solution arrays */
    int sizex = 0;     for (d=0; d<ndims; d++) sizex += dims[d];
    int sizeu = nvars; for (d=0; d<ndims; d++) sizeu *= dims[d];
    double *U,*x;
    x = (double*) calloc (sizex,sizeof(double));
    U = (double*) calloc (sizeu,sizeof(double));

    printf("Size of x vector: %d\nSize of U vector: %d\n",sizex,sizeu);

    /* read grid and solution */
    fread(x,sizeof(double),sizex,in);
    fread(U,sizeof(double),sizeu,in);
    /* done reading */
    fclose(in);

    /* write initial solution file */
    out = fopen("solution.inp","wb");
    fwrite(x,sizeof(double),sizex,out);
    fwrite(U,sizeof(double),sizeu,out);
    fclose(out);

    /* clean up */
    free(U);
    free(x);
  }

  return(0);
}
