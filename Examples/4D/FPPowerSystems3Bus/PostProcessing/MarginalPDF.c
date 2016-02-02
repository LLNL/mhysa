/*

  This code extracts the marginal PDFs from the 4D 
  PDF output and writes them to text files.
  (Note: output file must be binary.)

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

int WriteText(int N,double *x,double *u,char *f)
{
  FILE *out;
  out = fopen(f,"w");
  if (!out) {
    fprintf(stderr,"Error: could not open %s for writing.\n",f);
    return(1);
  }
  int i;
  for (i=0; i<N; i++) fprintf(out, "%4d\t%1.16E\t%1.16E\n",i,x[i],u[i]);
  fclose(out);
  return(0);
}

int MarginalPDF(char *filename,int overwrite)
{
  FILE *in;
  in = fopen(filename,"rb");
  if (!in) return(1);
  
  printf("Processing file %s.\n",filename);
  int ndims, nvars, dims[4];
  double *U,*x;
  
  /* read the file headers */
  fread(&ndims,sizeof(int),1,in);
  fread(&nvars,sizeof(int),1,in);
  
  /* some checks */
  if (ndims != 4) {
    printf("Error: ndims in %s not equal to 2!\n",filename);
    return(-1);
  }
  if (nvars != 1) {
    printf("Error: nvars in %s not equal to 1!\n",filename);
    return(-1);
  }
  
  /* read dimensions */
  fread(dims,sizeof(int),ndims,in);
  printf("Dimensions: %d x %d x %d x %d\n",dims[0],dims[1],dims[2],dims[3]);

  /* allocate grid and solution arrays */
  int sizex = dims[0] + dims[1] + dims[2] + dims[3];
  int sizeu = dims[0] * dims[1] * dims[2] * dims[3];
  x = (double*) calloc (sizex,sizeof(double));
  U = (double*) calloc (sizeu,sizeof(double));
  
  /* read grid and solution */
  fread(x,sizeof(double),sizex,in);
  fread(U,sizeof(double),sizeu,in);
  /* done reading */
  fclose(in);

  /* set output filenames */
  char txtfilet1[500], txtfileo1[500], txtfilet2[500], txtfileo2[500];
  strcpy(txtfilet1,filename); txtfilet1[0]='t'; txtfilet1[1]='1';
  strcpy(txtfileo1,filename); txtfileo1[0]='o'; txtfileo1[1]='1';
  strcpy(txtfilet2,filename); txtfilet2[0]='t'; txtfilet2[1]='2';
  strcpy(txtfileo2,filename); txtfileo2[0]='o'; txtfileo2[1]='2';
  if (overwrite) {
    txtfilet1[9]  = txtfilet2[9]  = 'd';
    txtfilet1[10] = txtfilet2[10] = 'a';
    txtfilet1[11] = txtfilet2[11] = 't';
    txtfileo1[9]  = txtfileo2[9]  = 'd';
    txtfileo1[10] = txtfileo2[10] = 'a';
    txtfileo1[11] = txtfileo2[11] = 't';
  } else {
    txtfilet1[3] = txtfilet2[3] = 'd';
    txtfilet1[4] = txtfilet2[4] = 'a';
    txtfilet1[5] = txtfilet2[5] = 't';
    txtfileo1[3] = txtfileo2[3] = 'd';
    txtfileo1[4] = txtfileo2[4] = 'a';
    txtfileo1[5] = txtfileo2[5] = 't';
  }
  
  /* variables for the 1D solution */
  double *U1d, integral;
  int    i[4], j;
  double *T1 = x;
  double *T2 = T1 + dims[0];
  double *O1 = T2 + dims[1];
  double *O2 = O1 + dims[2];
  
  /* calculate integral of the 4D PDF */
  integral = 0.0;
  for (i[0]=0; i[0]<dims[0]; i[0]++) {
    for (i[1]=0; i[1]<dims[1]; i[1]++) {
      for (i[2]=0; i[2]<dims[2]; i[2]++) {
        for (i[3]=0; i[3]<dims[3]; i[3]++) {

          int q = i[0] + dims[0] * (i[1] + dims[1] * (i[2] + dims[2] * i[3]));

          double dx[4];
          if      (i[0] == 0)          dx[0] = (T1[i[0]+1]-T1[i[0]]  );
          else if (i[0] == dims[0]-1)  dx[0] = (T1[i[0]]  -T1[i[0]-1]);
          else                         dx[0] = (T1[i[0]+1]-T1[i[0]-1])/2.0;
          if      (i[1] == 0)          dx[1] = (T2[i[1]+1]-T2[i[1]]  );
          else if (i[1] == dims[1]-1)  dx[1] = (T2[i[1]]  -T2[i[1]-1]);
          else                         dx[1] = (T2[i[1]+1]-T2[i[1]-1])/2.0;
          if      (i[2] == 0)          dx[2] = (O1[i[2]+1]-O1[i[2]]  );
          else if (i[2] == dims[2]-1)  dx[2] = (O1[i[2]]  -O1[i[2]-1]);
          else                         dx[2] = (O1[i[2]+1]-O1[i[2]-1])/2.0;
          if      (i[3] == 0)          dx[3] = (O2[i[3]+1]-O2[i[3]]  );
          else if (i[3] == dims[3]-1)  dx[3] = (O2[i[3]]  -O2[i[3]-1]);
          else                         dx[3] = (O2[i[3]+1]-O2[i[3]-1])/2.0;

          integral += U[q] * dx[0]* dx[1] * dx[2] * dx[3];
        }
      }
    }
  }
  printf("\tIntegral of 4D PDF: %1.16E\n",integral);
  
  /* Theta1 */
  printf("\tExtracting marginal PDF for theta1.\n");
  U1d = (double*) calloc (dims[0],sizeof(double));
  
  for (i[0]=0; i[0]<dims[0]; i[0]++) {

    U1d[i[0]] = 0.0;

    for (i[1]=0; i[1]<dims[1]; i[1]++) {
      for (i[2]=0; i[2]<dims[2]; i[2]++) {
        for (i[3]=0; i[3]<dims[3]; i[3]++) {

          int q = i[0] + dims[0] * (i[1] + dims[1] * (i[2] + dims[2] * i[3]));

          double dx[4];
          if      (i[0] == 0)          dx[0] = (T1[i[0]+1]-T1[i[0]]  );
          else if (i[0] == dims[0]-1)  dx[0] = (T1[i[0]]  -T1[i[0]-1]);
          else                         dx[0] = (T1[i[0]+1]-T1[i[0]-1])/2.0;
          if      (i[1] == 0)          dx[1] = (T2[i[1]+1]-T2[i[1]]  );
          else if (i[1] == dims[1]-1)  dx[1] = (T2[i[1]]  -T2[i[1]-1]);
          else                         dx[1] = (T2[i[1]+1]-T2[i[1]-1])/2.0;
          if      (i[2] == 0)          dx[2] = (O1[i[2]+1]-O1[i[2]]  );
          else if (i[2] == dims[2]-1)  dx[2] = (O1[i[2]]  -O1[i[2]-1]);
          else                         dx[2] = (O1[i[2]+1]-O1[i[2]-1])/2.0;
          if      (i[3] == 0)          dx[3] = (O2[i[3]+1]-O2[i[3]]  );
          else if (i[3] == dims[3]-1)  dx[3] = (O2[i[3]]  -O2[i[3]-1]);
          else                         dx[3] = (O2[i[3]+1]-O2[i[3]-1])/2.0;

          U1d[i[0]] += U[q] * dx[1] * dx[2] * dx[3];
        }
      }
    }
  }
  
  /* calculating integral of marginal PDF */
  integral = 0;
  for (j=0; j<dims[0]; j++) {
    double dx;
    if      (j == 0)         dx = (T1[j+1]-T1[j]  );
    else if (j == dims[0]-1) dx = (T1[j]  -T1[j-1]);
    else                     dx = (T1[j+1]-T1[j-1])/2.0;
    integral += U1d[j] * dx;
  }
  printf("\tIntegral of marginal PDF for theta1: %1.16E\n",integral);
  
  /* write to file */
  printf("\tWriting file %s.\n",txtfilet1);
  WriteText(dims[0],T1,U1d,txtfilet1);

  /* free arrays */
  free(U1d);

  /* Theta2 */
  printf("\tExtracting marginal PDF for theta2.\n");
  U1d = (double*) calloc (dims[1],sizeof(double));
  
  for (i[1]=0; i[1]<dims[1]; i[1]++) {

    U1d[i[1]] = 0.0;

    for (i[0]=0; i[0]<dims[0]; i[0]++) {
      for (i[2]=0; i[2]<dims[2]; i[2]++) {
        for (i[3]=0; i[3]<dims[3]; i[3]++) {

          int q = i[0] + dims[0] * (i[1] + dims[1] * (i[2] + dims[2] * i[3]));

          double dx[4];
          if      (i[0] == 0)          dx[0] = (T1[i[0]+1]-T1[i[0]]  );
          else if (i[0] == dims[0]-1)  dx[0] = (T1[i[0]]  -T1[i[0]-1]);
          else                         dx[0] = (T1[i[0]+1]-T1[i[0]-1])/2.0;
          if      (i[1] == 0)          dx[1] = (T2[i[1]+1]-T2[i[1]]  );
          else if (i[1] == dims[1]-1)  dx[1] = (T2[i[1]]  -T2[i[1]-1]);
          else                         dx[1] = (T2[i[1]+1]-T2[i[1]-1])/2.0;
          if      (i[2] == 0)          dx[2] = (O1[i[2]+1]-O1[i[2]]  );
          else if (i[2] == dims[2]-1)  dx[2] = (O1[i[2]]  -O1[i[2]-1]);
          else                         dx[2] = (O1[i[2]+1]-O1[i[2]-1])/2.0;
          if      (i[3] == 0)          dx[3] = (O2[i[3]+1]-O2[i[3]]  );
          else if (i[3] == dims[3]-1)  dx[3] = (O2[i[3]]  -O2[i[3]-1]);
          else                         dx[3] = (O2[i[3]+1]-O2[i[3]-1])/2.0;

          U1d[i[1]] += U[q] * dx[0] * dx[2] * dx[3];
        }
      }
    }
  }
  
  /* calculating integral of marginal PDF */
  integral = 0;
  for (j=0; j<dims[1]; j++) {
    double dx;
    if      (j == 0)         dx = (T2[j+1]-T2[j]  );
    else if (j == dims[1]-1) dx = (T2[j]  -T2[j-1]);
    else                     dx = (T2[j+1]-T2[j-1])/2.0;
    integral += U1d[j] * dx;
  }
  printf("\tIntegral of marginal PDF for theta2: %1.16E\n",integral);
  
  /* write to file */
  printf("\tWriting file %s.\n",txtfilet2);
  WriteText(dims[1],T2,U1d,txtfilet2);

  /* free arrays */
  free(U1d);

  /* Omega1 */
  printf("\tExtracting marginal PDF for omega1.\n");
  U1d = (double*) calloc (dims[2],sizeof(double));
  
  for (i[2]=0; i[2]<dims[2]; i[2]++) {

    U1d[i[2]] = 0.0;

    for (i[1]=0; i[1]<dims[1]; i[1]++) {
      for (i[0]=0; i[0]<dims[0]; i[0]++) {
        for (i[3]=0; i[3]<dims[3]; i[3]++) {

          int q = i[0] + dims[0] * (i[1] + dims[1] * (i[2] + dims[2] * i[3]));

          double dx[4];
          if      (i[0] == 0)          dx[0] = (T1[i[0]+1]-T1[i[0]]  );
          else if (i[0] == dims[0]-1)  dx[0] = (T1[i[0]]  -T1[i[0]-1]);
          else                         dx[0] = (T1[i[0]+1]-T1[i[0]-1])/2.0;
          if      (i[1] == 0)          dx[1] = (T2[i[1]+1]-T2[i[1]]  );
          else if (i[1] == dims[1]-1)  dx[1] = (T2[i[1]]  -T2[i[1]-1]);
          else                         dx[1] = (T2[i[1]+1]-T2[i[1]-1])/2.0;
          if      (i[2] == 0)          dx[2] = (O1[i[2]+1]-O1[i[2]]  );
          else if (i[2] == dims[2]-1)  dx[2] = (O1[i[2]]  -O1[i[2]-1]);
          else                         dx[2] = (O1[i[2]+1]-O1[i[2]-1])/2.0;
          if      (i[3] == 0)          dx[3] = (O2[i[3]+1]-O2[i[3]]  );
          else if (i[3] == dims[3]-1)  dx[3] = (O2[i[3]]  -O2[i[3]-1]);
          else                         dx[3] = (O2[i[3]+1]-O2[i[3]-1])/2.0;

          U1d[i[2]] += U[q] * dx[1] * dx[0] * dx[3];
        }
      }
    }
  }
  
  /* calculating integral of marginal PDF */
  integral = 0;
  for (j=0; j<dims[2]; j++) {
    double dx;
    if      (j == 0)         dx = (O1[j+1]-O1[j]  );
    else if (j == dims[2]-1) dx = (O1[j]  -O1[j-1]);
    else                     dx = (O1[j+1]-O1[j-1])/2.0;
    integral += U1d[j] * dx;
  }
  printf("\tIntegral of marginal PDF for omega1: %1.16E\n",integral);
  
  /* write to file */
  printf("\tWriting file %s.\n",txtfileo1);
  WriteText(dims[2],O1,U1d,txtfileo1);

  /* free arrays */
  free(U1d);

  /* Omega2 */
  printf("\tExtracting marginal PDF for omega2.\n");
  U1d = (double*) calloc (dims[3],sizeof(double));
  
  for (i[3]=0; i[3]<dims[3]; i[3]++) {

    U1d[i[3]] = 0.0;

    for (i[1]=0; i[1]<dims[1]; i[1]++) {
      for (i[2]=0; i[2]<dims[2]; i[2]++) {
        for (i[0]=0; i[0]<dims[0]; i[0]++) {

          int q = i[0] + dims[0] * (i[1] + dims[1] * (i[2] + dims[2] * i[3]));

          double dx[4];
          if      (i[0] == 0)          dx[0] = (T1[i[0]+1]-T1[i[0]]  );
          else if (i[0] == dims[0]-1)  dx[0] = (T1[i[0]]  -T1[i[0]-1]);
          else                         dx[0] = (T1[i[0]+1]-T1[i[0]-1])/2.0;
          if      (i[1] == 0)          dx[1] = (T2[i[1]+1]-T2[i[1]]  );
          else if (i[1] == dims[1]-1)  dx[1] = (T2[i[1]]  -T2[i[1]-1]);
          else                         dx[1] = (T2[i[1]+1]-T2[i[1]-1])/2.0;
          if      (i[2] == 0)          dx[2] = (O1[i[2]+1]-O1[i[2]]  );
          else if (i[2] == dims[2]-1)  dx[2] = (O1[i[2]]  -O1[i[2]-1]);
          else                         dx[2] = (O1[i[2]+1]-O1[i[2]-1])/2.0;
          if      (i[3] == 0)          dx[3] = (O2[i[3]+1]-O2[i[3]]  );
          else if (i[3] == dims[3]-1)  dx[3] = (O2[i[3]]  -O2[i[3]-1]);
          else                         dx[3] = (O2[i[3]+1]-O2[i[3]-1])/2.0;

          U1d[i[3]] += U[q] * dx[1] * dx[2] * dx[0];
        }
      }
    }
  }
  
  /* calculating integral of marginal PDF */
  integral = 0;
  for (j=0; j<dims[3]; j++) {
    double dx;
    if      (j == 0)         dx = (O2[j+1]-O2[j]  );
    else if (j == dims[3]-1) dx = (O2[j]  -O2[j-1]);
    else                     dx = (O2[j+1]-O2[j-1])/2.0;
    integral += U1d[j] * dx;
  }
  printf("\tIntegral of marginal PDF for omega2: %1.16E\n",integral);
  
  /* write to file */
  printf("\tWriting file %s.\n",txtfileo2);
  WriteText(dims[3],O2,U1d,txtfileo2);

  /* free arrays */
  free(U1d);

  return(0);
}

int main()
{
  FILE *out1, *out2, *in, *inputs;
  char filename[50], op_file_format[50], op_overwrite[50];

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
   			else if (!strcmp(word, "op_overwrite"     ))  fscanf(inputs,"%s" ,op_overwrite  );
      }
    }
    fclose(inputs);
  }
  if (strcmp(op_file_format,"binary") && strcmp(op_file_format,"bin")) {
    printf("Error: solution output needs to be in binary files.\n");
    return(0);
  }

  int err;
  if (!strcmp(op_overwrite,"no")) {
    strcpy(filename,"op_00000.bin");
    while(1) {
      err = MarginalPDF(filename,1);
      if (err == 1) {
        printf("No more files found, exiting.\n");
        break;
      } else if (err == -1) {
        printf("Invalid format detected in file. Exiting.\n");
        break;
      }
      IncrementFilename(filename);
    }
  } else {
    strcpy(filename,"op.bin");
    err = MarginalPDF(filename,0);
    if (err == -1) printf("Invalid format detected in file.\n");
  }

  return(0);
}
