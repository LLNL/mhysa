/*! @file PetscRegisterTIMethods.c
    @brief Register a custom time-integration method
    @author Debojyoti Ghosh
*/

#ifdef with_petsc

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <petscinterface.h>

#undef __FUNCT__
#define __FUNCT__ "PetscRegisterTIMethods"

/*!
  This function allows the registering of custom multi-stage time integration 
  methods and using it with PETSc. More than one method may be provided, each 
  with a unique name (the name \b must \b not conflict with the names of existing 
  methods in PETSc). The methods should be provided in a file \a "time_method.inp", 
  and the following must be provided:
  + Name
  + Class (i.e. RK or ARK)
  + Number of stages
  + Theoretical order
  + Butcher tableaux

  The method can then be used by invoking it by its name. For example, if a 
  custom ARKIMEX method \a "foo" is provided through \a "time_method.inp", it can be
  used by specifying \a "-ts_type arkimex -ts_arkimex_type foo" in the command
  line or \a ".petscrc" file.

  See \a /Examples/PETScInputs for examples of the input files required to 
  provide a time integration method.

  Currently supported classes of time integrators for which custom methods
  can be provided:
  + RK (Runge-Kutta): See TSRKRegister 
    (http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSRKRegister.html)
  + ARK (Additive Runge-Kutta): See TSARKIMEXRegister 
    (http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSRKRegister.html)

  To do:
  + Add support for TSGLEE when it gets merged to PETSc's master.
*/
int PetscRegisterTIMethods(
                            int rank /*!< MPI rank */
                          )
{
  PetscErrorCode ierr;
  int            ierr2;

  PetscFunctionBegin;

  /* Note: all processors read and register the custom methods */
  /* instead of root doing it and broadcasting.                */
  FILE *in;
  in = fopen("time_method.inp","r");
  if (in) {
    int n, N; /* N = total number of methods specified in the file */
    ierr2 = fscanf(in,"%d",&N); if (ierr2 != 1) return(1);
    for (n = 0; n < N; n++) {
      char      name[_MAX_STRING_SIZE_];      /* name of the scheme                           */
      char      type[_MAX_STRING_SIZE_];      /* type of scheme - ARKIMEX or RK               */
      PetscInt  s, order;                     /* number of stages and order                   */
      PetscReal *A , *b,  *c;                 /* Butcher tableaux entries for non-stiff terms */
      PetscReal *At, *bt, *ct;                /* Butcher tableaux entries for the stiff terms */
      PetscReal *bemb, *bembt;                /* Embedded method coefficients                 */
      PetscInt  pinterp;                      /* order of interpolation scheme                */
      PetscReal *bint,*bintt;                 /* Dense output interpolation coefficients      */

      /* Initializations */
      strcpy(name,"");
      s = order = pinterp = 0;
      A  = b  = c  = NULL;
      At = bt = ct = NULL;
      bemb = bembt = NULL;
      bint = bintt = NULL;

      /* Read the method */
      char word[_MAX_STRING_SIZE_];
      ierr2 = fscanf(in,"%s",word); if (ierr2 != 1) return(1);
      if (!strcmp(word,"begin")) {
        while (strcmp(word, "end")) {
          ierr2 = fscanf(in,"%s",word); if (ierr2 != 1) return(1);
          if      (!strcmp(word,"name"))    { ierr2 = fscanf(in,"%s",name);    if (ierr2 != 1) return(1); }
          else if (!strcmp(word,"class"))   { ierr2 = fscanf(in,"%s",type);    if (ierr2 != 1) return(1); }
          else if (!strcmp(word,"nstages")) { ierr2 = fscanf(in,"%d",&s);      if (ierr2 != 1) return(1); }
          else if (!strcmp(word,"order"))   { ierr2 = fscanf(in,"%d",&order);  if (ierr2 != 1) return(1); }
          else if (!strcmp(word,"pinterp")) { ierr2 = fscanf(in,"%d",&pinterp);if (ierr2 != 1) return(1); }
          else if (!strcmp(word,"At")) {
            if (s == 0) { 
              if (!rank) fprintf(stderr,"Error in PetscRegisterTIMethods(): nstages must be defined ");
              if (!rank) fprintf(stderr,"before specifying the Butcher tableaux entries.\n"              );
              return(1);
            } else {
              At = (PetscReal*) calloc (s*s, sizeof(PetscReal));
              int i, j;
              for (i = 0; i < s; i++) {
                for (j = 0; j < s; j++) {
                  ierr2 = fscanf(in,"%lf",&At[i*s+j]); if (ierr2 != 1) return(1);
                }
              }
            }
          } else if (!strcmp(word,"A")) {
            if (s == 0) { 
              if (!rank) fprintf(stderr,"Error in PetscRegisterTIMethods(): nstages must be defined ");
              if (!rank) fprintf(stderr,"before specifying the Butcher tableaux entries.\n"              );
              return(1);
            } else {
              A = (PetscReal*) calloc (s*s, sizeof(PetscReal));
              int i, j;
              for (i = 0; i < s; i++) {
                for (j = 0; j < s; j++) {
                  ierr2 = fscanf(in,"%lf",&A[i*s+j]); if (ierr2 != 1) return(1);
                }
              }
            }
          } else if (!strcmp(word,"bt")) {
            if (s == 0) { 
              if (!rank) fprintf(stderr,"Error in PetscRegisterTIMethods(): nstages must be defined ");
              if (!rank) fprintf(stderr,"before specifying the Butcher tableaux entries.\n"              );
              return(1);
            } else {
              bt = (PetscReal*) calloc (s, sizeof(PetscReal));
              int i;
              for (i = 0; i < s; i++) ierr2 = fscanf(in,"%lf",&bt[i]); if (ierr2 != 1) return(1);
            }
          } else if (!strcmp(word,"b")) {
            if (s == 0) { 
              if (!rank) fprintf(stderr,"Error in PetscRegisterTIMethods(): nstages must be defined ");
              if (!rank) fprintf(stderr,"before specifying the Butcher tableaux entries.\n"              );
              return(1);
            } else {
              b = (PetscReal*) calloc (s, sizeof(PetscReal));
              int i;
              for (i = 0; i < s; i++) ierr2 = fscanf(in,"%lf",&b[i]); if (ierr2 != 1) return(1);
            }
          } else if (!strcmp(word,"ct")) {
            if (s == 0) { 
              if (!rank) fprintf(stderr,"Error in PetscRegisterTIMethods(): nstages must be defined ");
              if (!rank) fprintf(stderr,"before specifying the Butcher tableaux entries.\n"              );
              return(1);
            } else {
              ct = (PetscReal*) calloc (s, sizeof(PetscReal));
              int i;
              for (i = 0; i < s; i++) ierr2 = fscanf(in,"%lf",&ct[i]); if (ierr2 != 1) return(1);
            }
          } else if (!strcmp(word,"c")) {
            if (s == 0) { 
              if (!rank) fprintf(stderr,"Error in PetscRegisterTIMethods(): nstages must be defined ");
              if (!rank) fprintf(stderr,"before specifying the Butcher tableaux entries.\n"              );
              return(1);
            } else {
              c = (PetscReal*) calloc (s, sizeof(PetscReal));
              int i;
              for (i = 0; i < s; i++) ierr2 = fscanf(in,"%lf",&c[i]); if (ierr2 != 1) return(1);
            }
          } else if (!strcmp(word,"bembt")) {
            if (s == 0) { 
              if (!rank) fprintf(stderr,"Error in PetscRegisterTIMethods(): nstages must be defined ");
              if (!rank) fprintf(stderr,"before specifying the Butcher tableaux entries.\n"              );
              return(1);
            } else {
              bembt = (PetscReal*) calloc (s, sizeof(PetscReal));
              int i;
              for (i = 0; i < s; i++) ierr2 = fscanf(in,"%lf",&bembt[i]); if (ierr2 != 1) return(1);
            }
          } else if (!strcmp(word,"bemb")) {
            if (s == 0) { 
              if (!rank) fprintf(stderr,"Error in PetscRegisterTIMethods(): nstages must be defined ");
              if (!rank) fprintf(stderr,"before specifying the Butcher tableaux entries.\n"              );
              return(1);
            } else {
              bemb = (PetscReal*) calloc (s, sizeof(PetscReal));
              int i;
              for (i = 0; i < s; i++) ierr2 = fscanf(in,"%lf",&bemb[i]); if (ierr2 != 1) return(1);
            }
          } else if (!strcmp(word,"bintt")) {
            if (s == 0 || pinterp == 0) { 
              if (!rank) fprintf(stderr,"Error in PetscRegisterTIMethods(): nstages and pinterp must be " );
              if (!rank) fprintf(stderr,"defined as positive values before specifying interpolation coeffs.\n");
              return(1);
            } else {
              bintt = (PetscReal*) calloc (s*pinterp, sizeof(PetscReal));
              int i, j;
              for (i = 0; i < s; i++) {
                for (j = 0; j < pinterp; j++) {
                  ierr2 = fscanf(in,"%lf",&bintt[i*s+j]); if (ierr2 != 1) return(1);
                }
              }
            }
          } else if (!strcmp(word,"bint")) {
            if (s == 0 || pinterp == 0) { 
              if (!rank) fprintf(stderr,"Error in PetscRegisterTIMethods(): nstages and pinterp must be " );
              if (!rank) fprintf(stderr,"defined as positive values before specifying interpolation coeffs.\n");
              return(1);
            } else {
              bint = (PetscReal*) calloc (s*pinterp, sizeof(PetscReal));
              int i, j;
              for (i = 0; i < s; i++) {
                for (j = 0; j < pinterp; j++) {
                  ierr2 = fscanf(in,"%lf",&bint[i*s+j]); if (ierr2 != 1) return(1);
                }
              }
            }
          }
        }
      } else {
    		if (!rank) fprintf(stderr,"Error: Illegal format in file \"time_method.inp\" (expected keyword \"begin\").\n");
        return(1);
      }
    
      /* Register the method */
      if (!strcmp(type,"arkimex")) {
        if (A && At) {
          ierr = TSARKIMEXRegister(name,order,s,At,bt,ct,A,b,c,bembt,bemb,pinterp,bintt,bint); CHKERRQ(ierr);
          if (!rank) {
            printf("\nRegistered custom ARKIMEX scheme \"%s\" with the following Butcher tableaux:-\n",name);
            int i,j;
            for (i = 0; i < s; i++) {
              if (c)  printf("  %+1.5lf |",c[i]);
              else    printf("           |");
              for (j = 0; j < s; j++) printf (" %+1.5lf :",A[i*s+j]);
              printf("          ");
              if (ct)  printf("%+1.5lf |",ct[i]);
              else     printf("           |");
              for (j = 0; j < s; j++) printf (" %+1.5lf :",At[i*s+j]);
              printf("\n");
            }
            printf("  ---------|");
            for (j = 0; j < s; j++) printf("-----------");
            printf("            ");
            printf("---------|");
            for (j = 0; j < s; j++) printf("-----------");
            printf("\n");
            printf("           |");
            if (b)   for (j = 0; j < s; j++) printf(" %+1.5lf :",b[j]);
            else     for (j = 0; j < s; j++) printf("          :");
            printf("          ");
            printf("           |");
            if (bt)  for (j = 0; j < s; j++) printf(" %+1.5lf :",bt[j]);
            else     for (j = 0; j < s; j++) printf("          :");
            printf("\n\n");
          }
        } else {
          if (!rank) {
            fprintf(stderr,"Warning in PetscRegisterTIMethods(): Failed to register method ");
            fprintf(stderr,"(A or At not defined).\n");
          }
        }
      } else if (!strcmp(type,"rk")) {
        if (A) {
          ierr = TSRKRegister(name,order,s,A,b,c,bemb,pinterp,bint); CHKERRQ(ierr);
          if (!rank) {
            printf("\nRegistered custom RK scheme \"%s\" with the following Butcher tableaux:-\n",name);
            int i,j;
            for (i = 0; i < s; i++) {
              if (c)  printf("  %+1.5lf |",c[i]);
              else    printf("           |");
              for (j = 0; j < s; j++) printf (" %+1.5lf :",A[i*s+j]);
              printf("\n");
            }
            printf("  ---------|");
            for (j = 0; j < s; j++) printf("-----------");
            printf("\n");
            printf("           |");
            if (b)   for (j = 0; j < s; j++) printf(" %+1.5lf :",b[j]);
            else     for (j = 0; j < s; j++) printf("          :");
            printf("\n\n");
          }
        } else {
          if (!rank) {
            fprintf(stderr,"Warning in PetscRegisterTIMethods(): Failed to register method ");
            fprintf(stderr,"(A not defined).\n");
          }
        }
      } else {
        if (!rank){
          fprintf(stderr,"Error in PetscRegisterTIMethods():  %s class of time-integration schemes ",type);
          fprintf(stderr,"does not support custom method registration and usage.\n");
        }
      }

      /* Free the arrays */
      if (At)     free(At);
      if (bt)     free(bt);
      if (ct)     free(ct);
      if (A)      free(A);
      if (b)      free(b);
      if (c)      free(c);
      if (bembt)  free(bembt);
      if (bemb)   free(bemb);
      if (bintt)  free(bintt);
      if (bint)   free(bint);
    
    }
  }
  PetscFunctionReturn(0);
}

#endif
