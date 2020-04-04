/*! @file NavierStokes3DIBForces.c
    @brief Compute aerodynamic forces on immersed body, if present
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <immersedboundaries.h>
#include <physicalmodels/navierstokes3d.h>
#include <mpivars.h>
#include <hypar.h>

int NavierStokes3DComputePressure(double*, const double* const, void*);
int NavierStokes3DComputeTemperature(double*, const double* const, void*);

/*!
  Calculate the shear forces on the immersed body surface: for each
  "local facet", i.e., facets (of the immersed body surface) that lie within
  the local computational domain of this MPI rank, compute the shear forces 
  at the facet centroid from the flow variables at the grid points surrounding 
  the facet.

  The array to hold the computed shear forces *must* be NULL when this
  function is called. If the current subdomain contains a part of the
  immersed body, they will point to arrays with the local data. Otherwise,
  they will remain NULL.

  The shear force array will be of the size (4 X #ImmersedBoundary::nfacets_local),
  where the ordering of the data is: x-component, y-component, z-component, magnitude.
*/
static int ComputeShear(void *s,              /*!< Solver object of type #HyPar */
                        void *m,              /*!< MPI object of type #MPIVariables */
                        const double* const u,/*!< Array containing the conserved flow variables */ 
                        double** const sf     /*!< Array for (x,y,z)-components & magnitude of shear */
                       )
{
  HyPar             *solver  = (HyPar*)          s;
  MPIVariables      *mpi     = (MPIVariables*)   m; 
  NavierStokes3D    *physics = (NavierStokes3D*) solver->physics;
  ImmersedBoundary  *IB      = (ImmersedBoundary*) solver->ib;

  int nvars = solver->nvars;
  int ns = physics->n_species;
  int nv = physics->n_vibeng;

  if ((*sf) != NULL) {
    fprintf(stderr, "Error in ComputeShear()\n");
    fprintf(stderr, " shear force array is not NULL!\n");
    return 1;
  }

  if (!solver->flag_ib) return(0);

  int      nfacets_local = IB->nfacets_local;
  FacetMap *fmap = IB->fmap;
  
  double v[solver->nvars];

  int nshear = 4;

  if (nfacets_local > 0) {

    (*sf) = (double*) calloc (nshear*nfacets_local, sizeof(double));
  
    if (physics->Re > 0) {

      for (int n = 0; n < nfacets_local; n++) {
    
        double *alpha;
        int    *nodes, j, k;
    
        alpha = &(fmap[n].interp_coeffs[0]);
        nodes = &(fmap[n].interp_nodes[0]);
        _ArraySetValue_(v,solver->nvars,0.0);
        for (j=0; j<_IB_NNODES_; j++) {
          for (k=0; k<solver->nvars; k++) {
            v[k] += ( alpha[j] * u[solver->nvars*nodes[j]+k] );
          }
        }
        double rho_s_c[ns], rho_t_c, uvel_c, vvel_c, wvel_c, E_c, E_v_c[nv], pressure_c, T_c;
        _NavierStokes3DGetFlowVar_(v,rho_s_c,rho_t_c,uvel_c,vvel_c,wvel_c,E_c,E_v_c,pressure_c,T_c,physics);

        alpha = &(fmap[n].interp_coeffs_ns[0]);
        nodes = &(fmap[n].interp_nodes_ns[0]);
        _ArraySetValue_(v,solver->nvars,0.0);
        for (j=0; j<_IB_NNODES_; j++) {
          for (k=0; k<solver->nvars; k++) {
            v[k] += ( alpha[j] * u[solver->nvars*nodes[j]+k] );
          }
        }
        double rho_s_ns[ns], rho_t_ns, uvel_ns, vvel_ns, wvel_ns, E_ns, E_v_ns[nv], pressure_ns, T_ns;
        _NavierStokes3DGetFlowVar_(v,rho_s_ns,rho_t_ns,uvel_ns,vvel_ns,wvel_ns,E_ns,E_v_ns,pressure_ns,T_ns,physics);
        
        double u_x = (uvel_ns - uvel_c) / fmap[n].dx;
        double v_x = (vvel_ns - vvel_c) / fmap[n].dx;
        double w_x = (wvel_ns - wvel_c) / fmap[n].dx;
        
        double u_y = (uvel_ns - uvel_c) / fmap[n].dy;
        double v_y = (vvel_ns - vvel_c) / fmap[n].dy;
        double w_y = (wvel_ns - wvel_c) / fmap[n].dy;
        
        double u_z = (uvel_ns - uvel_c) / fmap[n].dz;
        double v_z = (vvel_ns - vvel_c) / fmap[n].dz;
        double w_z = (wvel_ns - wvel_c) / fmap[n].dz;
        
        double nx = fmap[n].facet->nx;
        double ny = fmap[n].facet->ny;
        double nz = fmap[n].facet->nz;
        
        double T      = physics->gamma*pressure_c/rho_t_c;
        double mu     = raiseto(T, 0.76);
        double inv_Re = 1.0/physics->Re;
  
        double tau_x = (mu*inv_Re) * (2*u_x*nx + (u_y+v_x)*ny + (u_z+w_x)*nz);
        double tau_y = (mu*inv_Re) * ((v_x+u_y)*nx + 2*v_y*ny + (v_z+w_y)*nz);
        double tau_z = (mu*inv_Re) * ((w_x+u_z)*nx + (w_y+v_z)*ny + 2*w_z*nz);
  
        (*sf)[n*nshear+_XDIR_] = tau_x;
        (*sf)[n*nshear+_YDIR_] = tau_y;
        (*sf)[n*nshear+_ZDIR_] = tau_z;

        (*sf)[n*nshear+_ZDIR_+1] = sqrt(tau_x*tau_x + tau_y*tau_y + tau_z*tau_z);
      }
  
    } else {

      _ArraySetValue_((*sf), nshear*nfacets_local, 0.0);
  
    }

  }

  return 0;
}

/*! Write the surface data on the immersed body to a ASCII Tecplot file. */
static int WriteSurfaceData(  void*               m,              /*!< MPI object of type #MPIVariables */
                              void*               ib,             /*!< Immersed body object of type #ImmersedBoundary */
                              const double* const p_surface,      /*!< array with local surface pressure data */
                              const double* const T_surface,      /*!< array with local surface temperature data */
                              const double* const ngrad_p_surface,/*!< array with local normal gradient of surface pressure data */
                              const double* const ngrad_T_surface,/*!< array with local normal gradient of surface temperature data */
                              const double* const shear,          /*!< array with local shear data */
                              char*               filename        /*!< Name of file to write */
                            )
{
  MPIVariables *mpi = (MPIVariables*) m;
  ImmersedBoundary *IB  = (ImmersedBoundary*) ib;
  int ierr;

#ifndef serial
  MPI_Status status;
#endif

  /* collect the surface data into global arrays */
  double* p_surface_g = NULL;
  ierr = IBAssembleGlobalFacetData(mpi, IB, p_surface, &p_surface_g, 1);
  if (ierr) {
    fprintf(stderr,"IBAssembleGlobalFacetData() returned with an error.\n");
    return 1;
  }
  double* T_surface_g = NULL;
  ierr = IBAssembleGlobalFacetData(mpi, IB, T_surface, &T_surface_g, 1);
  if (ierr) {
    fprintf(stderr,"IBAssembleGlobalFacetData() returned with an error.\n");
    return 1;
  }
  double* ngrad_p_surface_g = NULL;
  ierr = IBAssembleGlobalFacetData(mpi, IB, ngrad_p_surface, &ngrad_p_surface_g, 1);
  if (ierr) {
    fprintf(stderr,"IBAssembleGlobalFacetData() returned with an error.\n");
    return 1;
  }
  double* ngrad_T_surface_g = NULL;
  ierr = IBAssembleGlobalFacetData(mpi, IB, ngrad_T_surface, &ngrad_T_surface_g, 1);
  if (ierr) {
    fprintf(stderr,"IBAssembleGlobalFacetData() returned with an error.\n");
    return 1;
  }
  double* shear_g = NULL;
  ierr = IBAssembleGlobalFacetData(mpi, IB, shear, &shear_g, 4);
  if (ierr) {
    fprintf(stderr,"IBAssembleGlobalFacetData() returned with an error.\n");
    return 1;
  }

  /* Rank 0 writes the file */
	if (!mpi->rank) {

    int nfacets_global = IB->body->nfacets;
    const Facet3D* const facets = IB->body->surface;

		FILE *out;
    out = fopen(filename,"w");
    fprintf(out,"TITLE = \"Surface data created by HyPar.\"\n");
    fprintf(out,"VARIABLES = \"X\", \"Y\", \"Z\", ");
    fprintf(out,"\"Surface_Pressure\", ");
    fprintf(out,"\"Surface_Temperature\", ");
    fprintf(out,"\"Normal_Grad_Surface_Pressure\", ");
    fprintf(out,"\"Normal_Grad_Surface_Temperature\", ");
    fprintf(out,"\"Shear_x\", ");
    fprintf(out,"\"Shear_y\", ");
    fprintf(out,"\"Shear_z\", ");
    fprintf(out,"\"Shear_magn\"");
    fprintf(out,"\n");
    fprintf(out,"ZONE N = %d, E = %d, DATAPACKING = POINT, ZONETYPE = FETRIANGLE\n",3*nfacets_global,nfacets_global);

		for (int n = 0; n < nfacets_global; n++) {
      fprintf(  out, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                facets[n].x1,
                facets[n].y1,
                facets[n].z1,
                p_surface_g[n],
                T_surface_g[n],
                ngrad_p_surface_g[n],
                ngrad_T_surface_g[n],
                shear_g[4*n+_XDIR_],
                shear_g[4*n+_YDIR_],
                shear_g[4*n+_ZDIR_],
                shear_g[4*n+_ZDIR_+1] );
      fprintf(  out, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                facets[n].x2,
                facets[n].y2,
                facets[n].z2,
                p_surface_g[n],
                T_surface_g[n],
                ngrad_p_surface_g[n],
                ngrad_T_surface_g[n],
                shear_g[4*n+_XDIR_],
                shear_g[4*n+_YDIR_],
                shear_g[4*n+_ZDIR_],
                shear_g[4*n+_ZDIR_+1] );
      fprintf(  out, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                facets[n].x3,
                facets[n].y3,
                facets[n].z3,
                p_surface_g[n],
                T_surface_g[n],
                ngrad_p_surface_g[n],
                ngrad_T_surface_g[n],
                shear_g[4*n+_XDIR_],
                shear_g[4*n+_YDIR_],
                shear_g[4*n+_ZDIR_],
                shear_g[4*n+_ZDIR_+1] );
		}
		for (int n = 0; n < nfacets_global; n++) fprintf(out,"%d %d %d\n",3*n+1,3*n+2,3*n+3);
		fclose(out);
	}

  if (p_surface_g) free(p_surface_g);
  if (T_surface_g) free(T_surface_g);
  if (ngrad_p_surface_g) free(ngrad_p_surface_g);
  if (ngrad_T_surface_g) free(ngrad_T_surface_g);
  if (shear_g) free(shear_g);

  return 0;
}


/*! Calculate the aerodynamic forces on the immersed body surface and write them
    to file
*/
int NavierStokes3DIBForces( void *s, /*!< Solver object of type #HyPar */
                            void *m  /*!< MPI object of type #MPIVariables */
                          )
{
  HyPar             *solver  = (HyPar*)          s;
  MPIVariables      *mpi     = (MPIVariables*)   m; 
  NavierStokes3D    *physics = (NavierStokes3D*) solver->physics;
  ImmersedBoundary  *IB      = (ImmersedBoundary*) solver->ib;
  int ierr;

  if (!solver->flag_ib) return(0);

  int npts = solver->npoints_local_wghosts;

  double* pressure = (double*) calloc (npts, sizeof(double));
  ierr = NavierStokes3DComputePressure(pressure, solver->u, solver);
  if (ierr) {
    fprintf(stderr,"Error in NavierStokes3DIBForces()\n");
    fprintf(stderr,"  NavierStokes3DComputePressure() returned with error.\n");
    return 1;
  }
  double* temperature = (double*) calloc(npts, sizeof(double));
  ierr = NavierStokes3DComputeTemperature(temperature, solver->u, solver);
  if (ierr) {
    fprintf(stderr,"Error in NavierStokes3DIBForces()\n");
    fprintf(stderr,"  NavierStokes3DComputeTemperature() returned with error.\n");
    return 1;
  }

  /* Compute surface pressure */
  double* p_surface = NULL;
  ierr = IBComputeFacetVar(solver, mpi, pressure, 1, &p_surface);
  if (ierr) {
    fprintf(stderr,"Error in NavierStokes3DIBForces()\n");
    fprintf(stderr,"  IBComputeFacetVar() returned with error.\n");
    return 1;
  }
  /* Compute surface temperature */
  double* T_surface = NULL;
  ierr = IBComputeFacetVar(solver, mpi, temperature, 1, &T_surface);
  if (ierr) {
    fprintf(stderr,"Error in NavierStokes3DIBForces()\n");
    fprintf(stderr,"  IBComputeFacetVar() returned with error.\n");
    return 1;
  }
  /* Compute normal surface pressure gradient */
  double *ngrad_p_surface = NULL;
  ierr = IBComputeNormalGradient(solver, mpi, pressure, 1, &ngrad_p_surface);
  if (ierr) {
    fprintf(stderr,"Error in NavierStokes3DIBForces()\n");
    fprintf(stderr,"  IBComputeNormalGradient() returned with error.\n");
    return 1;
  }
  /* Compute normal temperature gradient */
  double *ngrad_T_surface = NULL;
  ierr = IBComputeNormalGradient(solver, mpi, temperature, 1, &ngrad_T_surface);
  if (ierr) {
    fprintf(stderr,"Error in NavierStokes3DIBForces()\n");
    fprintf(stderr,"  IBComputeNormalGradient() returned with error.\n");
    return 1;
  }
  /* Compute shear forces */
  double *shear = NULL;
  ierr = ComputeShear(solver, mpi, solver->u, &shear);
  if (ierr) {
    fprintf(stderr,"Error in NavierStokes3DIBForces()\n");
    fprintf(stderr,"  ComputeShear() returned with error.\n");
    return 1;
  }

  char surface_filename[_MAX_STRING_SIZE_] = "surface";
  if (!strcmp(solver->op_overwrite,"no")) {
    strcat(surface_filename,solver->filename_index);
  }
  strcat(surface_filename,".dat");
  if (!mpi->rank) {
    printf("Writing immersed body surface data file %s.\n",surface_filename);
  }
  ierr = WriteSurfaceData(  mpi,
                            IB,
                            p_surface,
                            T_surface,
                            ngrad_p_surface,
                            ngrad_T_surface,
                            shear,
                            surface_filename );
  if (ierr) {
    fprintf(stderr,"Error in NavierStokes3DIBForces()\n");
    fprintf(stderr,"  WriteSurfaceData() returned with error\n");
    return 1;
  }

  free(pressure);
  free(temperature);
  if (p_surface) free(p_surface);
  if (T_surface) free(T_surface);
  if (ngrad_p_surface) free(ngrad_p_surface);
  if (ngrad_T_surface) free(ngrad_T_surface);
  if (shear) free(shear);

  return 0;
}
