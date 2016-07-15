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

/*!
  Write the aerodynamic forces on the immersed body to a ASCII Tecplot file.
*/
static void WriteSurfaceData(
                              void      *m,             /*!< MPI object of type #MPIVariables */
                              int       nfacets_local,  /*!< Number of "local" facets on this MPI rank */
                              int       nfacets_global, /*!< Total number of facets in immersed body */
                              FacetMap  *fmap,          /*!< List of local facets and their info */
                              Facet3D   *facets,        /*!< Array of facets that define the immersed body */
                              double    *sp,            /*!< Array of size nfacets_local with the surface pressure */
                              double    *tx,            /*!< Array of size nfacets_local with the shear force in x */
                              double    *ty,            /*!< Array of size nfacets_local with the shear force in y */
                              double    *tz,            /*!< Array of size nfacets_local with the shear force in z */ 
                              double    *tm,            /*!< Array of size nfacets_local with the shear force magnitude */
                              char      *filename       /*!< Name of file to write */
                            )
{
  MPIVariables  *mpi = (MPIVariables*) m;

#ifndef serial
  MPI_Status status;
#endif

  /* Rank 0 collects and writes the file */
	if (!mpi->rank) {

		int n;

    /* allocate arrays for whole domain */
		double *sp_wd = (double*) calloc (nfacets_global, sizeof(double));
		double *tx_wd = (double*) calloc (nfacets_global, sizeof(double));
		double *ty_wd = (double*) calloc (nfacets_global, sizeof(double));
		double *tz_wd = (double*) calloc (nfacets_global, sizeof(double));
		double *tm_wd = (double*) calloc (nfacets_global, sizeof(double));

    /* check array - to avoid double counting of facets */
		int *check = (int*) calloc (nfacets_global, sizeof(int));
		for (n = 0; n < nfacets_global; n++) check[n] = 0;
		int check_total_facets = 0;

    /* local data */
		for (n = 0; n < nfacets_local; n++) {

			if (!check[fmap[n].index]) {

				sp_wd[fmap[n].index] = sp[n];
				tx_wd[fmap[n].index] = tx[n];
				ty_wd[fmap[n].index] = ty[n];
				tz_wd[fmap[n].index] = tz[n];
				tm_wd[fmap[n].index] = tm[n];
				check[fmap[n].index] = 1;

			} else {

        fprintf(stderr,"Error in WriteSurfaceData() (in NavierStokes3DIBForces.c): ");
        fprintf(stderr,"Facet %d has already been assigned a value. Double counting of facet!\n",fmap[n].index);

			}

		}
		check_total_facets += nfacets_local;

#ifndef serial
		for (int proc = 1; proc < mpi->nproc; proc++) {

			int nf_incoming;
			MPI_Recv(&nf_incoming, 1, MPI_INT, proc, 98927, MPI_COMM_WORLD, &status);
			check_total_facets += nf_incoming;

      if (nf_incoming > 0) {

			  int     *indices_incoming  = (int*)     calloc(nf_incoming, sizeof(int));
			  double  *sp_incoming       = (double*)  calloc(nf_incoming, sizeof(double));
			  double  *tx_incoming       = (double*)  calloc(nf_incoming, sizeof(double));
			  double  *ty_incoming       = (double*)  calloc(nf_incoming, sizeof(double));
			  double  *tz_incoming       = (double*)  calloc(nf_incoming, sizeof(double));
			  double  *tm_incoming       = (double*)  calloc(nf_incoming, sizeof(double));

			  MPI_Recv(indices_incoming , nf_incoming, MPI_INT    , proc, 98928, mpi->world, &status);
			  MPI_Recv(sp_incoming      , nf_incoming, MPI_DOUBLE , proc, 98929, mpi->world, &status);
			  MPI_Recv(tx_incoming      , nf_incoming, MPI_DOUBLE , proc, 98930, mpi->world, &status);
			  MPI_Recv(ty_incoming      , nf_incoming, MPI_DOUBLE , proc, 98931, mpi->world, &status);
			  MPI_Recv(tz_incoming      , nf_incoming, MPI_DOUBLE , proc, 98932, mpi->world, &status);
			  MPI_Recv(tm_incoming      , nf_incoming, MPI_DOUBLE , proc, 98933, mpi->world, &status);

			  for (n = 0; n < nf_incoming; n++) {
				  if (!check[indices_incoming[n]]) {
					  sp_wd[indices_incoming[n]] = sp_incoming[n];
					  tx_wd[indices_incoming[n]] = tx_incoming[n];
					  ty_wd[indices_incoming[n]] = ty_incoming[n];
					  tz_wd[indices_incoming[n]] = tz_incoming[n];
					  tm_wd[indices_incoming[n]] = tm_incoming[n];
					  check[indices_incoming[n]] = 1;
				  } else {
            fprintf(stderr,"Error in WriteSurfaceData() (in NavierStokes3DIBForces.c): ");
            fprintf(stderr,"Facet %d has already been assigned a value. Double counting of facet!\n",indices_incoming[n]);
				  }
			  }
        
        free(indices_incoming);
        free(sp_incoming);
        free(tx_incoming);
        free(ty_incoming);
        free(tz_incoming);
        free(tm_incoming);	

      }
		}
#endif

		if (check_total_facets != nfacets_global)	{
      fprintf(stderr,"Error in WriteSurfaceData() (in NavierStokes3DIBForces.c): mismatch in total facet count.\n");
    }

		FILE *out;
    out = fopen(filename,"w");
    fprintf(out,"TITLE = \"Surface data created by HyPar.\"\n");
    fprintf(out,"VARIABLES = \"X\", \"Y\", \"Z\", \"Surface_Pressure\", \"Shear_x\", \"Shear_y\", \"Shear_z\", \"Shear_magn\"\n");
    fprintf(out,"ZONE N = %d, E = %d, DATAPACKING = POINT, ZONETYPE = FETRIANGLE\n",3*nfacets_global,nfacets_global);

		for (n = 0; n < nfacets_global; n++) {
      fprintf(out, "%lf %lf %lf %lf %lf %lf %lf %lf\n",facets[n].x1,facets[n].y1,facets[n].z1,sp_wd[n],tx_wd[n],ty_wd[n],tz_wd[n],tm_wd[n]);
      fprintf(out, "%lf %lf %lf %lf %lf %lf %lf %lf\n",facets[n].x2,facets[n].y2,facets[n].z2,sp_wd[n],tx_wd[n],ty_wd[n],tz_wd[n],tm_wd[n]);
      fprintf(out, "%lf %lf %lf %lf %lf %lf %lf %lf\n",facets[n].x3,facets[n].y3,facets[n].z3,sp_wd[n],tx_wd[n],ty_wd[n],tz_wd[n],tm_wd[n]);
		}
		for (n = 0; n < nfacets_global; n++) fprintf(out,"%d %d %d\n",3*n+1,3*n+2,3*n+3);
		fclose(out);

		free(sp_wd);
		free(tx_wd);
		free(ty_wd);
		free(tz_wd);
		free(tm_wd);
    free(check);

	} else {
#ifndef serial
		MPI_Send(&nfacets_local, 1, MPI_INT, 0, 98927, MPI_COMM_WORLD);

    if (nfacets_local > 0) {

      int i, *indices = (int*) calloc (nfacets_local, sizeof(int));
      for (i = 0; i < nfacets_local; i++) indices[i] = fmap[i].index;

		  MPI_Send(indices, nfacets_local, MPI_INT    , 0, 98928, mpi->world);
		  MPI_Send(sp     , nfacets_local, MPI_DOUBLE , 0, 98929, mpi->world);
		  MPI_Send(tx     , nfacets_local, MPI_DOUBLE , 0, 98930, mpi->world);
		  MPI_Send(ty     , nfacets_local, MPI_DOUBLE , 0, 98931, mpi->world);
		  MPI_Send(tz     , nfacets_local, MPI_DOUBLE , 0, 98932, mpi->world);
		  MPI_Send(tm     , nfacets_local, MPI_DOUBLE , 0, 98933, mpi->world);

      free(indices);
    }
#endif
	}
}


/*!
  Calculate the aerodynamic forces on the immersed body surface: for each
  "local facet", i.e., facets (of the immersed body surface) that lie within
  the local computational domain of this MPI rank, compute the surface pressure
  and the shear forces at the facet centroid from the flow variables at the
  grid points surrounding the facet.
*/
int NavierStokes3DIBForces(
                            void *s, /*!< Solver object of type #HyPar */
                            void *m  /*!< MPI object of type #MPIVariables */
                          )
{
  HyPar             *solver  = (HyPar*)          s;
  MPIVariables      *mpi     = (MPIVariables*)   m; 
  NavierStokes3D    *physics = (NavierStokes3D*) solver->physics;
  ImmersedBoundary  *IB      = (ImmersedBoundary*) solver->ib;

  if (!solver->flag_ib) return(0);

  int           nfacets_local = IB->nfacets_local, n;
  FacetMap      *fmap = IB->fmap;
  static double v[_MODEL_NVARS_];
  double        *surface_pressure,
                *shear_force_x,
                *shear_force_y,
                *shear_force_z,
                *shear_force_magn;

  if (nfacets_local > 0) {

    surface_pressure = (double*) calloc (nfacets_local, sizeof(double));
    shear_force_x    = (double*) calloc (nfacets_local, sizeof(double));
    shear_force_y    = (double*) calloc (nfacets_local, sizeof(double));
    shear_force_z    = (double*) calloc (nfacets_local, sizeof(double));
    shear_force_magn = (double*) calloc (nfacets_local, sizeof(double));

  } else {

    surface_pressure = NULL;
    shear_force_x    = NULL;
    shear_force_y    = NULL;
    shear_force_z    = NULL;
    shear_force_magn = NULL;

  }

  for (n = 0; n < nfacets_local; n++) {

    double *alpha;
    int    *nodes, j, k;

    alpha = &(fmap[n].interp_coeffs[0]);
    nodes = &(fmap[n].interp_nodes[0]);
    _ArraySetValue_(v,_MODEL_NVARS_,0.0);
    for (j=0; j<_IB_NNODES_; j++) {
      for (k=0; k<_MODEL_NVARS_; k++) {
        v[k] += ( alpha[j] * solver->u[_MODEL_NVARS_*nodes[j]+k] );
      }
    }
    double rho_c, uvel_c, vvel_c, wvel_c, energy_c, pressure_c;
    _NavierStokes3DGetFlowVar_(v,rho_c,uvel_c,vvel_c,wvel_c,energy_c,pressure_c,physics);

    alpha = &(fmap[n].interp_coeffs_ns[0]);
    nodes = &(fmap[n].interp_nodes_ns[0]);
    _ArraySetValue_(v,_MODEL_NVARS_,0.0);
    for (j=0; j<_IB_NNODES_; j++) {
      for (k=0; k<_MODEL_NVARS_; k++) {
        v[k] += ( alpha[j] * solver->u[_MODEL_NVARS_*nodes[j]+k] );
      }
    }
    double rho_ns, uvel_ns, vvel_ns, wvel_ns, energy_ns, pressure_ns;
    _NavierStokes3DGetFlowVar_(v,rho_ns,uvel_ns,vvel_ns,wvel_ns,energy_ns,pressure_ns,physics);

    surface_pressure[n] = pressure_c;

    if (physics->Re > 0) {
      
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
      
      double T      = physics->gamma*pressure_c/rho_c;
      double mu     = raiseto(T, 0.76);
      double inv_Re = 1.0/physics->Re;

      double tau_x = (mu*inv_Re) * (2*u_x*nx + (u_y+v_x)*ny + (u_z+w_x)*nz);
      double tau_y = (mu*inv_Re) * ((v_x+u_y)*nx + 2*v_y*ny + (v_z+w_y)*nz);
      double tau_z = (mu*inv_Re) * ((w_x+u_z)*nx + (w_y+v_z)*ny + 2*w_z*nz);

      shear_force_x[n] = tau_x;
      shear_force_y[n] = tau_y;
      shear_force_z[n] = tau_z;

    } else {

      shear_force_x[n] = 0.0;
      shear_force_y[n] = 0.0;
      shear_force_z[n] = 0.0;

    }
    shear_force_magn[n] = sqrt(   (shear_force_x[n])*(shear_force_x[n]) 
                                + (shear_force_y[n])*(shear_force_y[n]) 
                                + (shear_force_z[n])*(shear_force_z[n]) );

  }

  char surface_filename[_MAX_STRING_SIZE_] = "surface";
  if (!strcmp(solver->op_overwrite,"no")) strcat(surface_filename,solver->filename_index);
  strcat(surface_filename,".dat");

  if (!mpi->rank) {
    printf("Writing immersed body surface data file %s.\n",surface_filename);
  }
  WriteSurfaceData(mpi,nfacets_local,IB->body->nfacets,fmap,IB->body->surface,
                   surface_pressure,
                   shear_force_x,
                   shear_force_y,
                   shear_force_z,
                   shear_force_magn,
                   surface_filename);

  if (nfacets_local > 0) {
    free(surface_pressure);
    free(shear_force_x);
    free(shear_force_y);
    free(shear_force_z);
    free(shear_force_magn);
  }

  return(0);
}
