typedef struct solution_physics {

  /* CFL computation function */
  double (*ComputeCFL) ();

} SolverPhysics;
