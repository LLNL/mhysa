/* 

  Parallel direct solver for tridiagonal systems 

  tridiagLU  (a,b,c,x,n,ns,r,m) - Parallel tridiagonal solver
    
    Solves the tridiagonal system in parallel by reordering the
    points such that the first point of each subdomain is placed
    at the end.

    The interior points are eliminated in parallel, resulting in
    a reduced system consisting of the first point of each sub-
    domain.

    This reduced system is solved either by the gather-and-
    solve (tridiagLUGS) or the recursive-doubling (tridiagLURD)
    algorithms.

  tridiagLUGS(a,b,c,x,n,ns,r,m) - Tridiagonal solver based on
                                  "gather and solve"

    Each of the "ns" systems is gathered on one processor, 
    solved in serial, and the solution scattered back. The
    parallelism lies in solving the "ns" different systems 
    by multiple processors (i.e., each processor solves 
    ~ns/nproc number of systems in serial).

  Arguments:-
    a   [0,ns-1]x[0,n-1] double*          subdiagonal entries
    b   [0,ns-1]x[0,n-1] double*          diagonal entries
    c   [0,ns-1]x[0,n-1] double*          superdiagonal entries
    x   [0,ns-1]x[0,n-1] double*          right-hand side (solution)
    n                    int              local size of the system
    ns                   int              number of systems to solve
    r                    TridiagLU*       structure containing paramters
                                          for the tridiagonal solve and
                                          the walltimes:
                                            total_time
                                            stage1_time
                                            stage2_time
                                            stage3_time
                                            stage4_time
                        ** Note that these are process-specific. Calling 
                           function needs to do something to add/average 
                           them to get some global value.
                        ** Can be NULL if runtimes are not needed.
    m                   MPI_Comm*       MPI Communicator

  For a,b,c, and x, [0,ns-1] is the inner loop, i.e., the i-th row of the
  d-th system is a[i*ns+d], b[i*ns+d], c[i*ns+d] and x[i*ns+d].

  Return value (int) -> 0 (successful solve), -1 (singular system)

  Note:-
    x contains the final solution (right-hand side is replaced)
    a,b,c are not preserved
    On rank=0,        a[0*ns+d] has to be zero for all d.
    On rank=nproc-1,  c[(n-1)*ns+d] has to be zero for all d.

  For a serial tridiagonal solver, compile with the flag "-Dserial"
  or send NULL as the argument for the MPI communicator.

*/

/* definitions */
#define _TRIDIAG_JACOBI_  "jacobi"
#define _TRIDIAG_GS_      "gather-and-solve"

/* Data structure containing the stage runtimes */
typedef struct _tridiagLU_ {

  /* Parameters for tridiagLU() */
  char reducedsolvetype[50]; /* solver to use for the reduced system */
                             /* gather-and-solve direct solver       */
                             /* or iterative jabobi solver           */

  /* Parameters for the iterative tridiagonal solver    */
  int     evaluate_norm;  /* calculate norm at each iteration?*/
  int     maxiter;        /* maximum number of iterations     */
  double  atol,rtol;      /* absolute and relative tolerace   */
  int     exititer;       /* number of iterations it ran for  */
  double  exitnorm;       /* error norm at exit               */
  int     verbose;        /* print iterations and norms       */

  /* Runtimes for the 4 stages of tridiagLU() */
  double  total_time;
  double  stage1_time;
  double  stage2_time;
  double  stage3_time;
  double  stage4_time;

} TridiagLU;

int tridiagLU         (double*,double*,double*,double*,int,int,void*,void*);
int tridiagLUGS       (double*,double*,double*,double*,int,int,void*,void*);
int tridiagIterJacobi (double*,double*,double*,double*,int,int,void*,void*);
int tridiagLUInit     (void*,void*);

/* Block solvers */
int blocktridiagLU         (double*,double*,double*,double*,int,int,int,void*,void*);
int blocktridiagIterJacobi (double*,double*,double*,double*,int,int,int,void*,void*);
