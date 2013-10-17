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
    a   [0,ns-1]x[0,n-1] double**         subdiagonal entries
    b   [0,ns-1]x[0,n-1] double**         diagonal entries
    c   [0,ns-1]x[0,n-1] double**         superdiagonal entries
    x   [0,ns-1]x[0,n-1] double**         right-hand side (solution)
    n                    int              local size of the system
    ns                   int              number of systems to solve
    r                    TridiagLUTime*   structure containing the runtimes
                                            total_time
                                            stage1_time
                                            stage2_time
                                            stage3_time
                                            stage4_time
                        ** Note that these are process-specific. Calling 
                           function needs to do something to add/average 
                           them to get some global value.
                        ** Can be NULL if runtimes are not needed.
    m                   MPIContext*       structure containing the MPI
                                          context
                                          **See below

  Return value (int) -> 0 (successful solve), -1 (singular system)

  Note:-
    x contains the final solution (right-hand side is replaced)
    a,b,c are not preserved
    On rank=0,        a[0] has to be zero.
    On rank=nproc-1,  c[n-1] has to be zero.

  For a serial tridiagonal solver, compile with the flag "-Dserial"
  or send NULL as the argument for the MPI communicator.

*/


/* Data structure containing the stage runtimes */
typedef struct _tridiagLUruntimes_ {
  double  total_time;
  double  stage1_time;
  double  stage2_time;
  double  stage3_time;
  double  stage4_time;
} TridiagLUTime;


/* 
  Data structure for the MPI context

  rank      rank of this process with respect to the processes parti-
            cipating in the tridiagonal solve
  nproc     number of processes participating in the tridiagonal solve
  comm      MPI communicator
  proc      an array of size nproc containing the actual rank in comm
            for each rank 0,...,nproc
*/

typedef struct _mpi_context_ {
  int   rank;
  int   nproc;
  void* comm;
  int*  proc;
} MPIContext;

int tridiagLU  (double**,double**,double**,double**,int,int,void*,void*);
int tridiagLUGS(double**,double**,double**,double**,int,int,void*,void*);
