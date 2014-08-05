#define _MAX_STRING_SIZE_ 500
#ifdef debug
  #define IERR ierr= 
  #define _DECLARE_IERR_ int ierr = 0; 
  #define CHECKERR(ierr) { if (ierr) return(ierr); }
#else
  #define IERR
  #define _DECLARE_IERR_
  #define CHECKERR(ierr)
#endif

/* Macro to get the coordinate value from the array containing the coordinate values
 * along all the dimensions.
 * Arguments:-
 *  dir       : coordinate dimension (x is 0, y is 1, z is 2, ...) (int)
 *  i         : grid index along that dimension (int)
 *  dim       : array of the grid sizes in each dimension (int[])
 *  ghosts    : number of ghost points 
 *  x         : array containing the coordinate values in all dimensions
 * Output:
 *  coord     : the coordinate value
*/
#define _GetCoordinate_(dir,i,dim,ghosts,x,coord) \
  { \
    int arraycounter,offset = 0; \
    for (arraycounter = 0; arraycounter < dir; arraycounter++) offset += (dim[arraycounter]+2*ghosts); \
    coord = (x[offset+ghosts+i]); \
  }


