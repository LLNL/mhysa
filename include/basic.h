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

