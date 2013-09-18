inline int    ArrayAXPY              (double*,double,double*,int);

inline int    ArrayCopy1D_double     (double*,double*,int);
inline int    ArrayCopy1D_int        (int*,int*,int);
inline int    ArrayCopynD            (int,double*,double*,int*,int,int,int*,int);

inline int    ArrayIncrementIndex    (int,int*,int*);
inline int    ArrayIndex1D           (int,int*,int*,int*,int);

inline int    ArraySetValue_double   (double*,int,double);
inline int    ArraySetValue_int      (int*,int,int);

inline int    ArrayAdd1D_double      (double*,double*,double*,int);
inline int    ArrayAdd1D_int         (int*,int*,int*,int);
inline int    ArraySubtract1D_double (double*,double*,double*,int);
inline int    ArraySubtract1D_int    (int*,int*,int*,int);

inline double ArraySumSquarenD       (int,int,int*,int,int*,double*);
inline double ArraySumAbsnD          (int,int,int*,int,int*,double*);
inline double ArrayMaxnD             (int,int,int*,int,int*,double*);
