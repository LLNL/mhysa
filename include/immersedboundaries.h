/*! @file immersedboundaries.h
    @brief Structures and function definitions for immersed boundaries
    @author Debojyoti Ghosh
*/

/*! Immersed boundaries are implemented only for 3-dimensional simulations. */
#define _IB_NDIMS_ 3
/*! Number of grid points surrounding a random point in space, i.e., number of
    corners of a cube. */
#define _IB_NNODES_ 8

/*! \def Facet3D
    \brief Structure defining a facet.
    
    A "facet" is the basic unit of a tessellated surface in 3D: It is a triangle
    in an unstructured triangulated surface, defined by its vertices and surface
    normal. \sa https://en.wikipedia.org/wiki/STL_(file_format)
*/
/*!\brief Structure defining a facet.
    
    A "facet" is the basic unit of a tessellated surface in 3D: It is a triangle
    in an unstructured triangulated surface, defined by its vertices and surface
    normal. \sa https://en.wikipedia.org/wiki/STL_(file_format)
*/
typedef struct _facet_3d_{
	double x1, /*!< x-coordinate of vertex 1 */
         x2, /*!< x-coordinate of vertex 2 */ 
         x3, /*!< x-coordinate of vertex 3 */
	       y1, /*!< y-coordinate of vertex 1 */
         y2, /*!< y-coordinate of vertex 2 */ 
         y3, /*!< y-coordinate of vertex 3 */
	       z1, /*!< z-coordinate of vertex 1 */
         z2, /*!< z-coordinate of vertex 2 */ 
         z3, /*!< z-coordinate of vertex 3 */
	       nx, /*!< x-component of surface normal */
         ny, /*!< y-component of surface normal */
         nz; /*!< z-component of surface normal */
} Facet3D;

/*! \def Body3D
    \brief Structure defining a body.

    A 3D body whose surface is represented as a collection 
    of faces of type #Facet3D.
    \sa https://en.wikipedia.org/wiki/STL_(file_format)
*/
/*! \brief Structure defining a body.

    A 3D body whose surface is represented as a collection 
    of faces of type #Facet3D.
    \sa https://en.wikipedia.org/wiki/STL_(file_format)
*/
typedef struct _body_3d_{
  int     nfacets;    /*!< number of surface facets   */
  Facet3D *surface;   /*!< array of surface facets    */
  /* coordinates of bounding box */
  double xmin, /*!< x-coordinate of lower end of bounding box */
         xmax, /*!< x-coordinate of higher end of bounding box */
         ymin, /*!< y-coordinate of lower end of bounding box */
         ymax, /*!< y-coordinate of higher end of bounding box */
         zmin, /*!< z-coordinate of lower end of bounding box */
         zmax; /*!< z-coordinate of higher end of bounding box */
} Body3D;

/*! \def IBNode
    \brief Structure defining an immersed boundary node

    An immersed boundary node is a grid point inside the immersed
    body but within stencil-width-distance of a point outside
    the body.
*/
/*! \brief Structure defining an immersed boundary node

    An immersed boundary node is a grid point inside the immersed
    body but within stencil-width-distance of a point outside
    the body.
*/
typedef struct _boundary_node_{
  int     i,      /*!< grid index along x */
          j,      /*!< grid index along y */
          k,      /*!< grid index along z */
          p;      /*!< array index of this point */
  double  x,      /*!< x-coordinate of the boundary node */
          y,      /*!< y-coordinate of the boundary node */
          z;      /*!< z-coordinate of the boundary node */
  Facet3D *face;  /*!< the nearest facet of #ImmersedBoundary::body */

  int     interp_nodes[_IB_NNODES_];   /*!< indices of interior nodes from which to extrapolate */
  double  interp_coeffs[_IB_NNODES_],  /*!< interpolation coefficients corresponding to #IBNode::interp_nodes */
          interp_node_distance,        /*!< Distance from the immersed body surface to the interior node from which to extrapolate */
          surface_distance;            /*!< Distance from this node to the immersed body surface */

} IBNode;

/*! \def ImmersedBoundary
    \brief Structure containing variables for immersed boundary implementation.

    This structure contains all the variables needed to implement immersed
    boundaries.
*/
/*! \brief Structure containing variables for immersed boundary implementation.

    This structure contains all the variables needed to implement immersed
    boundaries.
*/
typedef struct immersed_boundary{
  Body3D  *body;      /*!< immersed body */
  IBNode  *boundary;  /*!< immersed boundary nodes */

  double  tolerance; /*!< zero tolerance */
  int     itr_max;   /*!< maximum intersections in ray-tracing method */
  int     n_boundary_nodes; /*!< number of immersed boundary nodes */
} ImmersedBoundary;


int IBReadBodySTL (Body3D**,char*,void*,int*);
int IBWriteBodySTL(Body3D*,char*,void*,int,int*);

int IBCleanup           (void*);
int IBComputeBoundingBox(Body3D*);
int IBIdentifyBody      (void*,int*,int*,int,void*,double*,double*);
int IBIdentifyBoundary  (void*,void*,int*,int,double*);
int IBNearestFacetNormal(void*,void*,double*,double,int*,int);
int IBInterpCoeffs      (void*,void*,double*,int*,int,double*);
