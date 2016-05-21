/*! @file IBComputeBoundingBox.c
    @author Debojyoti Ghosh
    @brief Compute bounding box for an immersed body
*/

#include <immersedboundaries.h>

/*! Compute the bounding box for a given body. */
int IBComputeBoundingBox(Body3D *b /*!< The body */)
{
  b->xmin = b->xmax = b->surface[0].x1;
  b->ymin = b->ymax = b->surface[0].y1;
  b->zmin = b->zmax = b->surface[0].z1;

  int n;
  for (n = 0; n < b->nfacets; n++) {
    if (b->surface[n].x1 < b->xmin) b->xmin = b->surface[n].x1;
    if (b->surface[n].x2 < b->xmin) b->xmin = b->surface[n].x2;
    if (b->surface[n].x3 < b->xmin) b->xmin = b->surface[n].x3;

    if (b->surface[n].y1 < b->ymin) b->ymin = b->surface[n].y1;
    if (b->surface[n].y2 < b->ymin) b->ymin = b->surface[n].y2;
    if (b->surface[n].y3 < b->ymin) b->ymin = b->surface[n].y3;

    if (b->surface[n].z1 < b->zmin) b->zmin = b->surface[n].z1;
    if (b->surface[n].z2 < b->zmin) b->zmin = b->surface[n].z2;
    if (b->surface[n].z3 < b->zmin) b->zmin = b->surface[n].z3;

    if (b->surface[n].x1 > b->xmax) b->xmax = b->surface[n].x1;
    if (b->surface[n].x2 > b->xmax) b->xmax = b->surface[n].x2;
    if (b->surface[n].x3 > b->xmax) b->xmax = b->surface[n].x3;

    if (b->surface[n].y1 > b->ymax) b->ymax = b->surface[n].y1;
    if (b->surface[n].y2 > b->ymax) b->ymax = b->surface[n].y2;
    if (b->surface[n].y3 > b->ymax) b->ymax = b->surface[n].y3;

    if (b->surface[n].z1 > b->zmax) b->zmax = b->surface[n].z1;
    if (b->surface[n].z2 > b->zmax) b->zmax = b->surface[n].z2;
    if (b->surface[n].z3 > b->zmax) b->zmax = b->surface[n].z3;
  }
  return(0);
}
