double GetCoordinate(int dir,int i,int *dim,int ghosts,double *x)
{
  int d,offset = 0;
  for (d = 0; d < dir; d++) offset += (dim[d]+2*ghosts);
  return(x[offset+ghosts+i]);
}
