#include <stdio.h>

void IncrementFilename(char *f)
{
  if (f[7] == '9') {
    f[7] = '0';
    if (f[6] == '9') {
      f[6] = '0';
      if (f[5] == '9') {
        f[5] = '0';
        if (f[4] == '9') {
          f[4] = '0';
          if (f[3] == '9') {
            f[3] = '0';
            fprintf(stderr,"Warning: file increment hit max limit. Resetting to zero.\n");
          } else {
            f[3]++;
          }
        } else {
          f[4]++;
        }
      } else {
        f[5]++;
      }
    } else {
      f[6]++;
    }
  } else {
    f[7]++;
  }
}

void IncrementFilenameIndex(char *f,int len)
{
  int i;
  for (i=len-1; i>=0; i--) {
    if (f[i] == '9') {
      f[i] = '0';
      if (!i) fprintf(stderr,"Warning: file increment hit max limit. Resetting to zero.\n");
    } else {
      f[i]++;
      break;
    }
  }
}
