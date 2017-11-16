/*! @file IncrementFilename.c
    @author Debojyoti Ghosh
    @brief Functions for incrementing filename indices
*/

#include <stdio.h>

/*! Increment the output filename of the form op_nnnnn.dat by 1,
    i.e., the number represented by the string nnnnn is incremented 
    by 1. 
*/
void IncrementFilename(char *f /*!< filename */)
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

/*! Increment a character string representing an integer by 1. For example:
    "00002" -> "00003"; "3421934" -> "3421935"; "999" -> "000". The string
    can be of arbitrary length.
*/
void IncrementFilenameIndex(
                              char *f,  /*!< Character string representing the integer */
                              int len   /*!< Length of the string */
                           )
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
