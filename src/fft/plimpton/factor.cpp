/* parallel FFT functions - 1998, 1999

   Steve Plimpton, MS 1111, Dept 9221, Sandia National Labs
   (505) 845-7873
   sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory of the distribution.
*/

#include <math.h>

/* ------------------------------------------------------------------- */
/* recursively divide n into small factors, return them in list */

void factor(int n, int *num, int *list)

{
  if (n == 1) {
    return;
  } else if (n % 2 == 0) {
    *list = 2;
    (*num)++;
    factor(n/2,num,list+1);
  } else if (n % 3 == 0) {
    *list = 3;
    (*num)++;
    factor(n/3,num,list+1);
  } else if (n % 5 == 0) {
    *list = 5;
    (*num)++;
    factor(n/5,num,list+1);
  } else if (n % 7 == 0) {
    *list = 7;
    (*num)++;
    factor(n/7,num,list+1);
  } else if (n % 11 == 0) {
    *list = 11;
    (*num)++;
    factor(n/11,num,list+1);
  } else if (n % 13 == 0) {
    *list = 13;
    (*num)++;
    factor(n/13,num,list+1);
  } else {
    *list = n;
    (*num)++;
    return;
  }
}


/* ------------------------------------------------------------------- */
/* divide n into 2 factors of as equal size as possible */

void bifactor(int n, int *factor1, int *factor2)

{
  int n1,n2,facmax;

  facmax = sqrt((double) n);

  for (n1 = facmax; n1 > 0; n1--) {
    n2 = n/n1;
    if (n1*n2 == n) {
      *factor1 = n1;
      *factor2 = n2;
      return;
    }
  }
}
