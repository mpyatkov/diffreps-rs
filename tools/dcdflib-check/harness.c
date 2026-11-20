#include <stdio.h>
#include <stdlib.h>
#include "cdflib.h"
/* mode b: "a b x" -> cumbet(x, 1-x, a, b); mode n: "s xn pr" -> cdfnbn which=1 */
int main(void) {
  char mode; double u, v, w;
  while (scanf(" %c %lf %lf %lf", &mode, &u, &v, &w) == 4) {
    if (mode == 'b') {
      double x = w, y = 1.0 - w, cum, ccum;
      cumbet(&x, &y, &u, &v, &cum, &ccum);
      printf("%a %a\n", cum, ccum);
    } else {
      int which = 1, status; double p, q, s = u, xn = v, pr = w, ompr = 1.0 - w, bound;
      cdfnbn(&which, &p, &q, &s, &xn, &pr, &ompr, &status, &bound);
      if (status == 0) printf("%a\n", p); else printf("undef\n");
    }
  }
  return 0;
}
