#include "option.h"

static int m_init = 0;

/*----------------------------------------------------------*/
/* Faure sequence generator, to generate quasirandom        */
/* Normal deviates.                                         */
/* Usage:                                                   */
/* Generate n random variables and store them in out.       */
/* MaxN= maximum number of path simulations and out[]       */
/* represents one sample path. I.E. Out is one vector from  */
/* the Faure sequence of length s.                          */
/* This program is a direct translation of the FORTRAN      */
/* routine by Fox, P. "ALGORITHM 647: Implementation and    */
/* Relative Efficiency of Quasirandom Sequence Generators." */
/* ACM Trans. Math. Software, 12 (1986), p. 362-376         */
/*----------------------------------------------------------*/
 
#define NORMAL                      /* -- This causes getrvsf to return Normal
                                     * rvs by the inverse cum norm function */



short rvsf(double out[], int s, unsigned long MaxN, SEQTYPE sequence_type)
{

#define MAXCHOOSE 20

		  static int      maxdigits, r, choose[MAXCHOOSE][MAXCHOOSE],
								y[MAXCHOOSE], hisum;
		  static unsigned long StartN, nextn, digit_flip, n;
        static double   rinv;

        double          rtemp, *x;

        int             i, j, k, sumy;

#define MAXPRIME 59
        static int      primes[MAXPRIME] = {2, 2, 3, 5, 5, 7, 7, 11, 11, 11, 11, 13, 13, 17, 17, 17,
                17, 19, 19, 23, 23, 23, 23, 29, 29, 29, 29, 29, 29, 31, 31, 37, 37, 37, 37, 37,
        37, 41, 41, 41, 41, 43, 43, 47, 47, 47, 47, 53, 53, 53, 53, 53, 53, 59, 59, 59, 59, 59, 59};

		  ldiv_t           divn;


		  if (sequence_type == SEQ_EXACT)
					 s--;
		  else if (sequence_type == SEQ_INVALID)
					 return FAIL;

        if (m_init == 0) {

                if (s > MAXPRIME) {
								/* cerr<<"Error! Number of time steps"<<s<<" is too large!\n"; */
                        return FAIL;
                }

                r = primes[s - 1];      /* Find the smallest prime that is >=
                                         * to s */

					 if ( r > 3 )
						StartN = r * (r-1) * (r-2) ;
					 else
						StartN = r * r * r;

					 assert(MaxN + StartN > 0);
					 assert(r > 0);
					 maxdigits = log(MaxN + StartN) / log(r);

                if (maxdigits > MAXCHOOSE) { 
                  /* cerr<<"Error! Too many iterations requested of the Faure routine \n";  */
                  return FAIL;
                }

                /* Compute Choose Functions modulo r */
                choose[0][0] = 1;
                for (i = 1; i <= maxdigits; i++) {
                        choose[i][0] = 1;
                        choose[i][i] = 1;
                }

                /*
                 * Generate Remaining Coefficients via the binomial tree
                 * recursion choose(i,j)=choose(i-1,j)+choose(i-1,j-1) and
                 * using the fact that mod(b+c,d)=mod(mod(b,d)+mod(c,d),d)
                 */
                for (j = 1; j <= maxdigits; j++)
                        for (i = 1 + j; i <= maxdigits; i++) {
                                choose[i][j] = (choose[i - 1][j] + choose[i - 1][j - 1]) % r;
                        }



					 nextn = StartN - 1;     /* start out at r^3-1 to avoid
													  * clustering near zero */
					 assert(nextn > 0);
					 assert(r > 0);
                hisum = (int) floor(log(nextn)/log(r)); /* calc # of digits for nextn in base r */
                digit_flip = StartN/r;  /* Increment Hisum when last 
                                         * digit is zero */

                m_init = 1;
		  }                       /* End initialization */

		  if (sequence_type == SEQ_EXACT) {
				x = out;
 		  /* Note order of operations here to prevent underflow with unsigned long */
#ifdef NORMAL
				out[0] = Ninv(((nextn + 1.5)- StartN) / (double)MaxN);
#else
				out[0] = ((nextn + 1.5) - StartN) / (double)MaxN;
#endif
		  } else
				x = out - 1;         /* Offset to return array starting at zero
                                 * instead of 1 */

        /* Compute First Element in the nth Faure Sequence */
        n = nextn;
        for (j = 0; j <= hisum; j++) {
					 divn = ldiv(n, (long)r);
					 y[j] = (int) (divn.rem);
                n = divn.quot;
        }

        /* Compute x(1) via sum(i=0, hisum) { y(i)r^(-i-1) } */
        rinv = 1.0 / r;
        x[1] = y[hisum];
        for (i = hisum - 1; i >= 0; i--)
                x[1] = y[i] + x[1] * rinv;
        x[1] *= rinv;

#ifdef FAURE
        printf("faure[1]=%f ", x[1]);
#endif

#ifdef NORMAL
		  x[1] = Ninv(x[1]);
#endif

        /* Apply the map C via New y(j)= sum Old y(j) * choose(i,j) mod r  */

        for (k = 2; k <= s; k++) {

                        x[k] = 0;

                rtemp = rinv;
                for (j = 0; j <= hisum; j++) {
                        sumy = 0;
                        for (i = j; i <= hisum; i++)
                                sumy += choose[i][j] * y[i];
                        y[j] = sumy % r;
                        x[k] += ((double) y[j]) * rtemp;
                        rtemp *= rinv;
                }               /* end for j */

#ifdef FAURE
                printf("faure[%d]=%f ", k, x[k]);
#endif

#ifdef NORMAL
					 x[k] = Ninv(x[k]);
#endif

        }                       /* end for k */

        nextn++;
        if ((nextn % digit_flip) == 0) {
                hisum++;
                if (hisum > MAXCHOOSE) {
						 /* cerr<<"Error! Max number of digits for choose exceeded \n";  */
                   return FAIL;
                }
                digit_flip *= r;
        }

#ifdef FAURE
		  printf("\n");
#endif

		  return OK;               /* Everything went O.K. */

}                               /* end rvsf */


/*----------------------------------------------------------*/
/* Reset the Faure Sequence Generator to start over at the  */
/* beginning                                                */
/*----------------------------------------------------------*/
void rvsf_reset(void) {
	m_init = 0;
}


