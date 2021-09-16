#include "option.h"

/*------------------------------------------------------------------------*/
/* Standard normal density function                                       */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.227                                           */
/*------------------------------------------------------------------------*/
double dN(double x) {
	static const double RT2PI = 2.5066283; /* Sqrt (2 * PI ) */

	return exp(-x*x/2.0)/RT2PI;
}

/*------------------------------------------------------------------------*/
/* Cummulative normal distribution function                               */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.227                                           */
/*------------------------------------------------------------------------*/
double N(double x) {
	static const double 	g = 0.2316419,
								a1 = 0.319381530,
								a2 = -0.356563782,
								a3 = 1.781477937,
								a4 = -1.821255978,
								a5 = 1.330274429;
	double k;

	if (x > 5.0)
		return 1.0;

	if (x >= 0)  {
		k = 1.0 / (1.0 + g * x);
		return  1.0 - dN(x) * k * (a1 + k * ( a2 + k * ( a3 + k * ( a4 + k * a5 ) ) ) );
	}
	else
		return  1.0 - N(-x);
}

/*---------------------------------------------------------------------------*/
/* The following routine was written by Boris Moro of TMG financial products */
/* & and was published as "The Full Monte, Risk Vol 8, #2, Feb 95, p.57"     */
/* This routine calculates inverse of cumulative normal distribution         */
/* function for argument u from the interval (0,1)                           */
/*  Good accuracy is achieved for 1.0E-10 < u < 1.0-1.0E-10,                 */
/* however even if u=O(1.0E-15) the result will still have at least first    */
/* three figures correct                                                     */
/*---------------------------------------------------------------------------*/
double Ninv(double u){
    
static double xsum  = 4.2454686881376569;
static double xdif  = 0.4179886424926431;

static double cofa[4]={
    2.50662823884,
  -18.61500062529,
   41.39119773534,
  -25.44106049637
  };
static double cofb[4]={
   -8.47351093090,
   23.08336743743,
  -21.06224101826,
    3.13082909833
  };


static double cof8[9]={
     7.7108870705487895,
     2.7772013533685169,
     0.3614964129261002,
     0.0373418233434554,
     0.0028297143036967,
     0.0001625716917922,
     0.0000080173304740,
     0.0000003840919865,
     0.0000000129707170
     };

double ndev,x,r,z,zz,yy,y;
int i;

x=u-0.5;
zz=fabs(x);
if(zz < 0.42)
  {
    r=x*x;
    ndev=x*(((cofa[3]*r+cofa[2])*r+cofa[1])*r+cofa[0])/
           ((((cofb[3]*r+cofb[2])*r+cofb[1])*r+cofb[0])*r+1.0);
    return(ndev);
  }

r=log(-log(0.5-zz));

r=(r + r - xsum)*xdif;

z=0.0;
zz=0.0;
        y=2.0*r;
        for (i=8;i>=1;i--)
          {
                yy=z;
                z=y*z-zz+cof8[i];
                zz=yy;
        }
        ndev=r*z-zz+0.5*cof8[0];
    if(x < 0.0)ndev=-ndev;
return(ndev);

}

/*------------------------------------------------------------------------*/
/* Support fn for Cummulative bivariate normal distribution below         */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.245                                           */
/*------------------------------------------------------------------------*/
static inline double f(double x, double y, double p, double ap, double bp) {

	return exp(ap * (2.0*x - ap) + bp * (2*y - bp) + 2.0*p*(x - ap)*(y-bp));
}

/*------------------------------------------------------------------------*/
/* Cummulative bivariate normal distribution                              */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.245                                           */
/*------------------------------------------------------------------------*/
double M(double a, double b, double p) {
	static const double A[] = {0.3253030, 0.4211071, 0.1334425, 0.006374323};
	static const double B[] = {0.1337764, 0.6243247, 1.3425378, 2.2626645};

	static const double pi = 3.141592654;
	static const double sqrt2 = 1.414213562; /* sqrt(2) */

	double ap, bp, sqrt_1_min_pp, denom, sum;
	double p1, p2, delta, sqrt_var;
	int i, j;

	if ( a <= 0.0 && b <= 0.0 && p <= 0.0) {
		if (p <= -1.0)
			p = -1.0 + 1.0e-15;
		sqrt_1_min_pp = sqrt(1.0-p*p);
		denom = sqrt2 * sqrt_1_min_pp;
		ap = a / denom;
		bp = b / denom;
		sum = 0.0;
		for (i=0; i<4; i++)
			for (j=0; j<4; j++)
				sum += A[i] * A[j] * f(B[i], B[j], p, ap, bp);
		sum *= sqrt_1_min_pp / pi;
		return sum;
	}

	if ( a <= 0.0 && b >= 0.0 && p >= 0.0)
		return N(a) - M(a, -b, -p);

	if ( a >= 0.0 && b <= 0.0 && p >= 0.0)
		return N(b) - M(-a, b, -p);

	if ( a >= 0.0 && b >= 0.0 && p <= 0.0)
		return N(a) + N(b) - 1.0 + M(-a, -b, p);

	if ( a * b * p > 0.0) {
		if (p >= 1.0)
			p = 1.0 - 1.0e-15;
		sqrt_var = sqrt(a*a - 2.0*p*a*b + b*b);
		p1 = (p*a - b)*SGN(a)/sqrt_var;
		p2 = (p*b - a)*SGN(b)/sqrt_var;
		delta = (1.0 - SGN(a)*SGN(b)) / 4.0;
		return M(a, 0.0, p1) + M(b, 0.0, p2) - delta;
	}

	/* The above should cover all the cases, if not we have an error */
	/* This line should never be reached */
	assert(FALSE);

	return 0.0;
}

/*------------------------------------------------------------------------*/
/* Cholesky Decomposition from "Numerical Recipes in C, 2nd Ed.,
	Press, Teukolsky et al. p. 97
	Given a positive definite symmetric matrix a[0...n-1][0...n-1] this
	routine constructs its Cholesky decomposition, A = L L'. On input,
	only the upper triangle of a need be given; it is not modified.  The
	Cholesky factor L is returned in the lower triangle of a, except for its
	diagonal elements which are returned in p[0...n-1]
 ------------------------------------------------------------------------*/
static short double_choldc(double **a, int n, double p[])
{
	int i,j,k;
	float sum;

	for (i=0;i<n;i++) {
		for (j=i;j<n;j++) {
			for (sum=a[i][j],k=i-1;k>=0;k--) sum -= a[i][k]*a[j][k];
			if (i == j) {
				if (sum <= 0.0)
					return FALSE; /* a with rounding errors, not +ve definite */
				p[i]=sqrt(sum);
			} else a[j][i]=sum/p[i];
		}
	}
	return TRUE;
}

/*------------------------------------------------------------------------*/
/* Simplified interface to call above Cholesky decomposition routine.
	Given a symmetric positive definite matrix a, the routine performs a
	Cholesky decomposition of a and returns the new matrix.  In the event
	that the matrix could not be decomposed, NULL is returned.  The values
	in a are not changed.
  ------------------------------------------------------------------------*/
double **Cholesky(double **a, int rows){
	double **chol_a, *p;
	int i, j;
	short status;

	chol_a = (double **)newMatrix(rows, rows, sizeof(double));
	if(chol_a == NULL)
		return NULL;

	/* Initialize upper triangle of chol_a with values in a */
	for (i=0; i<rows; i++)
		for(j=i; j<rows; j++)
			chol_a[i][j] = a[i][j];

	p = (double *) calloc(rows, sizeof(double));
	if (p==NULL)
		return NULL;

	status = double_choldc(chol_a, rows, p);
	if (status == FALSE) {
		deleteMatrix(chol_a);
		free(p);
		return NULL;
	}

	/* A true Cholesky matrix only has nonzero elements in lower/upper half */
	/* Zero out upper half of Cholesky matrix */
	for (i=0; i<rows; i++)
		for(j=i+1; j<rows; j++)
			chol_a[i][j] = 0.0;

#if 0
	/* Copy lower triangle of chol_a to uppper triangle of chol_a */
	for (i=0; i<rows; i++)
		for(j=0; j<i; j++)
			chol_a[j][i] = chol_a[i][j];
#endif

	/* Copy diagonal elements from p into chol_a */
	for (i=0; i<rows; i++)
		chol_a[i][i] = p[i];

	free(p);

	return chol_a;
}

void **newMatrix(unsigned long nrow, unsigned long ncol, size_t sElem)
/* allocate a matrix with subscript range m[0..nrows-1][0..ncols-1]
	based off of Numerical Recipies routines nrutil.c */
{
	unsigned long i;
	char **m;

	assert(nrow > 0);
	assert(ncol > 0);
	assert(sElem > 0);

	/* allocate rows & pointers to rows */
	m=(char **)calloc((size_t)(nrow + nrow*ncol), sElem);
	if (!m)  {
		/* MessageBox(NULL, "Unable to alloc memory for matrix!",
						"Positron Option Pricer",
						MB_ICONEXCLAMATION | MB_OK);
		*/
		return NULL;
	}

/* Adding long to pointer below may cause problem for very large numbers */
#pragma warn -sig /* Turn off warning */
	/* set pointers to rows */
	*m =(char *)( m + nrow );
#pragma warn +sig /* Turn on warning */

/* Adding long to pointer below may cause problem for very large numbers */
#pragma warn -sig /* Turn off warning */
	for(i=1;i<nrow;i++)
		m[i]=m[i-1]+ ncol*sElem;
#pragma warn +sig /* Turn on warning */

	/* return pointer to array of pointers to rows */
	return (void **)m;
}

void deleteMatrix(void *m)
/* free a matrix allocated by newMatrix() */
{
	free(m);
}

