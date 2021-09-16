#ifndef _OPTION_H_
#define _OPTION_H_

#include <math.h>
#include <stdlib.h>
#include <assert.h>
#define FALSE 0
#define TRUE 1
#include "date.h"

#define FAIL FALSE
#define OK TRUE
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#define SGN(A) ((A) >= 0.0 ? 1.0 : -1.0)

/* constant used to determine when vol, and time parameters have underflowed */
#define INPUT_UNDERFLOW 0.00001

/* Machine precision for use in calculating derivatives via */
/* finite differences.  Note that this should be chosen to  */
/* be an exactly represantable binary # for the best        */
/* accuracy.  See "Numerical Recipes in C, 2nd Ed." for a   */
/* discussion of this issue                                 */
const unsigned int eps = 2 << 13;

enum OPTTYPE
{
    OPT_INVALID,
    OPT_CALL,
    OPT_PUT
};
enum SENSTYPE
{
    SENS_ALLSENS = -1,
    SENS_INVALID,
    SENS_PRICE,
    SENS_DELTA,
    SENS_GAMMA,
    SENS_VEGA,
    SENS_THETA,
    SENS_RHO,
    SENS_DELTA2,
    SENS_GAMMA2,
    SENS_VEGA2,
    SENS_CORR,
    SENS_MAX
};

enum SEQTYPE
{
    SEQ_INVALID,
    SEQ_EXACT,
    SEQ_EXTENDIBLE
};

/* --------------------------------------------------------------------*/
/* amer.cpp */
/* --------------------------------------------------------------------*/
short amer(double S, double X, double vol, double r, double q, double t,
           OPTTYPE call_put, SENSTYPE sens, double *result);
short amer_vol(double S, double X, double price, double r, double q, double t,
               OPTTYPE call_put, double *result);

/* --------------------------------------------------------------------*/
/* asian.cpp */
/* --------------------------------------------------------------------*/
short arithmetic_asian(double Fwd_Avg, double Past_Avg, double X, double vol,
                       double r, double q, double T, double t, double start,
                       double stop, OPTTYPE call_put, SENSTYPE sens,
                       double *result);

/* --------------------------------------------------------------------*/
/* basket.cpp */
/* --------------------------------------------------------------------*/
short basket(double *S, double X, double *vol, double *w, double **corr,
             double r, double *q, double *t, double exp_t, int nCom,
             unsigned long nSim, OPTTYPE call_put, SENSTYPE sens,
             double *result);

/* --------------------------------------------------------------------*/
/* binary.cpp */
/* --------------------------------------------------------------------*/
short binary(double S, double X, double Q, double vol, double r, double q,
             double t, OPTTYPE call_put, SENSTYPE sens, double *result);

/* --------------------------------------------------------------------*/
/* boot.cpp */
/* --------------------------------------------------------------------*/
short BootstrapYields(double *yields, date_t *dates, int num_yields,
                      int interval, double **ppdDisc, date_t **ppdtDisc_dates,
                      int *num_disc);

/* --------------------------------------------------------------------*/
/* cmpd.cpp */
/* --------------------------------------------------------------------*/
short cmpd(double S, double X1, double X2, double vol,
           double r, double q, double t1, double t2,
           OPTTYPE opt_overlying, OPTTYPE opt_underlying,
           SENSTYPE sens, double *result);

/* --------------------------------------------------------------------*/
/* euro.cpp */
/* --------------------------------------------------------------------*/
short euro(double S, double X, double vol, double r, double q, double t,
           OPTTYPE call_put, SENSTYPE sens, double *result);
short euro_vol(double S, double X, double price, double r, double q, double t,
               OPTTYPE call_put, double *result);
short fwd_euro(double S, double vol, double r, double q, double t1, double t2,
               OPTTYPE call_put, SENSTYPE sens, double *result);

/* Do not call the following four functions, for use by cmpd.cpp only */
short euro_call(double S, double X, double vol, double r, double q, double t,
                double *result); /* Needed by cmpd.cpp */
short euro_put(double S, double X, double vol, double r, double q, double t,
               double *result); /* Needed by cmpd.cpp */
short delta_euro_call(double S, double X, double vol, double r, double q, double t,
                      double *result); /* Needed by cmpd.cpp */
short delta_euro_put(double S, double X, double vol, double r, double q, double t,
                     double *result); /* Needed by cmpd.cpp */

/* --------------------------------------------------------------------*/
/* excg.cpp */
/* --------------------------------------------------------------------*/
short excg(double S1, double S2, double vol1, double vol2,
           double rho, double q1, double q2, double t,
           OPTTYPE call_put, SENSTYPE sens,
           double *result);

/* --------------------------------------------------------------------*/
/* faure.cpp */
/* --------------------------------------------------------------------*/
short rvsf(double out[], int s, unsigned long MaxN, SEQTYPE sequence_type);
void rvsf_reset(void);

/* --------------------------------------------------------------------*/
/* optutil.cpp */
/* --------------------------------------------------------------------*/
double dN(double x);
double N(double x);
double Ninv(double u);
double M(double a, double b, double p);
void **newMatrix(unsigned long nrow, unsigned long ncol, size_t sElem);
void deleteMatrix(void *m);
double **Cholesky(double **a, int rows);

/* --------------------------------------------------------------------*/
/* spread.cpp */
/* --------------------------------------------------------------------*/
short spread(double S1, double S2, double X, double vol1, double vol2,
             double rho, double r, double q1, double q2, double t,
             OPTTYPE call_put, SENSTYPE sens, double *result);

#endif
