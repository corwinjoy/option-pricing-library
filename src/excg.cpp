#include "option.h"

/*------------------------------------------------------------------------*/
/* Margrabe's Option to Exchange S2 for S1 = Max(S1-S2, 0)                */
/* Formula taken from "Options, Futures, and Other Derivative Securities, */
/* Hull, 2nd Ed., p. 423"                                                 */
/* See also "The Value of an Option to Exchange one Asset for Another,"   */
/* William Margrabe, Journal of Finance, Vol 33, #1, March 1978           */
/*------------------------------------------------------------------------*/
static short excg_call(double S1, double S2, double vol1,
                       double vol2, double rho,
                       double q1, double q2, double t,
                       double *result)
{
    double d1, d2, vol_sqr, vol_rt_t;

    vol_sqr = vol1 * vol1 + vol2 * vol2 - 2.0 * rho * vol1 * vol2;

    if (vol_sqr < INPUT_UNDERFLOW)
    {
        *result = MAX(S1 * exp(-q1 * t) - S2 * exp(-q2 * t), 0.0);
        return TRUE;
    }

    if (t <= INPUT_UNDERFLOW)
    {
        *result = MAX(S1 - S2, 0.0);
        return TRUE;
    }

    vol_rt_t = sqrt(vol_sqr * t);
    d1 = (log(S1 / S2) + (q2 - q1 + vol_sqr / 2.0) * t) / vol_rt_t;
    d2 = d1 - vol_rt_t;
    *result = S1 * exp(-q1 * t) * N(d1) - S2 * exp(-q2 * t) * N(d2);
    return TRUE;
}

/*------------------------------------------------------------------------*/
/* Margrabe's Option to Exchange S1 for S2 = Max(S2-S1, 0)                */
/* Formula taken from "Options, Futures, and Other Derivative Securities, */
/* Hull, 2nd Ed., p. 423"                                                 */
/* See also "The Value of an Option to Exchange one Asset for Another,"   */
/* William Margrabe, Journal of Finance, Vol 33, #1, March 1978           */
/*------------------------------------------------------------------------*/
static short excg_put(double S1, double S2, double vol1,
                      double vol2, double rho,
                      double q1, double q2, double t,
                      double *result)
{

    return excg_call(S2, S1, vol2, vol1, rho, q2, q1, t, result);
}

/*------------------------------------------------------------------------*/
/* Exchange Options                                                       */
/* Wrapper to call the above routines                                     */
/*------------------------------------------------------------------------*/
short excg(double S1, double S2, double vol1, double vol2,
           double rho, double q1, double q2, double t,
           OPTTYPE call_put, SENSTYPE sens,
           double *result)
{

    double dx, f1, f2, f3, f4, f5;
    short status;
    short (*val_fn)(double S1, double S2, double vol1,
                    double vol2, double rho,
                    double q1, double q2, double t,
                    double *result);

    /* Check Input Parameters */
    if (call_put != OPT_CALL && call_put != OPT_PUT)
    {
        *result = 0.0;
        return FALSE;
    }

    /* Check Inputs */
    assert(result != NULL);
    if (S1 <= 0 || S2 <= 0 || vol1 < 0 || vol2 < 0 || q1 < 0 || q2 < 0 || rho > 1 || rho < -1)
    {
        *result = 0.0;
        return FALSE;
    }
    if ((S1 / S2) < 0.1)
    {
        *result = 0.0;
        return FALSE;
    }

    if (call_put == OPT_CALL)
        val_fn = excg_call;
    else
        val_fn = excg_put;

    switch (sens)
    {
    case SENS_PRICE:
        return (*val_fn)(S1, S2, vol1, vol2, rho, q1, q2, t, result);
    case SENS_DELTA:
        dx = S1 / (double)eps;
        status = (*val_fn)(S1, S2, vol1, vol2, rho, q1, q2, t, &f1);
        status = status | (*val_fn)(S1 + dx, S2, vol1, vol2, rho, q1,
                                    q2, t, &f2);
        *result = (f2 - f1) / dx;
        return status;
    case SENS_GAMMA:
        /* Note use of 5 point method here for calculating Gamma rather
            than the usual 3 point method used in the other routines.
            We need to use a 5 point method here because 3 point Gamma
            is unstable as vol_rt_t --> d1. In particular the inputs
            S1	        S2     Vol1  Vol2  Corr  r      t
             $19.10 	 $19.00 	33%	33%	78%	4.91%	0.219028063
            give a  poor performance
         */
        dx = S1 / (double)eps;
        status = (*val_fn)(S1 + 2.0 * dx, S2, vol1, vol2, rho, q1, q2, t, &f1);
        status = status | (*val_fn)(S1 + dx, S2, vol1, vol2, rho, q1,
                                    q2, t, &f2);
        status = status | (*val_fn)(S1, S2, vol1, vol2, rho, q1,
                                    q2, t, &f3);
        status = status | (*val_fn)(S1 - dx, S2, vol1, vol2, rho, q1,
                                    q2, t, &f4);
        status = status | (*val_fn)(S1 - 2.0 * dx, S2, vol1, vol2, rho, q1,
                                    q2, t, &f5);
        *result = (((f1 - f2) - (f2 - f3)) + ((f2 - f3) - (f3 - f4)) + ((f3 - f4) - (f4 - f5))) / 3.0 / dx / dx;
        return status;
    case SENS_VEGA:
        dx = vol1 / (double)eps + 1.0 / 16.0 / (double)eps;
        status = (*val_fn)(S1, S2, vol1, vol2, rho, q1, q2, t, &f1);
        status = status | (*val_fn)(S1, S2, vol1 + dx, vol2, rho, q1, q2, t, &f2);
        *result = (f2 - f1) / dx;
        return status;
    case SENS_DELTA2:
        dx = S2 / (double)eps;
        status = (*val_fn)(S1, S2, vol1, vol2, rho, q1, q2, t, &f1);
        status = status | (*val_fn)(S1, S2 + dx, vol1, vol2, rho, q1,
                                    q2, t, &f2);
        *result = (f2 - f1) / dx;
        return status;
    case SENS_GAMMA2:
        dx = S2 / (double)eps;
        status = (*val_fn)(S1, S2 + dx, vol1, vol2, rho, q1, q2, t, &f1);
        status = status | (*val_fn)(S1, S2 - dx, vol1, vol2, rho, q1,
                                    q2, t, &f2);
        status = status | (*val_fn)(S1, S2, vol1, vol2, rho, q1,
                                    q2, t, &f3);
        *result = (f1 + f2 - 2.0 * f3) / dx / dx;
        return status;
    case SENS_VEGA2:
        dx = vol2 / (double)eps + 1.0 / 16.0 / (double)eps;
        status = (*val_fn)(S1, S2, vol1, vol2, rho, q1, q2, t, &f1);
        status = status | (*val_fn)(S1, S2, vol1, vol2 + dx, rho, q1, q2, t, &f2);
        *result = (f2 - f1) / dx;
        return status;
    case SENS_THETA:
        dx = t / (double)eps + 1.0 / 16.0 / (double)eps;
        status = (*val_fn)(S1, S2, vol1, vol2, rho, q1, q2, t, &f1);
        status = status | (*val_fn)(S1, S2, vol1, vol2, rho, q1, q2, t - dx, &f2);
        *result = (f2 - f1) / dx;
        return status;
    case SENS_RHO:
        if (q1 == q2)
        { /* Assume r = q1 = q2 for futures */
            dx = q1 / (double)eps + 1.0 / 16.0 / (double)eps;
            status = (*val_fn)(S1, S2, vol1, vol2, rho, q1, q2, t, &f1);
            status = status | (*val_fn)(S1, S2, vol1, vol2, rho,
                                        q1 + dx, q2 + dx, t, &f2);
            *result = (f2 - f1) / dx;
            return status;
        }
        else
        {
            *result = 0.0;
            return TRUE;
        }
    case SENS_CORR:
        dx = rho / (double)eps * 4.0 + 1.0 / 16.0 / (double)eps;
        status = (*val_fn)(S1, S2, vol1, vol2, rho, q1, q2, t, &f1);
        status = status | (*val_fn)(S1, S2, vol1, vol2, rho + dx, q1, q2, t, &f2);
        *result = (f2 - f1) / dx;
        return status;
    default: /* Invalid sensitivity */
        *result = 0.0;
        return FALSE;
    }
}
