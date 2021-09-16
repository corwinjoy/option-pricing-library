#include "option.h"

/*---------------------------------------------------------*/
/* Function to Compute Arithmetic Asian Call Value         */
/* This value is approximated using Levy's                 */
/* approach "Pricing European Average Rate Options,        */
/* Journal of International Money and Finance, 1992        */
/* Volume 11, p. 474-491                                   */
/* Checks were made vs ["Average Intelligence," RISK 5,2   */
/* (1992) p.53-58]. And Turnbull & Wakeman [JFQA, 26       */
/* (Sept. 1991), 377-89]                                   */
/* Generalizations of the above formulas for in period     */
/* averages can be found in "Garman on Exotic Options,     */
/* Risk Magazine conference notes, Garman: Asian-89.       */
/* Or see Hull p. 423. Also see Garman on Exotic Options   */
/* for in progress averaging Asain-31, Asian-89  and       */
/* Levy ["Pricing European average rate currency options," */
/* Journal of International Money and Finance (1992), 11   */
/* p. 489 ]                                                */
/*---------------------------------------------------------*/

static short arithmetic_asian_setup(
    /* This subroutine computes a modified strike and a modified */
    /* forward price for use in the Black-Scholes equation to    */
    /* value asian options - the general idea is to choose these */
    /* two variables by assuming that the average is lognormally */
    /* distributed and match the first two moments.  See Levy    */
    /* for details                                               */
    
    /* Input Parameters For Arithmetic Asian Option Valuation    */
    double Fwd_Avg,	 /* Current estimate for the average from the
                                start of the average period until the end of
                                the averaging period */
    double Past_Avg, /* Estimate of in-progress average from start until now (t) */
    double *X,		 /* Strike price for the average */
    double r,		 /* Domestic risk free rate */
    double q,		 /* Dividend Yield/ Foreign Yield  */
    double vol,		 /* Volatility */
    double T,		 /* Option Delivery Date */
    double t,		 /* Today */
    double start,	 /* Start Date for averaging */
    double stop,	 /* Stop Date for averaging */
    double *Sa,		 /* Sa in Levy Expression */
    double *d1,
    double *d2,
    double *beta, /* scaling factor for in progress averaging */
    OPTTYPE call_put,
    /* Output Parameters For Option Valuation */
    double *value)
{

    double alpha, oldX, g, D, Y, V, var;

    if (stop <= start || stop > T)
    {
        *value = 0.0;
        return FALSE;
    }

    if (stop <= t || T < INPUT_UNDERFLOW)
    { /* Option has finished averaging */
        if (call_put == OPT_CALL)
        {
            if (T > t) /* Still have to wait to be paid */
                *value = exp(-r * (T - t)) * MAX(Past_Avg - (*X), 0.0);
            else
                *value = MAX(Past_Avg - (*X), 0.0);
            return TRUE;
        }
        else
        {			   /* Put option */
            if (T > t) /* Still have to wait to be paid */
                *value = exp(-r * (T - t)) * MAX((*X) - Past_Avg, 0.0);
            else
                *value = MAX((*X) - Past_Avg, 0.0);
            return TRUE;
        }
    }

    if (start < t)
    { /* Averaging has already started, make
                       adjustment to strike price. C.F. Garman: Asian-31
                       or Turnbull & Wakeman p. 382
                    */
        alpha = (t - start) / (stop - start);
        *beta = 1.0 - alpha;
        oldX = *X;
        *X = ((*X) - alpha * Past_Avg) / (*beta);
        if ((*X) < 0.0)
        { /*This option is deeply in the money,
                         just return value of the average minus
                   the strike price discounted back to today
                  */
            if (call_put == OPT_CALL)
            {
                *value = exp(-r * (T - t)) * MAX(alpha * Past_Avg + (*beta) * Fwd_Avg - oldX, 0);
                return TRUE;
            }
            else
            {
                *value = exp(-r * (T - t)) * MAX(oldX - (alpha * Past_Avg + (*beta) * Fwd_Avg), 0);
                return TRUE;
            }
        }
        start = t;
    }
    else
    {
        alpha = 0.0;
        *beta = 1.0;
    }

    /* Levy's Approximation for an arithmetic asian option as
       extended by Garman */
    g = r - q;
    if (g < 0.0001)
        g = 0.0001;

    /* Set today=time zero */
    T -= t;
    start -= t;
    stop -= t;

    var = vol * vol;

    if (vol <= 0.01)
    { /* Handle zero vol case */
        if (call_put == OPT_CALL)
        {
            *value = exp(-r * (T - t)) * MAX(alpha * Past_Avg + (*beta) * Fwd_Avg - oldX, 0);
            return TRUE;
        }
        else
        {
            *value = exp(-r * (T - t)) * MAX(oldX - (alpha * Past_Avg + (*beta) * Fwd_Avg), 0);
            return TRUE;
        }
    }

    if (start == 0.0)
    { /* Compute Sa=E[M(t)] and D=E[M(t)^2] via <7b> and <8b> on
        p. 480 of Levy's "Pricing European currency options "
      */
        *Sa = exp(-r * stop) * Fwd_Avg / stop / g * (exp(g * stop) - 1.);
        Y = 2. * Fwd_Avg * Fwd_Avg / (g + var);
        Y *= ((exp((2. * g + var) * stop) - 1.) / (2. * g + var) -
              (exp(g * stop) - 1.) / g);
        D = Y / stop / stop;
    }
    else
    { /* Compute Sa=E[M(t)] and D=E[M(t)^2] via <7b> and <8b> on
        p. 480 of Levy's "Pricing European currency options "
      */

        *Sa = exp(-r * stop) * Fwd_Avg * exp(g * start) / (stop - start) / g * (exp(g * (stop - start)) - 1.);
        Y = 2. * Fwd_Avg * Fwd_Avg * exp((2. * g + var) * start) / (g + var) *
            ((exp((2. * g + var) * (stop - start)) - 1.) / (2. * g + var) -
             (exp(g * (stop - start)) - 1.) / g);
        D = Y / (stop - start) / (stop - start);
    }

    V = log(D) - 2 * (r * stop + log(*Sa)); /* This is v^2, eqn. 7b from "Average Intelligence" */

    if (V <= 0.0)
        V = 0.000001;

    *d1 = (0.5 * log(D) - log(*X)) / sqrt(V); /* alpha + v^2 = 0.5 E[M(t)^2] = 0.5 ln(D) */
                                              /* See (6),(7a), (7b) in "Avg. Intelligence" */
    *d2 = *d1 - sqrt(V);

    return TRUE;
}

/*---------------------------------------------------------*/
/* Function to Compute Arithmetic Asian Call Value         */
/*---------------------------------------------------------*/
static short arithmetic_asian_call(double Fwd_Avg, double Past_Avg, double X,
                                   double r, double q, double vol, double T, double t, double start,
                                   double stop, double *value)
{

    double Sa;
    double d1;
    double d2;
    double beta;

    short success;

    *value = -1.0;
    success = arithmetic_asian_setup(Fwd_Avg, Past_Avg, &X, r, q, vol,
                                     T, t, start, stop, &Sa,
                                     &d1, &d2, &beta, OPT_CALL, value);

    if (success == FALSE) /* Failure */
        return success;

    if (*value != -1.0) /* Option deep in money - return value from setup */
        return TRUE;

    *value = Sa * N(d1) - exp(-r * stop) * X * N(d2);
    *value *= beta; /* Account for adjusted strike if in averaging period */

    if (T > stop)
        (*value) *= exp(-r * (T - stop)); /* Discount if payment after settle */

    return TRUE;
}

/*---------------------------------------------------------*/
/* Function to Compute Arithmetic Asian Put Value          */
/* Uses put-call parity from "Average Intelligence"        */
/*---------------------------------------------------------*/
static short arithmetic_asian_put(double Fwd_Avg, double Past_Avg, double X,
                                  double r, double q, double vol, double T, double t, double start,
                                  double stop, double *value)
{

    double Sa;
    double d1;
    double d2;
    double beta;

    short success;

    *value = -1.0;
    success = arithmetic_asian_setup(Fwd_Avg, Past_Avg, &X, r, q, vol,
                                     T, t, start, stop, &Sa,
                                     &d1, &d2, &beta, OPT_PUT, value);

    if (success == FALSE) /* Failure */
        return success;

    if (*value != -1.0) /* Option deep in money - return value from setup */
        return TRUE;

    *value = Sa * N(d1) - exp(-r * stop) * X * N(d2);

    /* Value is now = value for asian call; apply put-call parity */
    *value -= Sa - exp(-r * stop) * X;

    *value *= beta; /* Account for adjusted strike if in averaging period */

    if (T > stop)
        (*value) *= exp(-r * (T - stop)); /* Discount if payment after settle */

    return TRUE;
}

/*------------------------------------------------------------------------*/
/* Asian Options                                                          */
/* Wrapper to call the above routines                                     */
/*------------------------------------------------------------------------*/
short arithmetic_asian(double Fwd_Avg, double Past_Avg, double X, double vol,
                       double r, double q, double T, double t, double start,
                       double stop, OPTTYPE call_put, SENSTYPE sens,
                       double *result)
{

    double dx, f1, f2, f3;
    short status;
    short (*val_fn)(double Fwd_Avg, double Past_Avg, double X,
                    double r, double q, double vol, double T,
                    double t, double start,
                    double stop, double *value);

    /* Check Input Parameters */
    if (call_put != OPT_CALL && call_put != OPT_PUT)
    {
        *result = 0.0;
        return FALSE;
    }

    /* Check Inputs */
    assert(result != NULL);
    if (Fwd_Avg <= 0 || X <= 0 || vol < 0 || r < 0 || q < 0 || (stop - start) < 0 || stop > T)
    {
        *result = 0.0;
        return FALSE;
    }
    if ((Fwd_Avg / X) < 0.1)
    {
        *result = 0.0;
        return FALSE;
    }
    if ((T - t) < 0.0001)
        if (stop == T)
        {
            T = t;
            stop = t;
        }
        else
            T = t;

    if ((stop - start) < 0.0001)
        start = stop - 0.0001;

    if (call_put == OPT_CALL)
        val_fn = arithmetic_asian_call;
    else
        val_fn = arithmetic_asian_put;

    switch (sens)
    {
    case SENS_PRICE:
        return (*val_fn)(Fwd_Avg, Past_Avg, X, r, q, vol, T, t, start,
                         stop, result);
    case SENS_DELTA:
        dx = Fwd_Avg / (double)eps;
        status = (*val_fn)(Fwd_Avg, Past_Avg, X, r, q, vol, T, t, start,
                           stop, &f1);
        status = status | (*val_fn)(Fwd_Avg + dx, Past_Avg, X, r, q, vol, T, t,
                                    start, stop, &f2);
        *result = (f2 - f1) / dx;
        return status;
    case SENS_GAMMA:
        dx = Fwd_Avg / (double)eps;
        status = (*val_fn)(Fwd_Avg + dx, Past_Avg, X, r, q, vol, T, t,
                           start, stop, &f1);
        status = status | (*val_fn)(Fwd_Avg - dx, Past_Avg, X, r, q, vol, T,
                                    t, start, stop, &f2);
        status = status | (*val_fn)(Fwd_Avg, Past_Avg, X, r, q, vol, T,
                                    t, start, stop, &f3);
        *result = (f1 + f2 - 2.0 * f3) / dx / dx;
        return status;
    case SENS_VEGA:
        dx = vol / (double)eps + 1.0 / 16.0 / (double)eps;
        status = (*val_fn)(Fwd_Avg, Past_Avg, X, r, q, vol, T, t, start,
                           stop, &f1);
        status = status | (*val_fn)(Fwd_Avg, Past_Avg, X, r, q, vol + dx, T, t, start,
                                    stop, &f2);
        *result = (f2 - f1) / dx;
        return status;
    case SENS_THETA:
        dx = T / (double)eps + 1.0 / 16.0 / (double)eps;
        status = (*val_fn)(Fwd_Avg, Past_Avg, X, r, q, vol, T, t, start,
                           stop, &f1);
        status = status | (*val_fn)(Fwd_Avg, Past_Avg, X, r, q, vol, T - dx, t, start - dx,
                                    stop - dx, &f2);
        *result = (f2 - f1) / dx;
        return status;
    case SENS_RHO:
        dx = r / (double)eps + 1.0 / 16.0 / (double)eps;
        status = (*val_fn)(Fwd_Avg, Past_Avg, X, r, q, vol, T, t, start,
                           stop, &f1);
        if (r == q)
            status = status | (*val_fn)(Fwd_Avg, Past_Avg, X, r + dx, q + dx,
                                        vol, T, t, start, stop, &f2);
        else
            status = status | (*val_fn)(Fwd_Avg, Past_Avg, X, r + dx, q,
                                        vol, T, t, start, stop, &f2);
        *result = (f2 - f1) / dx;
        return status;
    default: /* Invalid sensitivity */
        *result = 0.0;
        return FALSE;
    }
}
