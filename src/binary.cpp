#include "option.h"

/*------------------------------------------------------------------------*/
/* Binary Cash or Nothing Call                                            */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.420                                           */
/*------------------------------------------------------------------------*/
static short binary_call(double S, double X, double Q, double vol, double r,
                         double q, double t, double *result)
{
    double c, d1, d2, vol_rt_t;

    if (vol < INPUT_UNDERFLOW)
    {
        if (S > X)
            *result = Q * exp(-r * t);
        else
            *result = 0.0;
        return TRUE;
    }

    if (t <= INPUT_UNDERFLOW)
    {
        if (S > X)
            *result = Q;
        else
            *result = 0.0;
        return TRUE;
    }

    vol_rt_t = vol * sqrt(t);
    d1 = (log(S / X) + (r - q + vol * vol / 2) * t) / vol_rt_t;
    d2 = d1 - vol_rt_t;
    c = Q * exp(-r * t) * N(d2);
    *result = c;

    return TRUE;
}

/*------------------------------------------------------------------------*/
/* Binary Cash or Nothing Put                                             */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.420                                           */
/*------------------------------------------------------------------------*/
static short binary_put(double S, double X, double Q, double vol, double r,
                        double q, double t, double *result)
{
    double p, d1, d2, vol_rt_t;

    if (vol < INPUT_UNDERFLOW)
    {
        if (S < X)
            *result = Q * exp(-r * t);
        else
            *result = 0.0;
        return TRUE;
    }

    if (t <= INPUT_UNDERFLOW)
    {
        if (S < X)
            *result = Q;
        else
            *result = 0.0;
        return TRUE;
    }
    vol_rt_t = vol * sqrt(t);
    d1 = (log(S / X) + (r - q + vol * vol / 2) * t) / vol_rt_t;
    d2 = d1 - vol_rt_t;
    p = Q * exp(-r * t) * N(-d2);
    *result = p;

    return TRUE;
}

/*------------------------------------------------------------------------*/
/* Binary Cash or Nothing Options                                         */
/* Wrapper to call the above routines                                     */
/*------------------------------------------------------------------------*/
short binary(double S, double X, double Q, double vol, double r, double q,
             double t, OPTTYPE call_put, SENSTYPE sens, double *result)
{
    double dx, f1, f2, f3;
    short status;
    short (*val_fn)(double S, double X, double Q, double vol, double r,
                    double q, double t, double *result);

    /* Check Input Parameters */
    if (call_put != OPT_CALL && call_put != OPT_PUT)
    {
        *result = 0.0;
        return FALSE;
    }

    /* Check Inputs */
    assert(result != NULL);
    if (S <= 0 || X <= 0 || vol < 0 || r < 0 || q < 0)
    {
        *result = 0.0;
        return FALSE;
    }
    if ((S / X) < 0.1)
    {
        *result = 0.0;
        return FALSE;
    }

    if (call_put == OPT_CALL)
        val_fn = binary_call;
    else
        val_fn = binary_put;

    switch (sens)
    {
    case SENS_PRICE:
        return (*val_fn)(S, X, Q, vol, r, q, t, result);
    case SENS_DELTA:
        dx = S / (double)eps;
        status = (*val_fn)(S, X, Q, vol, r, q, t, &f1);
        status = status | (*val_fn)(S + dx, X, Q, vol, r, q, t, &f2);
        *result = (f2 - f1) / dx;
        return status;
    case SENS_GAMMA:
        dx = S / (double)eps;
        status = (*val_fn)(S + dx, X, Q, vol, r, q, t, &f1);
        status = status | (*val_fn)(S - dx, X, Q, vol, r, q, t, &f2);
        status = status | (*val_fn)(S, X, Q, vol, r, q, t, &f3);
        *result = (f1 + f2 - 2.0 * f3) / dx / dx;
        return status;
    case SENS_VEGA:
        dx = vol / (double)eps + 1.0 / 16.0 / (double)eps;
        status = (*val_fn)(S, X, Q, vol, r, q, t, &f1);
        status = status | (*val_fn)(S, X, Q, vol + dx, r, q, t, &f2);
        *result = (f2 - f1) / dx;
        return status;
    case SENS_THETA:
        dx = t / (double)eps + 1.0 / 16.0 / (double)eps;
        status = (*val_fn)(S, X, Q, vol, r, q, t, &f1);
        status = status | (*val_fn)(S, X, Q, vol, r, q, t - dx, &f2);
        *result = (f2 - f1) / dx;
        return status;
    case SENS_RHO:
        dx = r / (double)eps + 1.0 / 16.0 / (double)eps;
        status = (*val_fn)(S, X, Q, vol, r, q, t, &f1);
        if (r == q)
            status = status | (*val_fn)(S, X, Q, vol, r + dx, q + dx, t, &f2);
        else
            status = status | (*val_fn)(S, X, Q, vol, r + dx, q, t, &f2);
        *result = (f2 - f1) / dx;
        return status;
    default: /* Invalid sensitivity */
        *result = 0.0;
        return FALSE;
    }
}
