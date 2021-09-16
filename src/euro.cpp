#include "option.h"

/*------------------------------------------------------------------------*/
/* Black - Scholes European Call Option                                   */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.248                                           */
/*------------------------------------------------------------------------*/
short euro_call(double S, double X, double vol, double r, double q, double t,
                double *result)
{
    double c, d1, d2, vol_rt_t;

    if (t <= INPUT_UNDERFLOW)
    {
        *result = MAX(S - X, 0.0);
        return TRUE;
    }

    if (vol < INPUT_UNDERFLOW)
    {
        *result = MAX(S * exp(-q * t) - X * exp(-r * t), 0.0);
        return TRUE;
    }

    vol_rt_t = vol * sqrt(t);
    d1 = (log(S / X) + (r - q + vol * vol / 2) * t) / vol_rt_t;
    d2 = d1 - vol_rt_t;
    c = S * exp(-q * t) * N(d1) - X * exp(-r * t) * N(d2);
    *result = c;

    return TRUE;
}

/*------------------------------------------------------------------------*/
/* Black - Scholes European Put Option                                    */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.248                                           */
/*------------------------------------------------------------------------*/
short euro_put(double S, double X, double vol, double r, double q, double t,
               double *result)
{
    double p, d1, d2, vol_rt_t;

    if (t <= INPUT_UNDERFLOW)
    {
        *result = MAX(X - S, 0.0);
        return TRUE;
    }

    if (vol < INPUT_UNDERFLOW)
    {
        *result = MAX(X * exp(-q * t) - S * exp(-r * t), 0.0);
        return TRUE;
    }

    vol_rt_t = vol * sqrt(t);
    d1 = (log(S / X) + (r - q + vol * vol / 2) * t) / vol_rt_t;
    d2 = d1 - vol_rt_t;
    p = X * exp(-r * t) * N(-d2) - S * exp(-q * t) * N(-d1);
    *result = p;

    return TRUE;
}

/*------------------------------------------------------------------------*/
/* Black - Scholes European Call Option Delta                             */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.305                                           */
/*------------------------------------------------------------------------*/
short delta_euro_call(double S, double X, double vol, double r, double q, double t,
                      double *result)
{
    double d1, vol_rt_t;

    if (t <= INPUT_UNDERFLOW)
    {
        if (S > X)
            *result = 1.0;
        else
            *result = 0.0;
        return TRUE;
    }
    if (vol < INPUT_UNDERFLOW)
    {
        if (S > X)
            *result = exp(-q * t);
        else
            *result = 0.0;
        return TRUE;
    }
    vol_rt_t = vol * sqrt(t);
    d1 = (log(S / X) + (r - q + vol * vol / 2) * t) / vol_rt_t;
    *result = exp(-q * t) * N(d1);

    return TRUE;
}

/*------------------------------------------------------------------------*/
/* Black - Scholes European Put Option  Delta                             */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.305                                           */
/*------------------------------------------------------------------------*/
short delta_euro_put(double S, double X, double vol, double r, double q, double t,
                     double *result)
{
    double d1, vol_rt_t;

    if (t <= INPUT_UNDERFLOW)
    {
        if (S < X)
            *result = 1.0;
        else
            *result = 0.0;
        return TRUE;
    }
    if (vol < INPUT_UNDERFLOW)
    {
        if (S < X)
            *result = -exp(-q * t);
        else
            *result = 0.0;
        return TRUE;
    }

    vol_rt_t = vol * sqrt(t);
    d1 = (log(S / X) + (r - q + vol * vol / 2) * t) / vol_rt_t;
    *result = exp(-q * t) * (N(d1) - 1.0);

    return TRUE;
}

/*------------------------------------------------------------------------*/
/* Black - Scholes European Call Option Theta                             */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.308                                           */
/*------------------------------------------------------------------------*/
static short theta_euro_call(double S, double X, double vol, double r, double q, double t,
                             double *result)
{
    double theta, d1, d2, vol_rt_t;

    if (t <= INPUT_UNDERFLOW)
    {
        *result = 0.0;
        return TRUE;
    }
    if (vol < INPUT_UNDERFLOW)
    {
        if (S > X)
            *result = q * S * exp(-q * t) - r * X * exp(-r * t);
        else
            *result = 0.0;
        return TRUE;
    }

    vol_rt_t = vol * sqrt(t);
    d1 = (log(S / X) + (r - q + vol * vol / 2) * t) / vol_rt_t;
    d2 = d1 - vol_rt_t;
    theta = -(S * dN(d1) * vol * exp(-q * t)) / (2.0 * sqrt(t));
    theta += q * S * N(d1) * exp(-q * t) - r * X * exp(-r * t) * N(d2);
    *result = theta;

    return TRUE;
}

/*------------------------------------------------------------------------*/
/* Black - Scholes European Put Option Theta                              */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.308                                           */
/*------------------------------------------------------------------------*/
static short theta_euro_put(double S, double X, double vol, double r, double q, double t,
                            double *result)
{
    double theta, d1, d2, vol_rt_t;

    if (t <= INPUT_UNDERFLOW)
    {
        *result = 0.0;
        return TRUE;
    }
    if (vol < INPUT_UNDERFLOW)
    {
        if (S < X)
            *result = -q * S * exp(-q * t) + r * X * exp(-r * t);
        else
            *result = 0.0;
        return TRUE;
    }

    vol_rt_t = vol * sqrt(t);
    d1 = (log(S / X) + (r - q + vol * vol / 2) * t) / vol_rt_t;
    d2 = d1 - vol_rt_t;
    theta = -(S * dN(d1) * vol * exp(-q * t)) / (2.0 * sqrt(t));
    theta += -q * S * N(-d1) * exp(-q * t) + r * X * exp(-r * t) * N(-d2);
    *result = theta;

    return TRUE;
}

/*------------------------------------------------------------------------*/
/* Black - Scholes European Call Option Gamma                             */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.313                                           */
/*------------------------------------------------------------------------*/
static short gamma_euro_call(double S, double X, double vol, double r, double q, double t,
                             double *result)
{
    double gamma, d1, vol_rt_t;

    if (vol < INPUT_UNDERFLOW || t <= INPUT_UNDERFLOW)
    {
        *result = 0.0;
        return TRUE;
    }

    vol_rt_t = vol * sqrt(t);
    d1 = (log(S / X) + (r - q + vol * vol / 2) * t) / vol_rt_t;
    gamma = dN(d1) * exp(-q * t) / (S * vol_rt_t);
    *result = gamma;

    return TRUE;
}

/*------------------------------------------------------------------------*/
/* Black - Scholes European Put Option Gamma                              */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.313                                           */
/*------------------------------------------------------------------------*/
static short gamma_euro_put(double S, double X, double vol, double r, double q, double t,
                            double *result)
{

    return gamma_euro_call(S, X, vol, r, q, t, result);
}

/*------------------------------------------------------------------------*/
/* Black - Scholes European Call Option Vega                              */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.316                                           */
/*------------------------------------------------------------------------*/
static short vega_euro_call(double S, double X, double vol, double r, double q, double t,
                            double *result)
{
    double vega, d1, vol_rt_t;

    if (vol < INPUT_UNDERFLOW || t <= INPUT_UNDERFLOW)
    {
        *result = 0.0;
        return TRUE;
    }

    vol_rt_t = vol * sqrt(t);
    d1 = (log(S / X) + (r - q + vol * vol / 2) * t) / vol_rt_t;
    vega = S * sqrt(t) * dN(d1) * exp(-q * t);
    *result = vega;

    return TRUE;
}

/*------------------------------------------------------------------------*/
/* Black - Scholes European Put Option Vega                               */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.316                                           */
/*------------------------------------------------------------------------*/
static short vega_euro_put(double S, double X, double vol, double r, double q, double t,
                           double *result)
{
    return vega_euro_call(S, X, vol, r, q, t, result);
}

/*------------------------------------------------------------------------*/
/* Black - Scholes European Call Option Rho                               */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.317                                           */
/* *AND*                                                                  */
/* From Natenberg "Option Volatility & Pricing Strategies,"               */
/* Probus Publishing, Chicago Ill. p.336                                  */
/* Note special case for when q = r by direct differentiation             */
/*------------------------------------------------------------------------*/
static short rho_euro_call(double S, double X, double vol, double r, double q, double t,
                           double *result)
{
    double rho, d1, d2, vol_rt_t;

    if (t <= INPUT_UNDERFLOW)
    {
        *result = 0.0;
        return TRUE;
    }
    if (vol < INPUT_UNDERFLOW)
    {
        if (S > X)
            if (q != r)
                *result = t * X * exp(-r * t);
            else
                *result = t * exp(-r * t) * (X - S);
        else
            *result = 0.0;
        return TRUE;
    }

    vol_rt_t = vol * sqrt(t);
    d1 = (log(S / X) + (r - q + vol * vol / 2) * t) / vol_rt_t;
    d2 = d1 - vol_rt_t;
    if (q != r)
        rho = t * X * exp(-r * t) * N(d2);
    else
        rho = t * exp(-r * t) * (X * N(d2) - S * N(d1));

    *result = rho;

    return TRUE;
}

/*------------------------------------------------------------------------*/
/* Black - Scholes European Put Option                                    */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.248                                           */
/*------------------------------------------------------------------------*/
static short rho_euro_put(double S, double X, double vol, double r, double q, double t,
                          double *result)
{
    double rho, d1, d2, vol_rt_t;

    if (t <= INPUT_UNDERFLOW)
    {
        *result = 0.0;
        return TRUE;
    }
    if (vol < INPUT_UNDERFLOW)
    {
        if (S < X)
            if (q != r)
                *result = -t * X * exp(-r * t);
            else
                *result = -t * exp(-r * t) * (X - S);
        else
            *result = 0.0;
        return TRUE;
    }

    vol_rt_t = vol * sqrt(t);
    d1 = (log(S / X) + (r - q + vol * vol / 2) * t) / vol_rt_t;
    d2 = d1 - vol_rt_t;
    if (q != r)
        rho = -t * X * exp(-r * t) * N(-d2);
    else
        rho = -t * exp(-r * t) * (X * N(-d2) - S * N(-d1));

    *result = rho;

    return TRUE;
}

/*------------------------------------------------------------------------*/
/* Black - Scholes European Options                                       */
/* Wrapper to call the above routines                                     */
/*------------------------------------------------------------------------*/
short euro(double S, double X, double vol, double r, double q, double t,
           OPTTYPE call_put, SENSTYPE sens, double *result)
{

    short status;

    /* Check Input Parameters */
    if (call_put != OPT_CALL && call_put != OPT_PUT)
    {
        *result = 0.0;
        return FALSE;
    }

    /* Check Inputs */
    assert(result != NULL);
    if (S <= INPUT_UNDERFLOW || X <= INPUT_UNDERFLOW || vol < INPUT_UNDERFLOW || r < 0 || q < 0)
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
    {
        switch (sens)
        {
        case SENS_PRICE:
            return euro_call(S, X, vol, r, q, t, result);
        case SENS_DELTA:
            return delta_euro_call(S, X, vol, r, q, t, result);
        case SENS_GAMMA:
            return gamma_euro_call(S, X, vol, r, q, t, result);
        case SENS_VEGA:
            return vega_euro_call(S, X, vol, r, q, t, result);
        case SENS_THETA:
            return theta_euro_call(S, X, vol, r, q, t, result);
        case SENS_RHO:
            return rho_euro_call(S, X, vol, r, q, t, result);
        case SENS_ALLSENS:
            status = euro_call(S, X, vol, r, q, t, &result[SENS_PRICE]);
            status = status && delta_euro_call(S, X, vol, r, q, t, &result[SENS_DELTA]);
            status = status && gamma_euro_call(S, X, vol, r, q, t, &result[SENS_GAMMA]);
            status = status && vega_euro_call(S, X, vol, r, q, t, &result[SENS_VEGA]);
            status = status && theta_euro_call(S, X, vol, r, q, t, &result[SENS_THETA]);
            status = status && rho_euro_call(S, X, vol, r, q, t, &result[SENS_RHO]);
            return status;

        default: /* Invalid Sensitivity */
            *result = 0.0;
            return FALSE;
        }
    }
    else
    {
        switch (sens)
        {
        case SENS_PRICE:
            return euro_put(S, X, vol, r, q, t, result);
        case SENS_DELTA:
            return delta_euro_put(S, X, vol, r, q, t, result);
        case SENS_GAMMA:
            return gamma_euro_put(S, X, vol, r, q, t, result);
        case SENS_VEGA:
            return vega_euro_put(S, X, vol, r, q, t, result);
        case SENS_THETA:
            return theta_euro_put(S, X, vol, r, q, t, result);
        case SENS_RHO:
            return rho_euro_put(S, X, vol, r, q, t, result);
        case SENS_ALLSENS:
            status = euro_put(S, X, vol, r, q, t, &result[SENS_PRICE]);
            status = status && delta_euro_put(S, X, vol, r, q, t, &result[SENS_DELTA]);
            status = status && gamma_euro_put(S, X, vol, r, q, t, &result[SENS_GAMMA]);
            status = status && vega_euro_put(S, X, vol, r, q, t, &result[SENS_VEGA]);
            status = status && theta_euro_put(S, X, vol, r, q, t, &result[SENS_THETA]);
            status = status && rho_euro_put(S, X, vol, r, q, t, &result[SENS_RHO]);
            return status;

        default: /* Invalid Sensitivity */
            *result = 0.0;
            return FALSE;
        }
    }
}

/*------------------------------------------------------------------------*/
/* Black - Scholes European Implied Vol Calculator                        */
/* Implied vol is computed by using Newton-Raphson on the Black-Scholes   */
/* formula.  See "Numerical Recipes in C, 2nd Ed., p. 364" for details    */
/*------------------------------------------------------------------------*/
short euro_vol(double S, double X, double price, double r, double q, double t,
               OPTTYPE call_put, double *result)
{

    static const int max_iter = 40;
    static const double acc = 0.0001;
    short (*fn)(double S, double X, double vol, double r, double q, double t,
                double *result);
    short (*dfn)(double S, double X, double vol, double r, double q, double t,
                 double *result);
    double guess, f, df, dx;
    int i;
    short status;

    /* Check Input Parameters */
    assert(result != NULL);
    if (call_put != OPT_CALL && call_put != OPT_PUT)
    {
        *result = 0.0;
        return FALSE;
    }
    if (S <= INPUT_UNDERFLOW || X <= INPUT_UNDERFLOW || r < 0 || q < 0 ||
        t <= INPUT_UNDERFLOW)
    {
        *result = 0.0;
        return FALSE;
    }

    if (call_put == OPT_CALL)
    {
        fn = euro_call;
        dfn = vega_euro_call;
    }
    else
    {
        fn = euro_put;
        dfn = vega_euro_put;
    }

    /* We choose an initial guess that should guarantee convergence */
    /* See Manaster, S. and G. Koehler, Journal of Finance, 37 (1), 1982 p. 227 */
    guess = sqrt(2.0 * fabs(log(S / X) + (r - q) * t) / t);
    if (guess == 0.0) /* i.e. S = X and r = q */
        guess = 0.3;  /* Try a reasonable value instead */

    for (i = 0; i < max_iter; i++)
    {
        status = (*fn)(S, X, guess, r, q, t, &f);			 /* BS price under our vol guess */
        f -= price;											 /*Newton's method finds root at zero */
        status = status | (*dfn)(S, X, guess, r, q, t, &df); /* BS derivative wrt volatility */
        if (status == FALSE || fabs(df) < acc * acc)
        {
            *result = 0.0;
            return FALSE;
        }
        dx = f / df;
        guess -= dx;
        if (guess <= acc || guess > 8)
        {
            *result = 0.0;
            return FALSE;
        }
        if (fabs(dx) < acc)
        { /* Convergence */
            *result = guess;
            return TRUE;
        }
    }

    /* Max # of iterations exceeded */
    *result = guess;
    return FALSE;
}

/*------------------------------------------------------------------------*/
/* Forward start European option, from Hull p. 416                        */
/*------------------------------------------------------------------------*/
short fwd_euro(double S, double vol, double r, double q, double t1, double t2,
               OPTTYPE call_put, SENSTYPE sens, double *result)
{
    double val;
    double X;
    short success;

    if (t1 < 0 || q < 0)
    {
        *result = 0.0;
        return FALSE;
    }
    else
    {
        X = S;
        success = euro(S, X, vol, r, q, t2 - t1, (OPTTYPE)call_put,
                       (SENSTYPE)sens, result);
        switch (sens)
        {
        case SENS_PRICE:
        case SENS_DELTA:
        case SENS_GAMMA:
        case SENS_VEGA:
            *result *= exp(-q * t1);
            break;
        case SENS_THETA:
            success = success | euro(S, X, vol, r, q, t2 - t1, (OPTTYPE)call_put,
                                     (SENSTYPE)SENS_PRICE, &val);
            *result = (*result) * exp(-q * t1) + q * exp(-q * t1) * val;
            break;
        case SENS_RHO:
            if (r == q)
            {
                success = success | euro(S, X, vol, r, q, t2 - t1, (OPTTYPE)call_put,
                                         (SENSTYPE)SENS_PRICE, &val);
                *result = (*result) * exp(-q * t1) - t1 * exp(-q * t1) * val;
            }
            else
                *result *= exp(-q * t1);
            break;
        }
    }
    return success;
}
