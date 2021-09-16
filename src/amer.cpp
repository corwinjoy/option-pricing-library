/*------------------------------------------------------------------------*/
/* American Option Pricer using the approximation of                      */
/* Barone-Adesi and Whaley                                                */
/*------------------------------------------------------------------------*/

#include "option.h"

/*------------------------------------------------------------------------*/
/* Equation for d1                                                        */
/* Used in the critical price equation below                              */
/*------------------------------------------------------------------------*/
static double fn_d1(double S, double X, double drift, double vol_rt_t)
{
    return (log(S / X) + drift) / vol_rt_t;
}

/*------------------------------------------------------------------------*/
/* Critical Call Price Equation                                           */
/* Equation for the critical price above which the American Call option   */
/* should be exercised                                                    */
/* This equation is solved below by using bisection                       */
/*------------------------------------------------------------------------*/
static double call_crit_price(double S_star, double X,
                              double g2, double drift,
                              double vol_rt_t, double disc_q, double disc_r)
{
    double c;
    double d1, d2, N_d1;

    d1 = fn_d1(S_star, X, drift, vol_rt_t);
    d2 = d1 - vol_rt_t;
    N_d1 = N(d1);
    c = S_star * disc_q * N_d1 - X * disc_r * N(d2);
    return S_star - X - (c + (1.0 - disc_q * N_d1) * S_star / g2);
}

/*------------------------------------------------------------------------*/
/* Critical Put Price Equation                                            */
/* Equation for the critical price above which the American Call option   */
/* should be exercised                                                    */
/* This equation is solved below by using bisection                       */
/*------------------------------------------------------------------------*/
static double put_crit_price(double S_star, double X,
                             double g1, double drift,
                             double vol_rt_t, double disc_q, double disc_r)
{
    double p;
    double d1, d2, N_neg_d1;

    d1 = fn_d1(S_star, X, drift, vol_rt_t);
    d2 = d1 - vol_rt_t;
    N_neg_d1 = N(-d1);
    p = X * disc_r * N(-d2) - S_star * disc_q * N_neg_d1;
    return X - S_star - (p - (1.0 - disc_q * N_neg_d1) * S_star / g1);
}

/*------------------------------------------------------------------------*/
/* Routine to find critical price (zero root) of either call_crit_price   */
/* or put_crit_price depending on the parameter fn                        */
/* Based on bisection                                                     */
/*------------------------------------------------------------------------*/
static short crit_bisect(double S, double X,
                         double drift, double vol_rt_t,
                         double disc_q, double disc_r, double g,
                         double *result,
                         double (*fn)(double S_star, double X,
                                      double g, double drift,
                                      double vol_rt_t, double disc_q, double disc_r))
{

    double min = S / 2.0;
    double max = S * 2.0;
    double mid, temp;
    double f_min, f_max, f_mid;

    static const int max_iter = 40;
    double acc = S * 0.001;
    int i;

    f_min = (*fn)(min, X, g, drift,
                  vol_rt_t, disc_q, disc_r);
    f_max = (*fn)(max, X, g, drift,
                  vol_rt_t, disc_q, disc_r);

    if (f_min * f_max >= 0.0)
    { /* No root in this interval */
        *result = -1.0;
        return FALSE;
    }

    if (f_min > 0.0)
    { /* Switch min & max s.t. min represents lower bound */
        temp = max;
        max = min;
        min = temp;
    }

    for (i = 0; i < max_iter; i++)
    {
        mid = (min + max) / 2.0;
        f_mid = (*fn)(mid, X, g, drift,
                      vol_rt_t, disc_q, disc_r);
        if (f_mid < 0.0)
            min = mid;
        else
            max = mid;
        if (fabs(max - min) < acc)
        { /* Check for convergence */
            *result = (min + max) / 2.0;
            return TRUE;
        }
    }

    /* Failed to converge */
    *result = (min + max) / 2.0;
    return FALSE;
}

/*------------------------------------------------------------------------*/
/* Function A1, used below for American Put Option Price - See Hull       */
/*------------------------------------------------------------------------*/
static double A1(double S_star, double X, double g1, double disc_q,
                 double drift, double vol_rt_t)
{
    double d1;
    double N_d1;

    d1 = fn_d1(S_star, X, drift, vol_rt_t);
    N_d1 = N(-d1);

    return -S_star / g1 * (1.0 - disc_q * N_d1);
}

/*------------------------------------------------------------------------*/
/* Function A2, used below for American Call Option Price - See Hull      */
/*------------------------------------------------------------------------*/
static double A2(double S_star, double X, double g2, double disc_q,
                 double drift, double vol_rt_t)
{
    double d1;

    d1 = fn_d1(S_star, X, drift, vol_rt_t);

    return S_star / g2 * (1.0 - disc_q * N(d1));
}

/*------------------------------------------------------------------------*/
/* American Call Option                                                   */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.367                                           */
/*------------------------------------------------------------------------*/
short amer_call(double S, double X, double vol, double r, double q, double t,
                double *result, double *crit_price)
{
    double S_star;
    short status;

    double drift = (r - q + vol * vol / 2) * t;
    double vol_rt_t;
    double disc_q = exp(-q * t);
    double disc_r = exp(-r * t);
    double h = 1.0 - disc_r;
    double a;
    double b;
    double g2;
    double c;
    double d1, d2, N_d1;

    if (vol < INPUT_UNDERFLOW)
    {
        *result = MAX(S * exp(-q * t) - X * exp(-r * t), 0.0);
        return TRUE;
    }

    if (t < INPUT_UNDERFLOW)
    {
        *result = MAX(S - X, 0.0);
        return TRUE;
    }

    a = 2.0 * r / (vol * vol);
    b = 2.0 * (r - q) / (vol * vol);
    g2 = (-(b - 1.0) + sqrt((b - 1.0) * (b - 1.0) + 4.0 * a / h)) / 2.0;
    vol_rt_t = vol * sqrt(t);

    if (*crit_price == 0.0)
    { /* Critical Price unknown, must calc */
        status = crit_bisect(S, X, drift, vol_rt_t, disc_q, disc_r,
                             g2, &S_star, call_crit_price);
        if (status == FALSE && S_star == -1.0)
        {
            /* No value for early exercise - return Euro Call */
            d1 = fn_d1(S, X, drift, vol_rt_t);
            d2 = d1 - vol_rt_t;
            N_d1 = N(d1);
            c = S * disc_q * N_d1 - X * disc_r * N(d2);
            *result = c;
            return TRUE;
        }
        if (status == FALSE)
        { /* Failed to converge */
            *result = 0.0;
            return FALSE;
        }
        *crit_price = S_star; /* Save critical price */
    }
    else
        S_star = *crit_price;

    if (S < S_star)
    { /* Early exercise not optimal */
        d1 = fn_d1(S, X, drift, vol_rt_t);
        d2 = d1 - vol_rt_t;
        N_d1 = N(d1);
        c = S * disc_q * N_d1 - X * disc_r * N(d2);
        *result = c + A2(S_star, X, g2, disc_q, drift, vol_rt_t) *
                          pow(S / S_star, g2);
        return TRUE;
    }
    else /* Early exercise is optimal */
    {
        *result = S - X;
        return TRUE;
    }
}

/*------------------------------------------------------------------------*/
/* American Put Option                                                    */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.367                                           */
/*------------------------------------------------------------------------*/
short amer_put(double S, double X, double vol, double r, double q, double t,
               double *result, double *crit_price)
{
    double S_star;
    short status;

    double drift = (r - q + vol * vol / 2) * t;
    double vol_rt_t;
    double disc_q = exp(-q * t);
    double disc_r = exp(-r * t);
    double h = 1.0 - disc_r;
    double a;
    double b;
    double g1;

    double p;
    double d1, d2;

    double A1_val, pow_val;

    if (vol < INPUT_UNDERFLOW)
    {
        *result = MAX(X * exp(-r * t) - S * exp(-q * t), 0.0);
        return TRUE;
    }

    if (t < INPUT_UNDERFLOW)
    {
        *result = MAX(X - S, 0.0);
        return TRUE;
    }

    a = 2.0 * r / (vol * vol);
    b = 2.0 * (r - q) / (vol * vol);
    g1 = (-(b - 1.0) - sqrt((b - 1.0) * (b - 1.0) + 4.0 * a / h)) / 2.0;
    vol_rt_t = vol * sqrt(t);

    if (*crit_price == 0.0)
    { /* Critical price unknown must calc */
        status = crit_bisect(S, X, drift, vol_rt_t, disc_q, disc_r,
                             g1, &S_star, put_crit_price);
        if (status == FALSE && S_star == -1.0)
        {
            /* No value for early exercise - return Euro Put */
            d1 = fn_d1(S, X, drift, vol_rt_t);
            d2 = d1 - vol_rt_t;
            p = X * exp(-r * t) * N(-d2) - S * exp(-q * t) * N(-d1);
            *result = p;
            return TRUE;
        }
        if (status == FALSE)
        { /* Failed to converge */
            *result = 0.0;
            return FALSE;
        }
        *crit_price = S_star; /* Save critical price */
    }
    else
        S_star = *crit_price;

    if (S > S_star)
    { /* Early exercise not optimal */
        d1 = fn_d1(S, X, drift, vol_rt_t);
        d2 = d1 - vol_rt_t;
        p = X * exp(-r * t) * N(-d2) - S * exp(-q * t) * N(-d1);
        A1_val = A1(S_star, X, g1, disc_q, drift, vol_rt_t);
        pow_val = pow(S / S_star, g1);
        *result = p + A1_val * pow_val;
        return TRUE;
    }
    else /* Early exercise is optimal */
    {
        *result = X - S;
        return TRUE;
    }
}

/*------------------------------------------------------------------------*/
/* Standard American Options                                              */
/* Wrapper to call the above routines                                     */
/*------------------------------------------------------------------------*/
short amer(double S, double X, double vol, double r, double q, double t,
           OPTTYPE call_put, SENSTYPE sens, double *result)
{
    double crit_price = 0.0;

    double dx, f1, f2, f3;
    short status;
    short (*amer_fn)(double S, double X, double vol, double r, double q, double t,
                     double *result, double *crit_price);

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
        amer_fn = amer_call;
    else
        amer_fn = amer_put;

    switch (sens)
    {
    case SENS_PRICE:
        return (*amer_fn)(S, X, vol, r, q, t, result, &crit_price);
    case SENS_DELTA:
        dx = S / (double)eps + 1.0 / 16.0 / (double)eps;
        status = (*amer_fn)(S, X, vol, r, q, t, &f1, &crit_price);
        status = status | (*amer_fn)(S + dx, X, vol, r, q, t, &f2, &crit_price);
        *result = (f2 - f1) / dx;
        return status;
    case SENS_GAMMA:
        dx = S / (double)eps + 1.0 / 16.0 / (double)eps;
        status = (*amer_fn)(S + dx, X, vol, r, q, t, &f1, &crit_price);
        status = status | (*amer_fn)(S - dx, X, vol, r, q, t, &f2, &crit_price);
        status = status | (*amer_fn)(S, X, vol, r, q, t, &f3, &crit_price);
        *result = (f1 + f2 - 2.0 * f3) / dx / dx;
        return status;
    case SENS_VEGA:
        dx = vol / (double)eps + 1.0 / 16.0 / (double)eps;
        status = (*amer_fn)(S, X, vol, r, q, t, &f1, &crit_price);
        status = status | (*amer_fn)(S, X, vol + dx, r, q, t, &f2, &crit_price);
        *result = (f2 - f1) / dx;
        return status;
    case SENS_THETA:
        dx = t / (double)eps + 1.0 / 16.0 / (double)eps;
        status = (*amer_fn)(S, X, vol, r, q, t, &f1, &crit_price);
        status = status | (*amer_fn)(S, X, vol, r, q, t - dx, &f2, &crit_price);
        *result = (f2 - f1) / dx;
        return status;
    case SENS_RHO:
        dx = r / (double)eps + 1.0 / 16.0 / (double)eps;
        status = (*amer_fn)(S, X, vol, r, q, t, &f1, &crit_price);
        if (r == q)
            status = status | (*amer_fn)(S, X, vol, r + dx, q + dx, t, &f2, &crit_price);
        else
            status = status | (*amer_fn)(S, X, vol, r + dx, q, t, &f2, &crit_price);
        *result = (f2 - f1) / dx;
        return status;
    default: /* Invalid Sensitivity */
        *result = 0.0;
        return FALSE;
    }
}

/*------------------------------------------------------------------------*/
/* Standard American Implied Vol Calculator                               */
/* Implied vol is computed by using Newton-Raphson on the Barone-Adesi    */
/* formula.  See "Numerical Recipes in C, 2nd Ed., p. 364" for details    */
/*------------------------------------------------------------------------*/
short amer_vol(double S, double X, double price, double r, double q, double t,
               OPTTYPE call_put, double *result)
{

    static const int max_iter = 20;
    static const double acc = 0.0001;
    short (*fn)(double S, double X, double vol, double r, double q, double t,
                double *result, double *crit_price);
    double guess, f, f2, df, dx, crit_price;
    int i;
    short status;

    /* Check Input Parameters */
    assert(result != NULL);
    if (call_put != OPT_CALL && call_put != OPT_PUT)
    {
        *result = 0.0;
        return FALSE;
    }
    if (S <= 0 || X <= 0 || r < 0 || q < 0)
    {
        *result = 0.0;
        return FALSE;
    }

    if (call_put == OPT_CALL)
    {
        fn = amer_call;
    }
    else
    {
        fn = amer_put;
    }

    guess = 0.3; /* Set initial vol guess at 30% */
    for (i = 0; i < max_iter; i++)
    {
        crit_price = 0.0;
        dx = guess / (double)eps;
        status = (*fn)(S, X, guess, r, q, t, &f, &crit_price); /* amer price under our vol guess */
        status = status | (*fn)(S, X, guess + dx, r, q, t, &f2, &crit_price);
        if (status == FALSE)
        {
            *result = 0.0;
            return FALSE;
        }
        df = (f2 - f) / dx; /* Amer derivative wrt volatility */
        f -= price;			/*Newton's method finds root at zero */
        dx = f / df;
        guess -= dx;
        if (guess <= 0 || guess > 8)
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
