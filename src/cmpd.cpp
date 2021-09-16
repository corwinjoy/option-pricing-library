#include "option.h"

/*------------------------------------------------------------------------*/
/* Black - Scholes European Implied Stock Price Calculator                */
/* Calculates the critical stock price S at which the option value is     */
/* equal to the amount given in price                                     */
/* Stock price is computed by using Newton-Raphson on the Black-Scholes   */
/* formula.  See "Numerical Recipes in C, 2nd Ed., p. 364" for details    */
/*------------------------------------------------------------------------*/
static short euro_S(double price, double X, double vol, double r, double q, double t,
                    OPTTYPE call_put, double *result)
{

    static const int max_iter = 20;
    static const double acc = 0.0001;
    double price_acc;
    short (*fn)(double S, double X, double vol, double r, double q, double t,
                double *result);
    short (*dfn)(double S, double X, double vol, double r, double q, double t,
                 double *result);
    double guess, f, df, dx;
    int i;

    /* Check Input Parameters */
    assert(result != NULL);
    if (call_put != OPT_CALL && call_put != OPT_PUT)
    {
        *result = 0.0;
        return FALSE;
    }
    if (price <= 0 || X <= 0 || r < 0 || q < 0)
    {
        *result = 0.0;
        return FALSE;
    }

    if (call_put == OPT_CALL)
    {
        fn = euro_call;
        dfn = delta_euro_call;
    }
    else
    {
        fn = euro_put;
        dfn = delta_euro_put;
    }

    guess = X; /* Set initial price = to exercise price */
    price_acc = acc * guess;

    for (i = 0; i < max_iter; i++)
    {
        (*fn)(guess, X, vol, r, q, t, &f);	 /* BS price under our guess for S*/
        f -= price;							 /*Newton's method finds root at zero */
        (*dfn)(guess, X, vol, r, q, t, &df); /* BS derivative wrt S */
        dx = f / df;
        guess -= dx;
        if (guess <= 0 || guess > 4.0 * X)
        {
            *result = 0.0;
            return FALSE;
        }
        if (fabs(dx) < price_acc)
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
/* Compound Option: Setup                                                 */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.417                                           */
/*------------------------------------------------------------------------*/
short cmpd_setup(double S, double X1, double X2, double vol,
                 double r, double q, double t1, double t2,
                 double *mu, double *var1, double *var2, double *a1,
                 double *a2, double *b1, double *b2, double *sqrt_t1_t2,
                 double *S_star, OPTTYPE underlying_type)
{

    short status;

    /* S_star is the stock price at time T1 for which the option price at
        time T1 equals X1 */
    status = euro_S(X1, X2, vol, r, q, t2 - t1, underlying_type, S_star);
    if (status == FALSE)
        return FALSE;
    *mu = r - q + vol * vol / 2.0;
    *var1 = vol * sqrt(t1);
    *var2 = vol * sqrt(t2);
    *a1 = (log(S / (*S_star)) + (*mu) * t1) / (*var1);
    *a2 = *a1 - *var1;
    *b1 = (log(S / X2) + (*mu) * t2) / (*var2);
    *b2 = *b1 - *var2;
    *sqrt_t1_t2 = sqrt(t1 / t2);
    return TRUE;
}

/*------------------------------------------------------------------------*/
/* Compound Option: Call on a Call                                        */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.417                                           */
/*------------------------------------------------------------------------*/
static short euro_call_call(double S, double X1, double X2, double vol,
                            double r, double q, double t1, double t2, double *result)
{
    double c, a1, a2, b1, b2, mu, var1, var2, sqrt_t1_t2, S_star;
    short status;

    if (vol < INPUT_UNDERFLOW)
    {
        *result = MAX(S * exp(-q * t2) - X2 * exp(-r * t2) - X1 * exp(-r * t1), 0.0);
        return TRUE;
    }

    if (t1 <= INPUT_UNDERFLOW || t2 <= t1)
    {
        *result = 0.0;
        return FALSE;
    }

    status = cmpd_setup(S, X1, X2, vol, r, q, t1, t2, &mu, &var1, &var2,
                        &a1, &a2, &b1, &b2, &sqrt_t1_t2, &S_star, OPT_CALL);
    if (status == FALSE)
    {
        *result = 0.0;
        return FALSE;
    }
    c = S * exp(-q * t2) * M(a1, b1, sqrt_t1_t2) - X2 * exp(-r * t2) * M(a2, b2, sqrt_t1_t2) - exp(-r * t1) * X1 * N(a2);
    *result = c;

    return TRUE;
}

/*------------------------------------------------------------------------*/
/* Compound Option: Call on a Put                                         */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.417                                           */
/*------------------------------------------------------------------------*/
static short euro_call_put(double S, double X1, double X2, double vol,
                           double r, double q, double t1, double t2, double *result)
{
    double c, a1, a2, b1, b2, mu, var1, var2, sqrt_t1_t2, S_star;
    short status;

    if (vol < INPUT_UNDERFLOW)
    {
        *result = MAX(X2 * exp(-r * t2) - S * exp(-q * t2) - X1 * exp(-r * t1), 0.0);
        return TRUE;
    }

    if (t1 <= INPUT_UNDERFLOW || t2 <= t1)
    {
        *result = 0.0;
        return FALSE;
    }

    status = cmpd_setup(S, X1, X2, vol, r, q, t1, t2, &mu, &var1, &var2,
                        &a1, &a2, &b1, &b2, &sqrt_t1_t2, &S_star, OPT_PUT);
    if (status == FALSE)
    {
        *result = 0.0;
        return FALSE;
    }
    c = X2 * exp(-r * t2) * M(-a2, -b2, sqrt_t1_t2) - S * exp(-q * t2) * M(-a1, -b1, sqrt_t1_t2) - exp(-r * t1) * X1 * N(-a2);
    *result = c;

    return TRUE;
}

/*------------------------------------------------------------------------*/
/* Compound Option: Put on a Call                                         */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.417                                           */
/*------------------------------------------------------------------------*/
static short euro_put_call(double S, double X1, double X2, double vol,
                           double r, double q, double t1, double t2, double *result)
{
    double c, a1, a2, b1, b2, mu, var1, var2, sqrt_t1_t2, S_star;
    short status;

    if (vol < INPUT_UNDERFLOW)
    {
        *result = MAX(X2 * exp(-r * t2) - S * exp(-q * t2) + X1 * exp(-r * t1), 0.0);
        return TRUE;
    }

    if (t1 <= INPUT_UNDERFLOW || t2 <= t1)
    {
        *result = 0.0;
        return FALSE;
    }

    status = cmpd_setup(S, X1, X2, vol, r, q, t1, t2, &mu, &var1, &var2,
                        &a1, &a2, &b1, &b2, &sqrt_t1_t2, &S_star, OPT_CALL);
    if (status == FALSE)
    {
        *result = 0.0;
        return FALSE;
    }
    c = X2 * exp(-r * t2) * M(-a2, b2, -sqrt_t1_t2) - S * exp(-q * t2) * M(-a1, b1, -sqrt_t1_t2) + exp(-r * t1) * X1 * N(-a2);
    *result = c;

    return TRUE;
}

/*------------------------------------------------------------------------*/
/* Compound Option: Put on a Put                                         */
/* From Hull & White "Options, Futures, and other Derivative Securities," */
/* 2nd Ed, Prentice Hall, p.417                                           */
/*------------------------------------------------------------------------*/
static short euro_put_put(double S, double X1, double X2, double vol,
                          double r, double q, double t1, double t2, double *result)
{
    double c, a1, a2, b1, b2, mu, var1, var2, sqrt_t1_t2, S_star;
    short status;

    if (vol < INPUT_UNDERFLOW)
    {
        *result = MAX(S * exp(-q * t2) - X2 * exp(-r * t2) + X1 * exp(-r * t1), 0.0);
        return TRUE;
    }

    if (t1 <= INPUT_UNDERFLOW || t2 <= t1)
    {
        *result = 0.0;
        return FALSE;
    }

    status = cmpd_setup(S, X1, X2, vol, r, q, t1, t2, &mu, &var1, &var2,
                        &a1, &a2, &b1, &b2, &sqrt_t1_t2, &S_star, OPT_PUT);
    if (status == FALSE)
    {
        *result = 0.0;
        return FALSE;
    }
    c = S * exp(-q * t2) * M(a1, -b1, -sqrt_t1_t2) - X2 * exp(-r * t2) * M(a2, -b2, -sqrt_t1_t2) + exp(-r * t1) * X1 * N(a2);
    *result = c;

    return TRUE;
}

/*------------------------------------------------------------------------*/
/* Compound option pricer                                                 */
/* Check for the above routines                                           */
/* Works by performing brute force integration over the underlying        */
/* european options                                                       */
/*------------------------------------------------------------------------*/
short cmpd_chk(double S, double X1, double X2, double vol,
               double r, double q, double t1, double t2,
               OPTTYPE opt_overlying, OPTTYPE opt_underlying,
               double *result)
{
    double epsilon, S_t1;
    int i;
    const int panels = 300;
    double mu, sqrt_t, sum, euro_val;
    mu = r - q;
    sqrt_t = sqrt(t1);

    sum = 0.0;
    for (i = 0; i < panels; i++)
    {
        epsilon = Ninv((i + 0.5) / (double)panels);
        S_t1 = S * exp((mu - vol * vol / 2.0) * t1 + vol * epsilon * sqrt_t);
        euro(S_t1, X2, vol, r, q, t2 - t1,
             opt_underlying, SENS_PRICE, &euro_val);
        if (opt_overlying == OPT_CALL)
            sum += MAX(euro_val - X1, 0);
        else
            sum += MAX(X1 - euro_val, 0);
    }
    sum /= panels;
    sum *= exp(-r * t1);
    *result = sum;
    return TRUE;
}

/*------------------------------------------------------------------------*/
/* Compound option pricer                                                 */
/* Wrapper to call the above routines                                     */
/*------------------------------------------------------------------------*/
short cmpd(double S, double X1, double X2, double vol,
           double r, double q, double t1, double t2,
           OPTTYPE opt_overlying, OPTTYPE opt_underlying,
           SENSTYPE sens, double *result)
{
    double dx, f1, f2, f3;
    short status;
    short (*val_fn)(double S, double X1, double X2, double vol,
                    double r, double q, double t1, double t2, double *result);

    /* Check Inputs */
    assert(result != NULL);
    if (S <= 0 || X1 <= 0 || X2 <= 0 || vol < 0 || r < 0 || q < 0)
    {
        *result = 0.0;
        return FALSE;
    }
    if ((S / X2) < 0.1 || (t2 - t1) <= INPUT_UNDERFLOW)
    {
        *result = 0.0;
        return FALSE;
    }

    switch (opt_overlying)
    {
    case OPT_CALL:
        switch (opt_underlying)
        {
        case OPT_CALL:
            val_fn = euro_call_call;
            break;
        case OPT_PUT:
            val_fn = euro_call_put;
            break;
        default: /* Invalid type */
            *result = 0.0;
            return FALSE;
        }
        break;
    case OPT_PUT:
        switch (opt_underlying)
        {
        case OPT_CALL:
            val_fn = euro_put_call;
            break;
        case OPT_PUT:
            val_fn = euro_put_put;
            break;
        default: /* Invalid type */
            *result = 0.0;
            return FALSE;
        }
        break;
    default: /* Invalid type */
        *result = 0.0;
        return FALSE;
    }

    switch (sens)
    {
    case SENS_PRICE:
        return (*val_fn)(S, X1, X2, vol, r, q, t1, t2, result);

#if 0
        case SENS_DELTA:  /* Patch here to call option checking routine */
            return cmpd_chk(S, X1, X2, vol, r, q, t1, t2, opt_overlying,
                    opt_underlying, result);
#endif

    case SENS_DELTA:
        dx = S / (double)eps;
        status = (*val_fn)(S, X1, X2, vol, r, q, t1, t2, &f1);
        status = status | (*val_fn)(S + dx, X1, X2, vol, r, q, t1, t2,
                                    &f2);
        *result = (f2 - f1) / dx;
        return status;
    case SENS_GAMMA:
        dx = S / (double)eps;
        status = (*val_fn)(S + dx, X1, X2, vol, r, q, t1, t2, &f1);
        status = status | (*val_fn)(S - dx, X1, X2, vol, r, q, t1, t2, &f2);
        status = status | (*val_fn)(S, X1, X2, vol, r, q, t1, t2, &f3);
        *result = (f1 + f2 - 2.0 * f3) / dx / dx;
        return status;
    case SENS_VEGA:
        dx = vol / (double)eps + 1.0 / 16.0 / (double)eps;
        status = (*val_fn)(S, X1, X2, vol, r, q, t1, t2, &f1);
        status = status | (*val_fn)(S, X1, X2, vol + dx, r, q, t1, t2, &f2);
        *result = (f2 - f1) / dx;
        return status;
    case SENS_THETA:
        dx = t2 / (double)eps + 1.0 / 16.0 / (double)eps;
        status = (*val_fn)(S, X1, X2, vol, r, q, t1, t2, &f1);
        status = status | (*val_fn)(S, X1, X2, vol, r, q, t1 - dx, t2 - dx, &f2);
        *result = (f2 - f1) / dx;
        return status;
    case SENS_RHO:
        if (r == q)
        { /* Option on futures */
            dx = q / (double)eps + 1.0 / 16.0 / (double)eps;
            status = (*val_fn)(S, X1, X2, vol, r, q, t1, t2, &f1);
            status = status | (*val_fn)(S, X1, X2, vol, r + dx, q + dx,
                                        t1, t2, &f2);
            *result = (f2 - f1) / dx;
            return status;
        }
        else
        {
            dx = r / (double)eps + 1.0 / 16.0 / (double)eps;
            status = (*val_fn)(S, X1, X2, vol, r, q, t1, t2, &f1);
            status = status | (*val_fn)(S, X1, X2, vol, r + dx, q,
                                        t1, t2, &f2);
            *result = (f2 - f1) / dx;
            return status;
        }
    default: /* Invalid Type */
        *result = 0.0;
        return FALSE;
    }
}
