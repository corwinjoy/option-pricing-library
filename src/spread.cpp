#include "option.h"

/*------------------------------------------------------------------------*/
/* Spread Call: Max( (S1 - S2) - X, 0)                                    */
/* The method used here is due to Ravindran, "Low Fat Spreads," Risk,     */
/* Vol. 6, #10, Oct. 1993                                                 */
/*------------------------------------------------------------------------*/
static short spread_call(double S1, double S2, double X, double vol1,
                         double vol2, double rho,
                         double r, double q1, double q2, double t,
                         double *result)
{
    double mu, vol_sqr, half_vol_sqr, vol, mu1, mu2, g1, g2,
        prob, step, Kstar, var2, value, rv1v2, ba, ca, d1, d2,
        epsilon, log_S2T;
    double a = 1.0,
           b = -1.0,
           c = -X;
    const int panels = 100;

    if (t <= INPUT_UNDERFLOW)
    {
        *result = MAX(S1 - S2 - X, 0.0);
        return TRUE;
    }

    g1 = r - q1;
    g2 = r - q2;

    mu1 = log(S1) + (g1 - 0.5 * vol1 * vol1) * t;
    mu2 = log(S2) + (g2 - 0.5 * vol2 * vol2) * t;
    var2 = vol2 * sqrt(t);
    vol_sqr = vol1 * vol1 * (1.0 - rho * rho) * t;
    vol = sqrt(vol_sqr);
    half_vol_sqr = vol_sqr / 2.0;
    rv1v2 = rho * vol1 / vol2;
    ba = b / a;
    ca = c / a;

    value = 0.0;
    step = 1.0 / (double)panels;

    if (vol < INPUT_UNDERFLOW)
    {
        for (prob = 0.5 / (double)panels; prob < 1.0; prob += step)
        {
            epsilon = Ninv(prob);

            /* Valuation for a given level of epsilon */
            log_S2T = epsilon * var2;
            mu = mu1 + rv1v2 * log_S2T;
            log_S2T += mu2;
            Kstar = -ba * exp(log_S2T) - ca;

            value += MAX(exp(mu + half_vol_sqr) - Kstar, 0.0);
        }
    }
    else
    {
        for (prob = 0.5 / (double)panels; prob < 1.0; prob += step)
        {
            epsilon = Ninv(prob);

            /* Valuation for a given level of epsilon */
            log_S2T = epsilon * var2;
            mu = mu1 + rv1v2 * log_S2T;
            log_S2T += mu2;
            Kstar = -ba * exp(log_S2T) - ca;

            if (Kstar <= 0.0)
            { /* Handle negative Kstar for puts */
                value += exp(mu + half_vol_sqr) - Kstar;
            }
            else
            {
                d1 = (mu - log(Kstar)) / vol;
                d2 = d1 + vol;

                value += exp(mu + half_vol_sqr) * N(d2) - Kstar * N(d1);
            }
        }
    }

    value *= a;
    value /= (double)panels;

    *result = value * exp(-r * t);
    return TRUE;
}

/*------------------------------------------------------------------------*/
/* Spread Put: Max( X - (S1 - S2), 0)                                     */
/* The method used here is due to Ravindran, "Low Fat Spreads," Risk,     */
/* Vol. 6, #10, Oct. 1993                                                 */
/*------------------------------------------------------------------------*/
static short spread_put(double S1, double S2, double X, double vol1,
                        double vol2, double rho,
                        double r, double q1, double q2, double t,
                        double *result)
{

    /* Use fact that Max(X - (S1 - S2), 0) = Max( (S2 - S1) - (-X), 0 )
        to value spread put as a modified spread call
    */
    return spread_call(S2, S1, -X, vol2, vol1, rho, r, q2, q1, t, result);
}

/*------------------------------------------------------------------------*/
/* Spread Options                                                         */
/* Wrapper to call the above routines                                     */
/*------------------------------------------------------------------------*/
short spread(double S1, double S2, double X, double vol1, double vol2,
             double rho, double r, double q1, double q2, double t,
             OPTTYPE call_put, SENSTYPE sens,
             double *result)
{

    double dx, f1, f2, f3;
    short status;
    short (*val_fn)(double S1, double S2, double X, double vol1,
                    double vol2, double rho,
                    double r, double q1, double q2, double t,
                    double *result);

    /* Check Input Parameters */
    if (call_put != OPT_CALL && call_put != OPT_PUT)
    {
        *result = 0.0;
        return FALSE;
    }

    /* Check Inputs */
    assert(result != NULL);
    if (S1 <= 0 || S2 <= 0 || vol1 < INPUT_UNDERFLOW || vol2 < INPUT_UNDERFLOW || r < 0 || q1 < 0 || q2 < 0 || rho < -1 || rho > 1)
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
        val_fn = spread_call;
    else
        val_fn = spread_put;

    switch (sens)
    {
    case SENS_PRICE:
        return (*val_fn)(S1, S2, X, vol1, vol2, rho, r, q1, q2, t, result);
    case SENS_DELTA:
        dx = S1 / (double)eps;
        status = (*val_fn)(S1, S2, X, vol1, vol2, rho, r, q1, q2, t, &f1);
        status = status | (*val_fn)(S1 + dx, S2, X, vol1, vol2, rho, r, q1,
                                    q2, t, &f2);
        *result = (f2 - f1) / dx;
        return status;
    case SENS_GAMMA:
        dx = S1 / (double)eps;
        status = (*val_fn)(S1 + dx, S2, X, vol1, vol2, rho, r, q1, q2, t, &f1);
        status = status | (*val_fn)(S1 - dx, S2, X, vol1, vol2, rho, r, q1,
                                    q2, t, &f2);
        status = status | (*val_fn)(S1, S2, X, vol1, vol2, rho, r, q1,
                                    q2, t, &f3);
        *result = (f1 + f2 - 2.0 * f3) / dx / dx;
        return status;
    case SENS_VEGA:
        dx = vol1 / (double)eps * 4.0 + 1.0 / 16.0 / (double)eps;
        status = (*val_fn)(S1, S2, X, vol1, vol2, rho, r, q1, q2, t, &f1);
        status = status | (*val_fn)(S1, S2, X, vol1 + dx, vol2, rho, r, q1,
                                    q2, t, &f2);
        *result = (f2 - f1) / dx;
        return status;
    case SENS_DELTA2:
        dx = S2 / (double)eps;
        status = (*val_fn)(S1, S2, X, vol1, vol2, rho, r, q1, q2, t, &f1);
        status = status | (*val_fn)(S1, S2 + dx, X, vol1, vol2, rho, r, q1,
                                    q2, t, &f2);
        *result = (f2 - f1) / dx;
        return status;
    case SENS_GAMMA2:
        dx = S2 / (double)eps;
        status = (*val_fn)(S1, S2 + dx, X, vol1, vol2, rho, r, q1, q2, t, &f1);
        status = status | (*val_fn)(S1, S2 - dx, X, vol1, vol2, rho, r, q1,
                                    q2, t, &f2);
        status = status | (*val_fn)(S1, S2, X, vol1, vol2, rho, r, q1,
                                    q2, t, &f3);
        *result = (f1 + f2 - 2.0 * f3) / dx / dx;
        return status;
    case SENS_VEGA2:
        dx = vol2 / (double)eps * 4.0 + 1.0 / 16.0 / (double)eps;
        status = (*val_fn)(S1, S2, X, vol1, vol2, rho, r, q1, q2, t, &f1);
        status = status | (*val_fn)(S1, S2, X, vol1, vol2 + dx, rho, r, q1, q2, t, &f2);
        *result = (f2 - f1) / dx;
        return status;
    case SENS_THETA:
        dx = t / (double)eps + 1.0 / 16.0 / (double)eps;
        status = (*val_fn)(S1, S2, X, vol1, vol2, rho, r, q1, q2, t, &f1);
        status = status | (*val_fn)(S1, S2, X, vol1, vol2, rho, r, q1, q2, t - dx, &f2);
        *result = (f2 - f1) / dx;
        return status;
    case SENS_RHO:
        dx = r / (double)eps + 1.0 / 16.0 / (double)eps;
        status = (*val_fn)(S1, S2, X, vol1, vol2, rho, r, q1, q2, t, &f1);
        if (r == q1 && r == q2)
            status = status | (*val_fn)(S1, S2, X, vol1, vol2, rho, r + dx,
                                        q1 + dx, q2 + dx, t, &f2);
        else
            status = status | (*val_fn)(S1, S2, X, vol1, vol2, rho, r + dx,
                                        q1, q2, t, &f2);
        *result = (f2 - f1) / dx;
        return status;
    case SENS_CORR:
        dx = rho / (double)eps * 4.0 + 1.0 / 16.0 / (double)eps;
        if (rho > 0.90)
            dx = -dx;
        status = (*val_fn)(S1, S2, X, vol1, vol2, rho, r, q1, q2, t, &f1);
        status = status | (*val_fn)(S1, S2, X, vol1, vol2, rho + dx, r, q1, q2, t, &f2);
        *result = (f2 - f1) / dx;
        return status;
    default: /* Invalid Sensitivity */
        *result = 0.0;
        return FALSE;
    }
}
