#include "option.h"
#include <limits.h>

/* Module level storage variables */
static double *m_epsilon = NULL;
static double *m_corr_epsilon = NULL;
static double *m_mean = NULL;
static double *m_var = NULL;
static double **m_chol_corr = NULL;

/*------------------------------------------------------------------------*/
/* Generate Correlated Settlement Prices S_t[i] based on the starting     */
/* prices S[i],                                                           */
/*	the Cholesky Decomposition of the correlation matrix chol_corr[i][j],  */
/* the lognormal means mean[i] = ((r - q[i]) - vol[i] ^ 2 / 2) * t[i]     */
/* and the lognormal variances var[i] = vol[i] * sqrt(t[i])               */
/* nCom = # of commodities in the basket                                  */
/* nSim = # of Quasi Monte Carlo simulations that will be performed       */
/* For a description of Monte Carlo see Hull, "Options, Futures and Other */
/* Derivative Securities, 2nd Ed., p. 331"                                */
/*------------------------------------------------------------------------*/
static short GeneratePrices(double *S,
							int nCom, unsigned long nSim, double *S_t)
{
	int i, j;
	short status;

	assert(m_epsilon != NULL);		/* Must have alloced storage space for rnd #s */
	assert(m_corr_epsilon != NULL); /* Must have alloced storage space for rnd #s */

	/* Generate Quasi-Random numbers and store them in m_epsilon */
	status = rvsf(m_epsilon, nCom, nSim, SEQ_EXACT);
	if (status == FALSE)
		return FALSE;

	/* Convert the iid normal deviates in into correlated random #s*/
	for (i = 0; i < nCom; i++)
	{
		m_corr_epsilon[i] = 0.0;
		for (j = 0; j <= i; j++)
			m_corr_epsilon[i] += m_epsilon[j] * m_chol_corr[i][j];
	}

	/* Convert the correlated random numbers into prices */
	for (i = 0; i < nCom; i++)
		S_t[i] = S[i] * exp(m_mean[i] + m_corr_epsilon[i] * m_var[i]);

	return TRUE;
}

/*------------------------------------------------------------------------*/
/* Routine to alloc memory and initialize module level variables for      */
/* basket option pricing                                                  */
/*------------------------------------------------------------------------*/
short GenPricesSetup(double *vol, double **corr,
					 double r, double *q, double *t, int nCom)
{
	int i;

	m_epsilon = (double *)calloc(nCom, sizeof(double));
	m_corr_epsilon = (double *)calloc(nCom, sizeof(double));
	if (m_epsilon == NULL || m_corr_epsilon == NULL)
		return FALSE;

	/* initialize the lognormal mean and variance for each commodity via */
	/* lognormal means mean[i] = ((r - q[i]) - vol[i] ^ 2 / 2) * t[i]    */
	/* lognormal variances var[i] = vol[i] * sqrt(t[i])                  */
	m_mean = (double *)calloc(nCom, sizeof(double));
	m_var = (double *)calloc(nCom, sizeof(double));
	if (m_mean == NULL || m_var == NULL)
		return FALSE;

	for (i = 0; i < nCom; i++)
	{
		if (t[i] < 0)
		{
			m_mean[i] = 0.0;
			m_var[i] = 0.0;
		}
		else
		{
			m_mean[i] = ((r - q[i]) - vol[i] * vol[i] / 2.0) * t[i];
			m_var[i] = vol[i] * sqrt(t[i]);
		}
	}

	m_chol_corr = Cholesky(corr, nCom);

	return TRUE;
}

/*------------------------------------------------------------------------*/
/* Routine to free memory allocated by basket_setup                       */
/*------------------------------------------------------------------------*/
void GenPricesCleanup(void)
{
	free(m_epsilon);
	m_epsilon = NULL;

	free(m_corr_epsilon);
	m_corr_epsilon = NULL;

	free(m_mean);
	m_mean = NULL;

	free(m_var);
	m_var = NULL;

	deleteMatrix(m_chol_corr);
	m_chol_corr = NULL;
}

/*------------------------------------------------------------------------*/
/* Basket Call Option: Max( w1 * S1 + w2 * S2 + ... + wn * Sn - X, 0)     */
/*------------------------------------------------------------------------*/
short basket_price(double *S, double X, double *w,
				   double r, double exp_t, int nCom,
				   unsigned long nSim, OPTTYPE call_put, double *result)
{
	unsigned long i;
	int j;
	double basket_index, sum;
	short status;
	double *S_t;

	/* Basket Setup should be called before this routine */
	assert(m_epsilon != NULL);
	assert(m_corr_epsilon != NULL);
	assert(m_mean != NULL);
	assert(m_var != NULL);

	S_t = (double *)calloc(nCom, sizeof(double));
	if (S_t == NULL)
		return FALSE;

	rvsf_reset(); /* Reset Faure Sequence Generator To Start At Beginning */
	sum = 0.0;
	for (i = 0; i < nSim; i++)
	{
		status = GeneratePrices(S, nCom, nSim, S_t);
		if (status == FALSE)
			return FALSE;
		basket_index = 0.0;
		for (j = 0; j < nCom; j++)
			basket_index += w[j] * S_t[j];
		if (call_put == OPT_CALL)
			sum += MAX(basket_index - X, 0);
		else
			sum += MAX(X - basket_index, 0);
	}

	sum /= nSim;
	sum *= exp(-r * exp_t);

	*result = sum;

	free(S_t);
	return TRUE;
}

/*------------------------------------------------------------------------*/
/* Basket Call Option: Max( w1 * S1 + w2 * S2 + ... + wn * Sn - X, 0)     */
/* or                                                                     */
/* Basket Put Option: Max( X - (w1 * S1 + w2 * S2 + ... + wn * Sn), 0)    */
/* Where S[i] = Price of the ith commodity                                */
/* X = strike price for the basket                                        */
/* vol[i] = volatilities of the commodities in the basket                 */
/* w[i] = weight in the basket for the ith commodity                      */
/* corr[i][j] = correlation matrix for the commodities in the basket      */
/* r = risk free discount rate                                            */
/* q[i] = dividend yield for the ith commodity                            */
/* t[i] = time at which price S[i] settles                                */
/* exp_t = time at which option expires                                   */
/* nCom = # of commodities in the basket                                  */
/* nSim = # of Quasi Monte Carlo simulations to perform                   */
/* Wrapper to call the above routines                                     */
/*------------------------------------------------------------------------*/
short basket(double *S, double X, double *vol, double *w, double **corr, double r,
			 double *q, double *t, double exp_t, int nCom,
			 unsigned long nSim, OPTTYPE call_put, SENSTYPE sens, double *result)
{
#if 0
	double dx, f1, f2, f3;
#endif
	short status;
	short (*val_fn)(double *S, double X, double *w,
					double r, double exp_t, int nCom,
					unsigned long nSim, OPTTYPE call_put, double *result);
	int i;

	/* Check for valid input pointers */
	if (result == NULL)
		return FALSE;

	if (S == NULL || vol == NULL || w == NULL || corr == NULL || q == NULL || t == NULL)
	{
		*result = 0.0;
		return FALSE;
	}

	if (*corr == NULL)
	{
		*result = 0.0;
		return FALSE;
	}

	/* Check Input Parameters */
	if (call_put != OPT_CALL && call_put != OPT_PUT)
	{
		*result = 0.0;
		return FALSE;
	}
	/* Check Inputs */
	if (r < 0)
	{
		*result = 0.0;
		return FALSE;
	}

	/* Check Inputs */
	for (i = 0; i < nCom; i++)
		if (S[i] <= 0 || vol[i] < INPUT_UNDERFLOW || q[i] < 0)
		{
			*result = 0.0;
			return FALSE;
		}

	val_fn = basket_price;
	status = GenPricesSetup(vol, corr, r, q, t, nCom);
	if (status == FALSE)
		return FALSE;

	switch (sens)
	{
	case SENS_PRICE:
		status = (*val_fn)(S, X, w, r, exp_t, nCom,
						   nSim, call_put, result);
		GenPricesCleanup();
		return status;
#if 0
		case SENS_DELTA:
			dx = S / (double) eps;
			status = (*val_fn)(S, X, Q, vol, r, q, t, &f1);
			status = status | (*val_fn)(S+dx, X, Q, vol, r, q, t, &f2);
			*result = (f2 - f1) / dx;
			return status;
		case SENS_GAMMA:
			dx = S / (double) eps;
			status = (*val_fn)(S+dx, X, Q, vol, r, q, t, &f1);
			status = status | (*val_fn)(S-dx, X, Q, vol, r, q, t, &f2);
			status = status | (*val_fn)(S, X, Q, vol, r, q, t, &f3);
			*result = (f1 + f2 - 2.0 * f3) / dx / dx;
			return status;
		case SENS_VEGA:
			dx = vol / (double) eps + 1.0 / 16.0 / (double) eps;
			status = (*val_fn)(S, X, Q, vol, r, q, t, &f1);
			status = status | (*val_fn)(S, X, Q, vol+dx, r, q, t, &f2);
			*result = (f2 - f1) / dx;
			return status;
		case SENS_THETA:
			dx = t / (double) eps;
			status = (*val_fn)(S, X, Q, vol, r, q, t, &f1);
			status = status | (*val_fn)(S, X, Q, vol, r, q, t-dx, &f2);
			*result = (f2 - f1) / dx;
			return status;
		case SENS_RHO:
			dx = r / (double) eps + 1.0 / 16.0 / (double) eps;
			status = (*val_fn)(S, X, Q, vol, r, q, t, &f1);
			if ( r == q)
				status = status | (*val_fn)(S, X, Q, vol, r+dx, q+dx, t, &f2);
			else
				status = status | (*val_fn)(S, X, Q, vol, r+dx, q, t, &f2);
			*result = (f2 - f1) / dx;
			return status;
#endif
	default: /* Invalid sensitivity */
		*result = 0.0;
		return FALSE;
	}
}
