#include "option.h"

/*------------------------------------------------------------------------*/
/* Convert a set of rates & dates to an interpolated rate curve with
/*		rates at specified intervals which are numbered starting from
/*		the first date, date[0]
/*
/* Inputs:
/*		rates = 1-D array of rates
/*		dates =  julian dates corresponding to each rate
/*		interval = # of months between each interpolated rate
/*
/* Outputs:
/*		TRUE = success, FALSE = fail
/*		interp_rates = array of rates at specified intervals
/*		interp_dates  = array of dates corresponding to the interpolated rates
/*------------------------------------------------------------------------*/
static short interpolate_rates(double *rates, date_t *dates, int num_rates,
							   int interval, date_t start_dt, double **interp_rates,
							   date_t **interp_dates,
							   int *num_dates)
{
	date_t last_dt, next_dt;
	ERRTYPE errcode;
	int i, dt_idx;
	date_t *new_dates;
	double *new_rates;

	assert(rates != NULL);
	assert(dates != NULL);
	assert(num_rates > 1);
	assert(interp_rates != NULL);
	assert(interp_dates != NULL);
	assert(num_dates != NULL);
	/* assert(interval >= 1); */
	assert(dates[num_rates - 1] > dates[0]);

	if (interval < 1)
		return FALSE;

	/* Count off new month intervals until we reach the last date in our
		input curve.  Place these dates (which seperated by the specified interval)
		in the array interp_dates
	*/

	/* Determine required # of dates */
	last_dt = start_dt;
	next_dt = last_dt;
	for (*num_dates = 0; last_dt <= dates[num_rates - 1]; (*num_dates)++)
	{
		errcode = addtime(interval, 0, 0, last_dt, &next_dt);
		if (!noerr(errcode))
			return FALSE;
		last_dt = next_dt;
	}

	/* Alloc space and store dates at specified new intervals */
	new_dates = (date_t *)calloc(*num_dates, sizeof(date_t));
	if (new_dates == NULL)
		return FALSE;
	new_dates[0] = start_dt;
	for (i = 1; i < *num_dates; i++)
		addtime(interval, 0, 0, new_dates[i - 1], &(new_dates[i]));

	/* Based on the dates in interp_dates interpolate from the original
		yield curve to find the yield corresponding to each date in
		interp_dates.
	*/
	new_rates = (double *)calloc(*num_dates, sizeof(double));
	if (new_rates == NULL)
		return FALSE;
	new_rates[0] = rates[0];
	dt_idx = 1;
	for (i = 1; i < *num_dates; i++)
	{
		/* Find next yield date > current interp date */
		while (new_dates[i] > dates[dt_idx] && dt_idx < num_rates - 1)
		{
			if (dates[dt_idx] <= dates[dt_idx - 1]) /* Check dates increasing */
				return FALSE;
			dt_idx++;
		}
		if (dt_idx >= num_rates) /* error */
			return FALSE;
		/* Interpolate rates */
		new_rates[i] = diffdate(new_dates[i], dates[dt_idx - 1]) / diffdate(dates[dt_idx], dates[dt_idx - 1]) * (rates[dt_idx] - rates[dt_idx - 1]) + rates[dt_idx - 1];
	}

	*interp_dates = new_dates;
	*interp_rates = new_rates;
	return TRUE;
}

/*------------------------------------------------------------------------*/
/*	Routine to bootstrap a treasury yield curve
/*
/*	Inputs:
/*		yields = 1-D array of bond equivalent treasury yields with
/*				the first yield = overnight rate
/*		dates =  julian dates corresponding to each yield
/*
/*	Outputs:
/*		TRUE = success, FALSE = fail
/*		ppdDisc = an array of discount factors at 6 month intervals out to the
/*					last yield date in dates;
/*    ppdtDisc_dates = array of dates corresponding to above discount factors
/*		num_disc = number of discount factors returned
/*
/* The algorithm for this routine is taken from two primary sources:
/* "By the Bootstraps, Risk Vol. 3, #6, June 1990, p.40-43"
/* "Options, Futures and Other Derivative Securities, Hull & White, 2nd ed,
/*  p. 84"
/*------------------------------------------------------------------------*/

static short BootstrapYieldsSixMo(double *yields, date_t *dates, int num_yields,
								  double **ppdDisc, date_t **ppdtDisc_dates,
								  int *num_disc)
{
	double *six_mo_yields;
	date_t *six_mo_dates;
	double *disc;
	double coupon, npv_coupons;
	int i, j, num_dates;
	short success;
	date_t today;

	/* By convention, first rate = overnight rate, first date = overnight date */
	addday(-1, dates[0], &today);

	/* Convert original set of yields to an interpolated curve with
		yields every six months
	*/
	success = interpolate_rates(yields, dates, num_yields, 6, today, &six_mo_yields,
								&six_mo_dates, &num_dates);
	if (success == FALSE)
		return FALSE;

	/* Override first yield = overnight yield */
	six_mo_dates[0] = dates[0];
	six_mo_yields[0] = yields[0];

	disc = (double *)calloc(num_dates, sizeof(double));
	if (disc == NULL)
		return FALSE;

	/* Overnight Discount */
	disc[0] = pow(1.0 + six_mo_yields[0] / 2.0, -2.0 / 365.25);

	/* Discount at 6 months */
	disc[1] = pow(1.0 + six_mo_yields[1] / 2.0,
				  -2.0 * (six_mo_dates[1] - today) / 365.25);

	/* For each date [i] find a new discount rate disc[i] s.t.
	$1 = NPV(coupons) + final_payment * disc[i]
		so
		since coupon = six_mo_yields[i]/2.0,
				final payment = $1 + coupon

	disc[i] = ($1 - NPV(coupons)) / final_payment
	where
				NPV(coupons) = SUM(j=1...i-1) (coupon * disc[j])
	*/
	for (i = 2; i < num_dates; i++)
	{
		npv_coupons = 0.0;
		coupon = six_mo_yields[i] / 2.0;
		for (j = 1; j < i; j++)
			npv_coupons += coupon * disc[j];
		disc[i] = (1.0 - npv_coupons) / (1.0 + coupon);
	}

	*ppdDisc = disc;
	*ppdtDisc_dates = six_mo_dates;
	*num_disc = num_dates;

	free(six_mo_yields);
	return TRUE;
}

/*------------------------------------------------------------------------*/
/* Convert a set of discounts to a set of discount rates with a specified
/* compounding frequency using an Actual/365.25 daycount convention
/*
/* Inputs:
/*		discs = 1-D array of discounts
/*		disc_dates =  julian dates corresponding to each discount
/*		num_disc = number of discount factors in array
/*		today = today date from which discounts are calculated
/*		freq = compounding frequency for use in calculating the rates
/*
/* Outputs:
/*		TRUE = success, FALSE = fail
/*		rates = discount rates under the specified compounding frequency
/*    rate_dates = dates corresponding to rates
/*------------------------------------------------------------------------*/
static short disc_to_rate(double *disc, date_t *disc_dates, int num_disc,
						  date_t today, int freq, double **ppdRates,
						  date_t **ppdRate_dates)
{
	double *rates, years;
	date_t *rate_dates;
	int i;

	rates = (double *)calloc(num_disc, sizeof(double));
	rate_dates = (date_t *)calloc(num_disc, sizeof(date_t));
	if (rates == NULL || rate_dates == NULL)
		return FALSE;

	for (i = 0; i < num_disc; i++)
	{
		if (disc[i] <= 0 || disc_dates[i] < today)
			return FALSE;
		if (disc_dates[i] != today)
		{
			years = diffdate(disc_dates[i], today) / 365.25;
			rates[i] = (pow(disc[i], -1.0 / (freq * years)) - 1.0) * freq;
		}
		else
		{
			rates[i] = 0.0;
		}
		rate_dates[i] = disc_dates[i];
	}

	*ppdRates = rates;
	*ppdRate_dates = rate_dates;
	return TRUE;
}

/*------------------------------------------------------------------------*/
/* Convert a set of discount rates with a specified compounding period
/* to a set of discount factors using an Actual/365.25 daycount convention
/*
/* Inputs:
/*		rates = discount rates under the specified compounding frequency
/*    rate_dates = dates corresponding to rates
/*		num_rate = number of rates in the array
/*		today = today date from which discounts are calculated
/*		freq = compounding frequency for use in calculating the rates
/* Outputs:
/*		TRUE = success, FALSE = fail
/*		discs = 1-D array of discounts
/*		disc_dates =  julian dates corresponding to each discount
/*
/*------------------------------------------------------------------------*/
static short rate_to_disc(double *rates, date_t *rate_dates, int num_rates,
						  date_t today, int freq, double **ppdDisc,
						  date_t **ppdDisc_dates)
{
	double *disc, years;
	date_t *disc_dates;
	int i;

	disc = (double *)calloc(num_rates, sizeof(double));
	disc_dates = (date_t *)calloc(num_rates, sizeof(date_t));
	if (disc == NULL || disc_dates == NULL)
		return FALSE;

	for (i = 0; i < num_rates; i++)
	{
		if (rates[i] < 0 || rate_dates[i] < today)
			return FALSE;
		if (rate_dates[i] > today)
		{
			years = diffdate(rate_dates[i], today) / 365.25;
			disc[i] = pow(1.0 + rates[i] / freq, -freq * years);
		}
		else
		{
			disc[i] = 1.0;
		}
		disc_dates[i] = rate_dates[i];
	}

	*ppdDisc = disc;
	*ppdDisc_dates = disc_dates;
	return TRUE;
}

/*------------------------------------------------------------------------*/
/*	Routine to bootstrap a treasury yield curve, and return discount
/* rates at the specified monthly intervals
/*
/*	Inputs:
/*		yields = 1-D array of bond equivalent treasury yields with
/*				the first yield = overnight rate & first date = today
/*		dates =  julian dates corresponding to each yield
/*		interval = # of months between each discount factor to be returnd
/*
/*	Outputs:
/*		TRUE = success, FALSE = fail
/*		ppdDisc = an array of discount factors at the specified intervals out
/*					to the last yield date in dates;
/*    ppdtDisc_dates = array of dates corresponding to above discount factors
/*		num_disc = number of discount factors returned
/*
/* The algorithm for this routine is taken from two primary sources:
/* "By the Bootstraps, Risk Vol. 3, #6, June 1990, p.40-43"
/* "Options, Futures and Other Derivative Securities, Hull & White, 2nd ed,
/*  p. 84"
/*------------------------------------------------------------------------*/
short BootstrapYields(double *yields, date_t *dates, int num_yields,
					  int interval, double **ppdDisc, date_t **ppdtDisc_dates,
					  int *num_disc)
{
	short success;
	double *disc_rates, *interp_disc_rates, *interp_disc;
	date_t *disc_rate_dates, *interp_disc_rate_dates, *interp_disc_dates;
	date_t today;
	int num_interp_disc_rates;

	/* Find discount factors at 6 month intervals */
	success = BootstrapYieldsSixMo(yields, dates, num_yields, ppdDisc,
								   ppdtDisc_dates, num_disc);
	if (success == FALSE)
		return FALSE;

	/* By convention, first rate = overnight rate, first date = overnight date */
	/* Set today = date[0] -1 */
	addday(-1, dates[0], &today);

	/* Convert discount factors to semiannual discount rates */
	success = disc_to_rate(*ppdDisc, *ppdtDisc_dates, *num_disc,
						   today, 2, &disc_rates, &disc_rate_dates);
	if (success == FALSE)
		return FALSE;

	free(*ppdDisc);
	free(*ppdtDisc_dates);

	/* Interpolate semmiannual rates to rates at the desired interval */
	success = interpolate_rates(disc_rates, disc_rate_dates, *num_disc,
								interval, today, &interp_disc_rates,
								&interp_disc_rate_dates,
								&num_interp_disc_rates);
	if (success == FALSE)
		return FALSE;

	/* Override first yield = overnight yield */
	interp_disc_rate_dates[0] = dates[0];
	interp_disc_rates[0] = yields[0];

	/* Convert interpolated semmiannual rates back into discount factors */
	success = rate_to_disc(interp_disc_rates, interp_disc_rate_dates,
						   num_interp_disc_rates,
						   today, 2, &interp_disc,
						   &interp_disc_dates);
	if (success == FALSE)
		return FALSE;

	/* Free uneeded intermediate results */
	free(disc_rates);
	free(disc_rate_dates);

	free(interp_disc_rates);
	free(interp_disc_rate_dates);

	/* Return the final set of discount factors */
	*ppdDisc = interp_disc;
	*ppdtDisc_dates = interp_disc_dates;
	*num_disc = num_interp_disc_rates;

	return TRUE;
}
