#include "date.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

/*------------------------------------------------------------------------*/
/* Check return status for errors, returns true if no error
/*------------------------------------------------------------------------*/
short noerr(ERRTYPE errcode)
{
    if (errcode == ERR_OK)
        return TRUE;
    else
        return FALSE;
}
#define IGREG (15 + 31L * (10 + 12L * 1582))
/*------------------------------------------------------------------------*/
/* This routine computes the Julian Day Number that begins at noon of the
/* calendar date specified by month mm, day id, and year iyyy.
/* The calculated date is returned in juldate.  The routine returns
/* TRUE if the calculation was successful and FALSE otherwisePositive
/* year sifnifies AD; negative, BC. Remember that the year after 1 B.C.
/* was 1 A.D.
/* Taken from "Numerical Recipes in C," 2nd Ed,
/*	Press, Teukolsky et al., p.11
/*------------------------------------------------------------------------*/

ERRTYPE julday(int mm, int id, int iyyy, date_t *juldate)
{
    date_t jul;
    int ja, jy = iyyy, jm;

    if (jy == 0)
        return ERR_FAIL;
    if (jy < 0)
        ++jy;
    if (mm > 2)
    {
        jm = mm + 1;
    }
    else
    {
        --jy;
        jm = mm + 13;
    }
    jul = (long)(floor(365.25 * jy) + floor(30.6001 * jm) + id + 1720995L);
    if (id + 31L * (mm + 12L * iyyy) >= IGREG)
    {
        ja = (int)(0.01 * jy);
        jul += 2 - ja + (int)(0.25 * ja);
    }
    *juldate = jul;

    return ERR_OK;
}
#undef IGREG

#define IGREG 2299161L
/*------------------------------------------------------------------------*/
/* Inverse of the function julday given above.  Here julian is input
/* as a Julian Day Number, and the routine outputs mm, id, and iyyy
/* as the month day and year on which the specified Julian Day
/* started at noon.
/* Taken from "Numerical Recipes in C," 2nd Ed,
/*	Press, Teukolsky et al., p.14
/*------------------------------------------------------------------------*/
ERRTYPE unjulday(date_t julian, int *mm, int *id, int *iyyy)
{
    date_t ja, jalpha, jb, jc, jd, je;

    if (julian >= IGREG)
    {
        jalpha = (long)(((float)(julian - 1867216L) - 0.25) / 36524.25);
        ja = julian + 1 + jalpha - (long)(0.25 * jalpha);
    }
    else
        ja = julian;
    jb = ja + 1524;
    jc = (long)(6680.0 + ((float)(jb - 2439870L) - 122.1) / 365.25);
    jd = (long)(365 * jc + (0.25 * jc));
    je = (long)((jb - jd) / 30.6001);
    *id = (int)(jb - jd - (long)(30.6001 * je));
    *mm = (int)(je - 1);
    if (*mm > 12)
        *mm -= 12;
    *iyyy = (int)(jc - 4715);
    if (*mm > 2)
        --(*iyyy);
    if (*iyyy <= 0)
        --(*iyyy);

    return ERR_OK;
}
#undef IGREG

/*------------------------------------------------------------------------*/
/* This routine takes the date stored in olddate and adds to it the
/* specified number of months mm, days id, and years iyyy to obtain
/* a new date which is returned in newdate
/*------------------------------------------------------------------------*/
ERRTYPE addtime(int mm, int id, int iyyy, date_t olddate, date_t *newdate)
{
    int old_mm, old_id, old_iyyy;
    int new_mm, new_id, new_iyyy;
    div_t month_div;
    ERRTYPE errcode;

    errcode = unjulday(olddate, &old_mm, &old_id, &old_iyyy);
    if (!noerr(errcode))
        return errcode;

    new_id = old_id + id; /* routine julday can handle id > 31 */

    /* add the new months to old_mm, if this adds up to more than 12
        add the quotient in years to old_yy and set new_mm equal
        to the remainder
    */
    new_mm = old_mm + mm;
    month_div = div(new_mm, 12);
    old_iyyy += month_div.quot;
    new_mm = month_div.rem;

    new_iyyy = old_iyyy + iyyy;

    return julday(new_mm, new_id, new_iyyy, newdate);
}

/*------------------------------------------------------------------------*/
/* This routine takes the date stored in olddate and adds to it the
/* specified number of days id, to obtain
/* a new date which is returned in newdate
/*------------------------------------------------------------------------*/
ERRTYPE addday(int id, date_t olddate, date_t *newdate)
{

    *newdate = olddate + id;
    return ERR_OK;
}

/*------------------------------------------------------------------------*/
/* Convert a Julian date to a string in the form MON-DD-YYYY
/*------------------------------------------------------------------------*/
char *jultostr(date_t juldate)
{
    static char szDate[12];
    static char szErr[] = "#DATE_ERR";
    static char mm_name[13][4] = {"", "JAN", "FEB", "MAR", "APR", "MAY", "JUN",
                                  "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};
    int mm, id, iyyy;
    ERRTYPE errcode;

    errcode = unjulday(juldate, &mm, &id, &iyyy);
    if (!noerr(errcode))
        return szErr;

    sprintf(szDate, "%3s-%2d-%4d", mm_name[mm], id, iyyy);
    return szDate;
}

/*------------------------------------------------------------------------*/
/* Returns date1 - date2 expressed in days
/* Should always use this routine to find the difference between two
/* Julian days - if we later change to a type that contains seconds
/* direct subtraction may no longer work
/*------------------------------------------------------------------------*/
double diffdate(date_t date1, date_t date2)
{
    date_t diff;
    double result;
    diff = date1 - date2;
    result = (double)diff;
    return result;
}

/*------------------------------------------------------------------------*/
/* Converts date to a double
/*------------------------------------------------------------------------*/
double jultod(date_t date)
{
    return (double)date;
}
