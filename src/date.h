#ifndef _DATE_H_
#define _DATE_H_
/* Header file for dates routines */
#define FALSE 0
#define TRUE 1

typedef long date_t;
enum ERRTYPE
{
    ERR_FAIL,
    ERR_OK
};

short noerr(ERRTYPE errcode);
ERRTYPE julday(int mm, int id, int iyyy, date_t *juldate);
ERRTYPE unjulday(date_t julian, int *mm, int *id, int *iyyy);
ERRTYPE addtime(int mm, int id, int iyyy, date_t olddate, date_t *newdate);
ERRTYPE addday(int id, date_t olddate, date_t *newdate);
char *jultostr(date_t juldate);
double diffdate(date_t date1, date_t date2);
double jultod(date_t date);

#endif