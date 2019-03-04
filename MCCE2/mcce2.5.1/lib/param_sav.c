/*******************************************************************************
 * NAME
 *        param_sav - save a value indexed by 3 key strings
 *
 * SYNOPSIS
 *        #include <gdbm.h>
 *        #include <mcce.h>
 *
 *        int param_sav(char *key1, char *key2, char *key3, void *value, int s);
 *
 * DESCRIPTION
 *        The param_sav() function saves the value indexed by 3 key strings.  The
 *        leading and ending spaces of these keys are stripped off in saving,  so
 *        it is the same to inquire a value by either " CB " or "CB".   The value
 *        can be a pointer to any data type or array, as long as the size of  the
 *        data entry, s, is passed in correctly. It can be retreived later by the
 *        function param_get().  Use "-lgdbm" option to compile these programs,
 *        just like use "-lm" to compile programs calling math functions.
 *
 * RETURN VALUE
 *        The returned integer is 0 upon success or -1 upon failure.
 *
 * SEE ALSO
 *        db_open, db_close, iatom, param_get
 *
 * EXAMPLE
 *       #include <stdio.h>
 *       #include <gdbm.h>
 *       #include "mcce.h"
 *
 *       int main()
 *       {  char  conformers[] = "ASP0, ASP1, ASP-";
 *          int   i = 10;
 *          float pKa = 5.7;
 *
 *          char  ret_conf[160];
 *          int   ret_i = 0;
 *          float ret_pKa = 0.0;
 *
 *          db_open();
 *
 *          param_sav("CONFLIST", "ASP",  "",   &conformers, sizeof(conformers));
 *          param_sav("IATOM   ", "ASP0", " CG ", &i, sizeof(int));
 *          param_sav("PKA     ", "ASP-", "",   &pKa, sizeof(float));
 *
 *          param_get("CONFLIST", "ASP", "", &ret_conf);
 *          printf("%s\n",ret_conf);
 *
 *          ret_i = iatom("ASP0", "CG");
 *          printf("%d\n",ret_i);
 *
 *          param_get("PKA", "ASP-", "", &ret_pKa);
 *          printf("%.2f\n",ret_pKa);
 *
 *          db_close();
 *          return 0;
 *       }
 *
 * DEPENDENCY
 *        strip_spc
 *
 * AUTHOR
 *        Junjun Mao, 06/03/2003
 *******************************************************************************/

#include <string.h>
#include <gdbm.h>
#include "mcce.h"

extern GDBM_FILE param_db;    /* defined by db_open() */

int param_sav(char *key1, char *key2, char *key3, void *value, int s)
{  datum pkey, pvalue;
   char key[MAXCHAR_LINE];
   char sbuff[MAXCHAR_LINE];

   /* convert 3 key strings to one key, leading and ending spaces stripped */
   strip(key, key1);
   strip(sbuff, key2); strcat(key, sbuff);
   strip(sbuff, key3); strcat(key, sbuff);

   /* condtruct key, value pair */
   pkey.dptr = key;
   pkey.dsize = strlen(key);
   pvalue.dptr = value;
   pvalue.dsize = s;    /* The data size is passed in rather measured by strlen().
                         * When a generic pointer was passed in, the variable type
                         * was lost.     It became not reliable to detect the data
                         * length by the terminating NULL character,  thought most
                         * time strlen(value) would return the right size. */

   return gdbm_store(param_db, pkey, pvalue, GDBM_REPLACE);
}
