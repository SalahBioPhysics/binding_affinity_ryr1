/*******************************************************************************
 * NAME
 *        param_get - retrieve the parameter value as a string by 3 key strings
 *
 * SYNOPSIS
 *        #include <gdbm.h>
 *        #include <mcce.h>
 *
 *        int param_get(char *key1, char *key2, char *key3, void *value);
 *
 * DESCRIPTION
 *        The  param_get()  function  retrieves  the value in a pointer  by 3 key
 *        strings. The 3 key strings are defined in parameter files and loaded to
 *        the database by param_sav(),  but  the keys  are not  sensitive to  the
 *        leading and ending spaces. For example,  atom name " CB "  and "CB" are
 *        the same. The value is of generic type, can be cast into the appropiate
 *        data type and array.  Refer to the example for details.   Use  "-lgdbm"
 *        option to compile these  programs,   like use "-lm" to compile programs
 *        calling math functions.
 *
 * RETURN VALUE
 *        The returned integer is 0 upon success or -1 upon failure.
 *
 * SEE ALSO
 *        db_open, db_close, iatom, param_sav
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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gdbm.h>
#include "mcce.h"

extern GDBM_FILE param_db;

int param_get(char *key1, char *key2, char *key3, void *value)
/* WARNING: if the expected value is a string, it has to be long enough
 *          (>=MAXCHAR_LINE) to accept the stored value to avoid over
 *          boundary writing.
 */
{  datum pkey, pvalue;
   char key[MAXCHAR_LINE];
   char sbuff[MAXCHAR_LINE];

   /* convert 3 key strings to one key, leading and ending spaces stripped */
   strip(key, key1);
   strip(sbuff, key2); strcat(key, sbuff);
   strip(sbuff, key3); strcat(key, sbuff);

   /* make the key string to be the database search key */
   pkey.dptr = key;
   pkey.dsize = strlen(key);

   /* get the value */
   if (!gdbm_exists(param_db, pkey)) return -1; /* failure */
   else {                                       /* success */
      pvalue = gdbm_fetch(param_db, pkey);
      memcpy(value, pvalue.dptr, pvalue.dsize);
      free(pvalue.dptr);
      return 0;
   }
}
