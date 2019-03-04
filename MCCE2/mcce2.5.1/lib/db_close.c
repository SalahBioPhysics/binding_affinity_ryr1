/*******************************************************************************
 * NAME
 *        db_close - close mcce parameter database
 *
 * SYNOPSIS
 *        #include <gdbm.h>
 *        #include <mcce.h>
 *
 *        int db_close();
 *
 * DESCRIPTION
 *        The db_close()  function clears and closes the mcce parameter database.
 *        You only need  to call  this function once when  you no longer need the
 *        parameter database. Use "-lgdbm" option to compile these programs, just
 *        like use "-lm" to compile programs calling math functions.
 *
 * RETURN VALUE
 *        The returned integer is 0.
 *
 * SEE ALSO
 *        db_open, iatom, param_get, param_sav
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
 *          param_sav("IATOM",    "ASP0", "CG", &i, sizeof(int));
 *          param_sav("PKA",      "ASP-", "",   &pKa, sizeof(float));
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
 * AUTHOR
 *        Junjun Mao, 06/02/2003
 *******************************************************************************/

#include <gdbm.h>
#include "mcce.h"

extern GDBM_FILE param_db; /* defined by db_open() */
extern char gdbm_file[256];

int db_close()
{  gdbm_close(param_db);   /* close the param database */
   remove(gdbm_file);     /* delete temporary file */
   return 0;
}
