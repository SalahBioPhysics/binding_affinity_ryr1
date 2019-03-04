/*******************************************************************************
 * NAME
 *        iatom - get the index number from conformer name and atom name
 *
 * SYNOPSIS
 *        #include <gdbm.h>
 *        #include <mcce.h>
 *
 *        int iatom(char *conf_name, char *atom_name);
 *
 * DESCRIPTION
 *        The iatom() function returns a unique  index number  predefined in  the
 *        paramter file by the conformer name and atom name.     The keys are not
 *        sensitive  to the leading and ending  spaces of the  conformer name and
 *        atom name. Use "-lgdbm" option to compile these programs, just like use
 *        "-lm" to compile programs calling math functions.
 *
 * RETURN VALUE
 *        The returned integer is the index number of an atom in a conformer.  If
 *        the atom index number is found in parameter files, -1 is returned.
 *
 * SEE ALSO
 *        db_open, db_close, param_get, param_sav
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
 *        Junjun Mao, 06/02/2003
 *******************************************************************************/
#include <string.h>
#include <gdbm.h>
#include <stdlib.h>
#include "mcce.h"

extern GDBM_FILE param_db;    /* defined by db_open() */

int iatom(char *conf_name, char *atom_name)
{  datum pkey, pvalue;
   char key[MAXCHAR_LINE];
   char sbuff[MAXCHAR_LINE];
   int i_atom;

   /* construct the atom key */
   strcpy(key, "IATOM");
   strip(sbuff, conf_name); strcat(key, sbuff);
   strip(sbuff, atom_name); strcat(key, sbuff);
   pkey.dptr = key;
   pkey.dsize = strlen(key);

   /* find the atom in the database */
   if (gdbm_exists(param_db, pkey)) {        /* existing key */
      pvalue = gdbm_fetch(param_db, pkey);
      /* cast generic pointer pvalue.dptr to integer pointer by (int *), then
       * dereference the pointer to an integer by *,   so the return value is
       * an integer */
       i_atom = * (int *) pvalue.dptr;
       free(pvalue.dptr);
      return  i_atom;
   }
   else return -1;                           /* not in the database */
}
