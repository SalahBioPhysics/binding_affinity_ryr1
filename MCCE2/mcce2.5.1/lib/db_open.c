/*******************************************************************************
 * NAME
 *        db_open - open mcce parameter database
 *
 * SYNOPSIS
 *        #include <gdbm.h>
 *        #include <mcce.h>
 *
 *        int db_open();
 *
 * DESCRIPTION
 *        The db_open() function opens a paramter database for fast access to the
 *        parameter entries with the support  of  standard gdbm library.  If this
 *        function is called  while the database is  already open,  it prints out
 *        a warning message and use the existing database.  Use "-lgdbm"   option
 *        to  compile  these programs,  just  like  use "-lm" to compile programs
 *        calling math functions.
 *
 * RETURN VALUE
 *        The returned integer is the status, 0 is success, other numbers are the
 *        error code defined by gdbm.h.
 *
 * SEE ALSO
 *        db_close, iatom, param_get, param_sav
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcce.h"

GDBM_FILE param_db;
char gdbm_file[256];

int db_open()
{
    GDBM_FILE temp_db;
    
    //remove(DUMMY_GDBM);

    strcpy(gdbm_file, DUMMY_GDBM);
    mkstemp(gdbm_file);
    
    /* open the param database */
    temp_db = gdbm_open(gdbm_file, 16*1024, GDBM_WRCREAT, 0666, 0);
    
    if (gdbm_errno == 0)          /* success, accept temp_db */
        param_db = temp_db;
    else if (gdbm_errno == 10)    /* open existing database, use old param_db */
        printf("   db_open(): Warning, database is already open.\n");
    
    return gdbm_errno;
}
