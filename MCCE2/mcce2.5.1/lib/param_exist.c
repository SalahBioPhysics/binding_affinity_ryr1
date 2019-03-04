#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gdbm.h>
#include "mcce.h"

extern GDBM_FILE param_db;

int param_exist(char *key1, char *key2, char *key3)
{
    datum pkey;
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
    if (!gdbm_exists(param_db, pkey)) {
        return 0; /* failure */
    }
    else {
        return 1;  /* success */
    }
}
