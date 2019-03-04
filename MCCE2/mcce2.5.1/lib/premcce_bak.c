/**************************************************************
* Premcce
*
* This module is the step 0 of mcce. It does these:
* 1. Identify NTR and CTR;
* 2. Rename atoms/residues;
* 3. Identify clashes;
* 4. Create conformers for atoms with altLoc;
* 5. Try loading to the data sturcture;
**************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mcce.h"
#include <errno.h>
#include <unistd.h>

#define MAXLINES_NAME 1000
#define BOND_THR   2.2
#define BOND_THR_H 1.8

char rules_from[MAXLINES_NAME][15];
char rules_to[MAXLINES_NAME][15];
int c;

FILE *premcce_rename(FILE *fp_in);
FILE *premcce_terminals(FILE *fp_in);
FILE *premcce_confname(FILE *fp);
int premcce_match(char *sfrom, char *pattern);
int premcce_replace(char *sto, char *pattern);
int premcce_hvatoms(PROT prot);
int premcce_clash(PROT prot);
int create_param(FILE *pdb_fp, int k_line);

int premcce()
{  FILE *fp1, *fp2, *fp;
    char sbuff[MAXCHAR_LINE];
    char *strp;
    PROT prot;
    
    /* read pdb file into stream */
    printf("   Read pdb file \"%s\"...\n", env.inpdb);   fflush(stdout);
    if ((fp1=fopen(env.inpdb, "r")) == NULL) {
        printf("   FATAL: premcce(): \"failed opening file \"%s\"\"\n", env.inpdb);
        return USERERR;
    }
    fp2 = tmpfile();
    while ((c=fgetc(fp1)) != EOF) fputc(c, fp2);
    fputc(EOF, fp2);
    fclose(fp1);
    rewind(fp2);
    printf("   Done\n");
    fflush(stdout);
    
    /* read rename rules */
    printf("   Read name rule file \"%s\"...\n", env.rename_rules);   fflush(stdout);
    if ((fp1=fopen(env.rename_rules, "r")) == NULL) {
        printf("   WARNNING: premcce(): \"failed opening file \"%s\"\"\n", env.rename_rules);
        printf("                        Skipping rename process\n");
        fp1 = fp2;
    }
    else {
        c = 0;
        while (fgets(sbuff, sizeof(sbuff), fp1)) {
            if (!(strp = strtok(sbuff, " \n"))) continue;
            if (strp[0] == '#') continue;
            if (strlen(strp) != 14) {
                printf("   FATAL: premcce():\"\"%s\" is not a valid naming rule.\"\n", strp);
                return USERERR;
            }
            strcpy(rules_from[c], strp); rules_from[c][14] = '\0';
            
            if (!(strp = strtok(NULL, " \n"))) {
                printf("   FATAL: premcce(): \"missing replace-to string.\"\n");
                return USERERR;
            }
            if (strlen(strp) != 14) {
                printf("   FATAL: premcce(): \"%s\" is not a valid naming rule.\"\n", strp);
                return USERERR;
            }
            strcpy(rules_to[c], strp); rules_to[c][14] = '\0';
            
            c++;
        }
        fclose(fp1);
        printf("   Done\n");
        fflush(stdout);
        /* pdb now in fp2 */
        
        /* rename residue/atom  based on name rules */
        printf("   Rename residue and atom names...\n");   fflush(stdout);
        fp1 = premcce_rename(fp2);
        fclose(fp2);
    }
    printf("   Done\n");
    fflush(stdout);
    /* pdb now in fp1 */
    
    /* Identify NTR and CTR */
    printf("   Identify NTR and CTR...\n");   fflush(stdout);
    fp2 = premcce_terminals(fp1);
    fclose(fp1); rewind(fp2);
    printf("   Done\n");
    fflush(stdout);
    /* pdb now in fp2 */
    
    
    /* label conformer by altLoc */
    printf("   Label backbone, sidechain and altLoc conformers...\n");
    remove(env.new_tpl);
    while (!(fp1 = premcce_confname(fp2))) {
        if (env.err_msg >= 0) {
            printf("   Creating temporary parameter file for unrecognized residue...\n");
            if (create_param(fp2, env.err_msg)) printf("   parameter file %s not made.\n",env.new_tpl),exit(-1);
            load_param(env.new_tpl);
            rewind(fp2);
            printf("   Trying labelling again...\n");
        }
        else {
            printf("   STOP: fatal errors reported by premcce_confname()\n");
            return USERERR;
        }
    }
    fclose(fp2);
    printf("   Done\n");
    fflush(stdout);
    
    /*
    while ((c=fgetc(fp1)) != EOF) {
        if (fputc(c, stdout)==EOF)
            printf("ERROR\n");
    }
    fputc(EOF, stdout); rewind(fp1);
    */
    
    /*
    while ((c=fgetc(fp1)) != EOF) fputc(c, stdout);
    fputc(EOF, stdout);
    */
    
    
    /* Load to data structure */
    printf("   Load pdb lines into data structure...\n");   fflush(stdout);
    prot = load_pdb(fp1);
    fclose(fp1);
    if (prot.n_res == 0) {
        printf("   There are errors in pdb file, quiting ...\n");
        return USERERR;
    }
    printf("   Done\n");
    fflush(stdout);
    
    
    /* Missing heavy atoms, complete altLoc conformer */
    printf("   Check missing heavy atoms and complete altLoc conformers...\n");   fflush(stdout);
    if (premcce_hvatoms(prot)) {
        printf("   Missing heavy atoms detected.\n");
        printf("   Ignore warning messages if they are in the terminal residues\n");
    }
    printf("   Done\n");
    fflush(stdout);
    
    
    /* Ligand finding and clash checking */
    printf("   Find distance clash (<%.3f)...\n", env.clash_distance);   fflush(stdout);
    if (premcce_clash(prot) == 0)
        printf("   No clash found.\n");
    printf("   Done\n");
    fflush(stdout);
    
    printf("   Make connectivity network ...\n");
    if (get_connect12(prot)) {
        printf("   Errors were detected when making connectivity.\n");
        printf("   Done.\n");
    }
    else printf("   Passed.\n");
    
    /* write out pdb file */
    if (!(fp = fopen(STEP1_OUT, "w"))) {
        printf("   FATAL: premcce(): \"Can't open file %s to write\n", STEP1_OUT);
        return USERERR;
    }
    write_pdb(fp, prot);
    fclose(fp);
    
    printf("   Output file is \"%s\", which will be the input of step 2.\n", STEP1_OUT);   fflush(stdout);
    fflush(stdout);
    return 0;
}

FILE *premcce_rename(FILE *fp_in)
{  char sbuff[MAXCHAR_LINE];
    char sfrom[15];
    char sto[15];
    int i;
    FILE *fp_out;
    
    fp_out = tmpfile();
    
    while (fgets(sbuff, sizeof(sbuff), fp_in)) {
        for (i=0; i<c; i++) {
            strncpy(sfrom, sbuff+12, 14); sfrom[14] = '\0';
            if ((!strncmp(sbuff, "ATOM  ", 6) || !strncmp(sbuff, "HETATM", 6)) &&
            premcce_match(sfrom, rules_from[i])) {
                strcpy(sto, sfrom);                 /* a copy */
                premcce_replace(sto, rules_to[i]);  /* update copy */
                printf("   Renaming \"%s\" to \"%s\"\n", sfrom, sto);
                strncpy(sbuff+12, sto, 14);
            }
        }
        fputs(sbuff, fp_out);
    }
    fputc(EOF, fp_out);
    rewind(fp_out); rewind(fp_in);
    
    return fp_out;
}

int premcce_match(char *sfrom, char *pattern)
{  int i;
    int match;
    
    match = 1;
    //printf("===%s===%s===\n", sfrom, pattern);
    for (i=0; i<14; i++) {
        if (pattern[i] == '*' || pattern[i] == sfrom[i]) continue;
        else {
            match = 0;
            break;
        }
    }
    return match;
}

int premcce_replace(char *sto, char *pattern)
{  int i;
    
    for (i=0; i<14; i++) {
        if (pattern[i] == '*') continue;
        else sto[i] = pattern[i];
    }
    return 0;
}

FILE *premcce_terminals(FILE *fp_in)
{  ATOM *Ns, *Cs;
    int  i_Ns, n_Ns, i_Cs, n_Cs;
    char sbuff[MAXCHAR_LINE];
    char stemp[MAXCHAR_LINE];
    char NTR_atoms[] = "1H  2H  3H   N   CA ";
    char CTR_atoms[] = " C   O   OXT";
    char ter;
    FILE *fp1, *fp2, *fp_out;
    int resSeq;
    int i;
    
    /* collect all N and C atoms */
    n_Ns = 0; n_Cs = 0;
    while (fgets(sbuff, sizeof(sbuff), fp_in)) {
        if (!strncmp(sbuff, "ATOM  ", 6)) {
            if (!strncmp(sbuff+12, " N  ", 4)) n_Ns++;
            else if (!strncmp(sbuff+12, " C  ", 4)) n_Cs++;
        }
    }
    rewind(fp_in);
    if ((Ns = (ATOM *) malloc(n_Ns*sizeof(ATOM))) == NULL) {
        printf("   FATAL: premcce_terminals(): \"Memory erreor.\"\n");
        return NULL;
    }
    if ((Cs = (ATOM *) malloc(n_Cs*sizeof(ATOM))) == NULL) {
        printf("   FATAL: premcce_terminals(): \"Memory erreor.\"\n");
        return NULL;
    }
    i_Ns = 0; i_Cs = 0;
    while (fgets(sbuff, sizeof(sbuff), fp_in)) {
        if (!strncmp(sbuff, "ATOM  ", 6)) {
            if (!strncmp(sbuff+12, " N  ", 4)) {
                Ns[i_Ns] = pdbline2atom(sbuff);
                i_Ns++;
            }
            else if (!strncmp(sbuff+12, " C  ", 4)) {
                Cs[i_Cs] = pdbline2atom(sbuff);
                i_Cs++;
            }
        }
    }
    rewind(fp_in);
    
    fp1 = tmpfile();
    while ((i=fgetc(fp_in)) != EOF) fputc(i, fp1);
    fputc(EOF, fp1); rewind(fp1); rewind(fp_in);
    
    /* Do Ntrs by distance checking */
    for (i_Ns=0; i_Ns<n_Ns; i_Ns++) {
        ter = 1;
        for (i_Cs=0; i_Cs<n_Cs; i_Cs++) {
            if (ddvv(Ns[i_Ns].xyz, Cs[i_Cs].xyz) < 3.0) {
                ter = 0;
                break;
            }
        }
        if (ter) {
            printf("   Labeling \"%s %c%4d\" as NTR\n", Ns[i_Ns].resName,
            Ns[i_Ns].chainID,
            Ns[i_Ns].resSeq);
            fp2 = tmpfile();
            while (fgets(sbuff, sizeof(sbuff), fp1)) {
                if (!strncmp(sbuff, "ATOM  ", 6) &&
                    !strncmp(sbuff+17, Ns[i_Ns].resName, 3) &&
                sbuff[21] == Ns[i_Ns].chainID) {
                    strncpy(stemp, sbuff+22, 4); stemp[4] = '\0';
                    resSeq = atoi(stemp);
                    strncpy(stemp, sbuff+12, 4); stemp[4] = '\0';
                    if(strstr(NTR_atoms, stemp) && resSeq==Ns[i_Ns].resSeq)
                        strncpy(sbuff+17, "NTR", 3);
                }
                fputs(sbuff, fp2);
            }
            fputc(EOF, fp2); rewind(fp2);
            fclose(fp1); fp1 = tmpfile();
            while ((c=fgetc(fp2)) != EOF) fputc(c, fp1);
            fputc(EOF, fp1); rewind(fp1); fclose(fp2);
        }
    }
    
    /* Do Ctrs by distance checking */
    for (i_Cs=0; i_Cs<n_Cs; i_Cs++) {
        ter = 1;
        for (i_Ns=0; i_Ns<n_Ns; i_Ns++) {
            if (ddvv(Ns[i_Ns].xyz, Cs[i_Cs].xyz) < 3.0) {
                ter = 0;
                break;
            }
        }
        if (ter) {
            printf("   Labeling \"%s %c%4d\" as CTR\n", Cs[i_Cs].resName,
            Cs[i_Cs].chainID,
            Cs[i_Cs].resSeq);
            fp2 = tmpfile();
            while (fgets(sbuff, sizeof(sbuff), fp1)) {
                if (!strncmp(sbuff, "ATOM  ", 6) &&
                    !strncmp(sbuff+17, Cs[i_Cs].resName, 3) &&
                sbuff[21] == Cs[i_Cs].chainID) {
                    strncpy(stemp, sbuff+22, 4); stemp[4] = '\0';
                    resSeq = atoi(stemp);
                    strncpy(stemp, sbuff+12, 4); stemp[4] = '\0';
                    if(strstr(CTR_atoms, stemp) && resSeq==Cs[i_Cs].resSeq)
                        strncpy(sbuff+17, "CTR", 3);
                }
                fputs(sbuff, fp2);
            }
            fputc(EOF, fp2); rewind(fp2);
            fclose(fp1); fp1 = tmpfile();
            while ((c=fgetc(fp2)) != EOF) fputc(c, fp1);
            fputc(EOF, fp1); rewind(fp1); fclose(fp2);
        }
    }
    
    fp_out = tmpfile();
    while ((c=fgetc(fp1)) != EOF) fputc(c, fp_out);
    fputc(EOF, fp_out); rewind(fp_out);
    
    rewind(fp_in);
    
    free(Ns); free(Cs);
    return fp_out;
}


FILE *premcce_confname(FILE *fp) {
    FILE *fp2;
    STRINGS pdb, checklist, residue;
    STRINGS conflist;
    char  sbuffer[MAXCHAR_LINE];   /* line buffer */
    char  confName[6], confID[3];
    int   fatal_err=0;
    int  i_check, i_conf, k_atom, iatom,i_pdb;
    ATOM atom;
    int i;
    
    memset(&pdb,        0, sizeof(STRINGS));
    memset(&checklist,  0, sizeof(STRINGS));
    memset(&residue,    0, sizeof(STRINGS));
    env.err_msg = -1;
    
    /* jmao */
    pdb.n = 0;
    
    while ( fgets(sbuffer, MAXCHAR_LINE, fp) ) pdb.n++;
    if (!(pdb.strings = (char **) malloc(pdb.n*sizeof(char *)))) {
        printf("   FATAL: premcce_confname(): memory error\n");
        return NULL;
    }
    
    rewind(fp);
    i = 0;
    while ( fgets(sbuffer, MAXCHAR_LINE, fp) ) {
        if (sbuffer[strlen(sbuffer) - 1] == '\n') sbuffer[strlen(sbuffer) - 1] = '\0';      /* delete new line character at the end of the line */
        if (!strncmp("ATOM  ", sbuffer, 6) || !strncmp("HETATM", sbuffer, 6)) {
            while (strlen(sbuffer)<84) strcat(sbuffer, " ");      /* Extend length of pdb line to 84 characters */
        }
        
        pdb.strings[i] = (char *)malloc(MAXCHAR_LINE*sizeof(char));
        memset(pdb.strings[i], 0, MAXCHAR_LINE*sizeof(char));
        strcpy(pdb.strings[i], sbuffer);
        i++;
        
        if (!strncmp("ATOM  ", sbuffer, 6) || !strncmp("HETATM", sbuffer, 6)) {
            
            atom = pdbline2atom(sbuffer);
            
            /* error checkings for conflist parameter */
            if ( param_get("CONFLIST",atom.resName, "",&conflist) ) {
                printf("   Error! premcce_confname(): Can't get conformer list of this residue %s\n", atom.resName);
                env.err_msg = i-1;
                return NULL;
            }
            
            if ( strcmp(conflist.strings[0]+3, "BK" ) )  {
                printf("   FATAL: first conformer in conformer list of %s is not named as \"BK\"\n", atom.resName);
                return NULL;
            }
            
            if ( strncmp(sbuffer+80,"  ",2) ) continue;         /* Conformer is given. */
            
            if ( !param_get("IATOM", conflist.strings[0], atom.name, &iatom) ) {
                strncpy(pdb.strings[i-1]+80, "BKOr", 4);
                strncpy(pdb.strings[i-1]+27, "00", 2);
                continue;
            }
            
            if (conflist.n>1) {
                if ( !param_get("IATOM", conflist.strings[1], atom.name, &iatom) ) continue;
                
                checklist.n++;
                checklist.strings = realloc(checklist.strings, checklist.n*sizeof(void *));
                checklist.strings[checklist.n-1] = pdb.strings[i-1];
            }
        }
    }
    rewind(fp);
    
    while (checklist.n) {
        atom = pdbline2atom(checklist.strings[0]);
        
        /* Collect all atoms belonging to the residue of first line in checklist */
        for (i_pdb = 0; i_pdb<pdb.n; i_pdb++) {
            if (
                !strncmp(pdb.strings[i_pdb]+16, checklist.strings[0]+16, 4) &&
            !strncmp(pdb.strings[i_pdb]+21, checklist.strings[0]+21, 8) &&
            !strncmp(pdb.strings[i_pdb]+80, checklist.strings[0]+80, 2) ) {
                residue.n++;
                residue.strings = realloc(residue.strings, residue.n*sizeof(void *));
                residue.strings[residue.n-1] = pdb.strings[i_pdb];
            }
        }
        
        /* Remove all atoms in this checklist */
        for (i_check = 0; i_check<checklist.n; i_check++) {
            if (
                !strncmp(checklist.strings[i_check]+16, residue.strings[0]+16, 4) &&
            !strncmp(checklist.strings[i_check]+21, residue.strings[0]+21, 8) &&
            !strncmp(checklist.strings[i_check]+80, residue.strings[0]+80, 2)
            ) {
                checklist.n--;
                memmove(&checklist.strings[i_check],&checklist.strings[i_check+1],(checklist.n-i_check)*sizeof(void *));
                i_check--;
            }
        }
        
        /* Assign confomer name and conformer ID to the side chain of this residue */
        
        /* looping over all side chain confomers in conformer list */
        param_get("CONFLIST",atom.resName, "",&conflist);
        for (i_conf=1; i_conf<conflist.n; i_conf++) {
            /* Check if all atoms can go into this conformer */
            for (k_atom=0; k_atom<residue.n; k_atom++) {
                atom = pdbline2atom(residue.strings[k_atom]);
                /* If can't find IATOM slot for one atom, break looping all atoms and go on to use next conformer name */
                if ( param_get("IATOM", conflist.strings[i_conf], atom.name, &iatom) ) {
                    break;
                }
            }
            if (k_atom == residue.n) {
                break;              /* If looping over all atoms is finished without broken, means all atoms can be found by current conformer name, then break looping over conformer list */
            }
        }
        
        if ( i_conf == conflist.n ) {
            printf("   Error! Can't find one side chain to fit all the atoms, residue %s %c%4d\n",atom.resName, atom.chainID, atom.resSeq);
            printf("          Atoms found from PDB file are: \n");
            for (k_atom=0; k_atom<residue.n; k_atom++) {
                atom = pdbline2atom(residue.strings[k_atom]);
                printf("          %s\n",atom.name);
            }
            fatal_err++;
            
            strcpy(confName,"?????");
            strcpy(confID, "??");
        }
        else {
            strcpy(confName,conflist.strings[i_conf]);
            
            if (!strcmp(atom.confID,"  ")) {
                if (atom.altLoc != ' ') {
                    confID[0] = '_';
                    confID[1] = atom.altLoc;
                }
                else {
                    sprintf(confID, "%02d", i_conf);
                }
            }
            else {
                strcpy(confID, atom.confID);
            }
        }
        
        for (k_atom=0; k_atom<residue.n; k_atom++) {
            strncpy(residue.strings[k_atom]+80, confName+3, 2);
            strncpy(residue.strings[k_atom]+82, "Or", 2);
            strncpy(residue.strings[k_atom]+27, confID, 2);
        }
        strncpy(sbuffer,residue.strings[0]+16,10);sbuffer[10]='\0';
        free(residue.strings);
        memset(&residue,0,sizeof(STRINGS));
    }
    
    for (i_pdb = 0; i_pdb<pdb.n; i_pdb++) {
        if ( strncmp(pdb.strings[i_pdb]+80,"  ",2) ) continue;         /* Conformer is given. */
        atom = pdbline2atom(pdb.strings[i_pdb]);
        param_get("CONFLIST",atom.resName, "",&conflist);
        if ( !param_get("IATOM", conflist.strings[1], atom.name, &iatom) ) {
            strncpy(pdb.strings[i_pdb]+80, conflist.strings[1]+3, 2);
            strncpy(pdb.strings[i_pdb]+82, "Or", 2);
            strncpy(pdb.strings[i_pdb]+27, "01", 2);
        }
        else {
            fatal_err++;
        }
    }
    
    
    if (fatal_err) return NULL;
    
    fp2 = tmpfile();
    for (i_pdb=0; i_pdb<pdb.n; i_pdb++) {
        fprintf(fp2, "%s\n",pdb.strings[i_pdb]);
    }
    rewind(fp2);
    
    free_strings(&pdb); /* jmao */
    return fp2;
}

int premcce_hvatoms(PROT prot)
{  int kr, kc, ka;
    char Missing, Alt;
    char sbuff[MAXCHAR_LINE], siatom[MAXCHAR_LINE];
    
    /* altLoc atoms */
    Alt = 0;
    for (kr=0; kr<prot.n_res; kr++) {
        for (kc=0; kc<prot.res[kr].n_conf; kc++) {
            for (ka=0; ka<prot.res[kr].conf[kc].n_atom; ka++) {
                if (prot.res[kr].conf[kc].atom[ka].altLoc!=' ') {
                    Alt = 1;
                    break;
                }
            }
            if (Alt && kc > 1) {
                for (ka=0; ka<prot.res[kr].conf[kc].n_atom; ka++) {
                    if (!prot.res[kr].conf[kc].atom[ka].on) {
                        prot.res[kr].conf[kc].atom[ka] = prot.res[kr].conf[1].atom[ka];
                    }
                }
                Alt = 0;
            }
        }
    }
    
    /* Missing atoms */
    Missing = 0;
    for (kr=0; kr<prot.n_res; kr++) {
        for (kc=0; kc<prot.res[kr].n_conf; kc++) {
            for (ka=0; ka<prot.res[kr].conf[kc].n_atom; ka++) {
                if (!prot.res[kr].conf[kc].atom[ka].on) {
                    sprintf(siatom, "%3d", ka);
                    if (param_get("ATOMNAME", prot.res[kr].conf[kc].confName, siatom, sbuff)) {
                        printf("   Missing ATOMNAME records for slot \"%d\" of conformer %s.\n",
                        ka, prot.res[kr].conf[kc].confName);
                        Missing++;
                    }
                    else if (sbuff[1] != 'H') {
                        printf("   Missing heavy atom %s of conf %s in \"%s %c %3d\".\n",
                        sbuff, prot.res[kr].conf[kc].confName, prot.res[kr].resName, prot.res[kr].chainID, prot.res[kr].resSeq);
                        Missing++;
                    }
                }
            }
        }
    }
    
    return Missing;
}

int premcce_clash(PROT prot)
{  int kr, kc, ka, ir, ic, ia;
    float limit = env.clash_distance * env.clash_distance;
    float dd;
    int n=0;
    
    for (kr=0; kr<prot.n_res; kr++) {
        for (ir=kr+1; ir<prot.n_res; ir++) {
            for (kc=0; kc<prot.res[kr].n_conf; kc++) {
                for (ka=0; ka<prot.res[kr].conf[kc].n_atom; ka++) {
                    if (prot.res[kr].conf[kc].atom[ka].on == 0 ||
                        prot.res[kr].conf[kc].atom[ka].name[1] == 'H') continue;
                    for (ic=0; ic<prot.res[ir].n_conf; ic++) {
                        for (ia=0; ia<prot.res[ir].conf[ic].n_atom; ia++) {
                            if (prot.res[ir].conf[ic].atom[ia].on == 0 ||
                                prot.res[ir].conf[ic].atom[ia].name[1] == 'H') continue;
                            if ((dd=ddvv(prot.res[kr].conf[kc].atom[ka].xyz, prot.res[ir].conf[ic].atom[ia].xyz)) < limit) {
                                /* exclude normal bonds */
                                if ((!strcmp(prot.res[kr].conf[kc].atom[ka].name, " C  ") && !strcmp(prot.res[ir].conf[ic].atom[ia].name, " N  ")) ||
                                    (!strcmp(prot.res[kr].conf[kc].atom[ka].name, " N  ") && !strcmp(prot.res[ir].conf[ic].atom[ia].name, " C  ")) ||
                                (!strcmp(prot.res[kr].conf[kc].atom[ka].name, " CA ") && !strcmp(prot.res[ir].conf[ic].atom[ia].name, " C  ")) ||
                                (!strcmp(prot.res[kr].conf[kc].atom[ka].name, " C  ") && !strcmp(prot.res[ir].conf[ic].atom[ia].name, " CA ")))
                                continue;
                                else {
                                    printf("   d=%5.2f: \"%s %s %c %3d\" to \"%s %s %c %3d\"\n", sqrt(dd),
                                    prot.res[kr].conf[kc].atom[ka].name,
                                    prot.res[kr].conf[kc].atom[ka].resName,
                                    prot.res[kr].conf[kc].atom[ka].chainID,
                                    prot.res[kr].conf[kc].atom[ka].resSeq,
                                    prot.res[ir].conf[ic].atom[ia].name,
                                    prot.res[ir].conf[ic].atom[ia].resName,
                                    prot.res[kr].conf[kc].atom[ka].chainID,
                                    prot.res[kr].conf[kc].atom[ka].resSeq);
                                    n++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    return n;
}


int create_param(FILE *pdb_fp, int k_line) {
    CONF conf;
    int cnt, ichar, i_atom, j_atom;
    char sbuffer[MAXCHAR_LINE], line[MAXCHAR_LINE];
    FILE *param_fp;
    
    /* Collect a list of atoms from the same residue. */
    rewind(pdb_fp);
    cnt = -1;
    while (cnt < k_line && fgets(sbuffer, MAXCHAR_LINE, pdb_fp) ) {
        cnt++;
    }
    if (cnt != k_line) {
        printf(" Error! create_param(): Error in reaching defined line number. line counter = %d, defined line number = %d\n",cnt,k_line);
        return USERERR;
    }
    if (strncmp("ATOM  ", sbuffer, 6) && strncmp("HETATM", sbuffer, 6)) {
        printf(" Error! create_param(): input line is not an atom line. \n");
        return USERERR;
    }
    
    strcpy(line, sbuffer);
    if (line[strlen(line)-1] == '\n') line[strlen(line)-1] = '\0';
    while (strlen(line)<84) strcat(line, " ");
    
    memset(&conf,0,sizeof(CONF));
    rewind(pdb_fp);
    while ( fgets(sbuffer, MAXCHAR_LINE, pdb_fp) ) {
        if (strncmp("ATOM  ", sbuffer, 6) && strncmp("HETATM", sbuffer, 6)) continue;
        if (sbuffer[strlen(sbuffer)-1] == '\n') sbuffer[strlen(sbuffer)-1] = '\0';
        while (strlen(sbuffer)<84) strcat(sbuffer, " ");
        if (!strncmp(sbuffer+16, line+16, 4) &&
            !strncmp(sbuffer+21, line+21, 8) &&
        !strncmp(sbuffer+80, line+80, 2) ) {
            conf.n_atom++;
            conf.atom = realloc(conf.atom, conf.n_atom*sizeof(ATOM));
            memset(&conf.atom[conf.n_atom-1],0,sizeof(ATOM));
            conf.atom[conf.n_atom-1] = pdbline2atom(sbuffer);
        }
    }
    
    for (ichar=0;ichar<3;ichar++) {
        if (conf.atom[0].resName[ichar] == ' ') {
            printf("   Error! create_param(): A space found in residue name %s\n",conf.atom[0].resName);
            return USERERR;
        }
    }
    
    
    param_fp = fopen(env.new_tpl,"a");
    fprintf(param_fp,"### This is a temporary parameter file made for residue %s ###\n", conf.atom[0].resName);
    fprintf(param_fp,"### Make sure that all the parameters are verified before using this file as a global parameter file ###\n\n");
    
    /* CONFLSIT */
    strcpy(sbuffer, "CONFLIST ");
    strcat(sbuffer, conf.atom[0].resName);
    strcat(sbuffer, "        ");
    
    strcat(sbuffer, conf.atom[0].resName);
    strcat(sbuffer, "BK ");
    
    fprintf(param_fp,"%s\n\n",sbuffer);
    
    /* NATOM */
    strcpy(sbuffer, "NATOM    ");
    strcat(sbuffer, conf.atom[0].resName);
    strcat(sbuffer, "BK      ");
    sprintf(sbuffer+20, "%d", conf.n_atom);
    fprintf(param_fp,"%s\n\n",sbuffer);
    
    /* IATOM */
    for (i_atom=0; i_atom<conf.n_atom; i_atom++) {
        strcpy(sbuffer, "IATOM    ");
        strcat(sbuffer, conf.atom[0].resName);
        strcat(sbuffer, "BK ");
        strcat(sbuffer, conf.atom[i_atom].name);
        sprintf(sbuffer+19, " %4d", i_atom);
        fprintf(param_fp,"%s\n",sbuffer);
    }
    fprintf(param_fp,"\n");
    
    /* ATOMNAME */
    for (i_atom=0; i_atom<conf.n_atom; i_atom++) {
        strcpy(sbuffer, "ATOMNAME ");
        strcat(sbuffer, conf.atom[0].resName);
        sprintf(sbuffer+12, "BK %4d", i_atom);
        strcat(sbuffer, " ");
        strcat(sbuffer, conf.atom[i_atom].name);
        fprintf(param_fp,"%s\n",sbuffer);
    }
    fprintf(param_fp,"\n");
    
    /* CONNECT */
    for (i_atom=0; i_atom<conf.n_atom; i_atom++) {
        strcpy(sbuffer, "CONNECT  ");
        strcat(sbuffer, conf.atom[0].resName);
        strcat(sbuffer, "BK ");
        strcat(sbuffer, conf.atom[i_atom].name);
        strcat(sbuffer, " ion      ");
        
        for (j_atom=0; j_atom<conf.n_atom; j_atom++) {
            if (j_atom == i_atom) continue;
            if (conf.atom[i_atom].name[1] == 'H' || conf.atom[j_atom].name[1] == 'H') {
                if (dvv(conf.atom[i_atom].xyz, conf.atom[j_atom].xyz) > BOND_THR_H) continue;
            }
            else {
                if (dvv(conf.atom[i_atom].xyz, conf.atom[j_atom].xyz) > BOND_THR) continue;
            }
            strcat(sbuffer, "  0   ");
            strcat(sbuffer, conf.atom[j_atom].name);
        }
        fprintf(param_fp,"%s\n",sbuffer);
    }
    fprintf(param_fp,"\n");
    
    fclose(param_fp);
    return 0;
}
