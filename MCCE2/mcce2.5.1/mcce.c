#include <stdio.h>
#include "mcce.h"
//#include <mpi.h>

void welcome();

int main(int argc, char *argv[])
{
   /* Welcome */
   welcome();


   /* Do step 0, initialization */
   db_open();
   printf("Step 0. Initialize enviroment\n"); fflush(stdout);
   if (init()) {
      db_close();
      printf("Help message: double check file \"run.prm\" in current directory.\n");
      return USERERR;
   }
   else printf("Step 0 Done.\n\n");


   /* Do step 1, premcce */
   if (env.do_premcce) {
      printf("Step 1. Test and format structral file\n"); fflush(stdout);
      if (premcce()) {db_close(); return USERERR;}
      else printf("Step 1 Done.\n\n");
   }
   else printf("Not doing \"Step 1. Test and format structral file\"\n\n");


   /* Do step 2. rotamers */
   if (env.do_rotamers) {
      printf("Step 2. Make multi side chain conformers\n"); fflush(stdout);
	if (rotamers()) {
		db_close(); return USERERR;
    }
    else printf("Step 2 Done.\n\n");
   }
   else printf("Not doing \"Step 2. Make multi side chain conformers\"\n\n");

   /* Do step 3. energies */
   if (env.do_energies) {
      printf("Step 3. Compute energy lookup table\n"); fflush(stdout);
      if (energies()) {db_close(); return USERERR;}
      else printf("Step 3 Done.\n\n");
   }
   else printf("Not doing \"Step 3. Compute energy lookup table\"\n\n");

   /* Do step 4. Monte Carlo */
   if (env.do_monte) {
      printf("Step 4. Monte Carlo Sampling\n"); fflush(stdout);
      if (!env.monte_adv_opt) {
      if (monte()) {db_close(); return USERERR;}
           else printf("Step 4 Done.\n\n");
       }
       else {
           if (monte2()) {db_close(); return USERERR;}
           else printf("Step 4 Done.\n\n");
       }
   }
   else printf("Not doing \"Step 4. Monte Carlo Sampling\"\n\n");


   db_close();
   return 0;
}

void welcome()
{  printf("===========================================================\n");
   printf("<<< MCCE Multi-Conformation Continuum Electrostatics >>>   \n");
   printf(" Marilyn Gunner's Lab at City College of New York, 2005    \n");
   printf("-----------------------------------------------------------\n");
   printf("Version:        2.5.1                                      \n");
   printf("MCCE Home Page: http://www.sci.ccny.cuny.edu/~gunner/mcce  \n");
   printf("Support:        gunner@sci.ccny.cuny.edu                   \n");
   printf("Developed by:   Junjun Mao, Yifan Song, Marilyn Gunner     \n");
   printf("Reference MCCE: If you publish data calculated with MCCE,  \n");
   printf("                you need to cite papers suggested in MCCE  \n");
   printf("                Home Page.                                 \n");
   printf("===========================================================\n\n");
   printf("Last Updates:                                              \n");
//   printf("   05/10, Step2, use relative energy to decide vdw clash   \n");
   printf("   04/10, Added Pascal Comte's GA/SA algorithm for sidechain\n");
   printf("          packing and sampling.     			      \n");
   printf("   04/10, Added APBS electrostatic calculation option      \n");
   printf("   08/08, Total pairwise is set to be 999.0 if vdw is 999.0\n");
   printf("   08/02, step 3, zero interaction bug fixed               \n");
   printf("   07/23, step 1 changes the policy on AltLoc backbone atoms\n");
   printf("          now it is warning instead of fatal error.        \n");
   printf("   06/22, step 1 recognizes FME and ACE as NTR caps        \n");
   printf("   05/27, step 2 reduces rotation steps for surface residues\n");
   printf("          The run.prm.default uses 0.25 SAS to decide suface\n");
   printf("          residues.                                        \n");
   printf("   05/04, Program doesn't stop at duplicate parameters,    \n");
   printf("          it issues a warning message. This fixes the deadlock\n");
   printf("          caused by multiple new cofactors                 \n");
   printf("   04/24, Energy lookup table compressed, size is reduced  \n");
   printf("          to about 1/12                                    \n");
   printf("   04/24, Opp files are compressed. Temporary files are    \n");
   printf("          stored under local directory.                    \n");
   printf("   04/14, Step 2 rotamer pruning by pairwise corrected.    \n");
   printf("   04/08, Step 2 repacking now ses better H bnd correction.\n");
   printf("   01/26, Step 3 fixed incorrect errno message.            \n");
   printf("   01/14, Step 3 dielectric boundary condition improved.   \n");
   printf("   12/07, Step 3 dielectric boundary condition is corrected.\n");
   printf("   11/17, pKa fitting writes \"pKa titration curve too sharp\"\n");
   printf("          when occ jumps from 0.015 to 0.985.              \n");
   printf("   11/17, Exposed conformer is marked as E in history.     \n");
   printf("   11/02, Optional quick step 3 included.                  \n");
   printf("   11/02, Pruning fuction is added at the end of step 2.   \n");
   printf("          Parameter lines were added to run.prm.           \n");
   printf("   10/29, Temporary gdbm files have file names.            \n");
   printf("   10/11, Use biggest rxn energy from the last 3 delphi    \n");
   printf("          focusing runs.                                   \n");
   printf("===========================================================\n\n");
   fflush(stdout);

   return;
}
