#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main()
{  FILE *fp;
   float potential[70][70][70];
   char sbuff[128];

  
   fp = fopen("fort.10", "w");
   fprintf(fp, "gsize=%d\n", 65);
   fprintf(fp, "scale=%.2f\n", 0.5);
   fprintf(fp, "in(unpdb,file=\"fort.13\")\n");
   fprintf(fp, "indi=%.1f\n", 4.0);
   fprintf(fp, "exdi=%.1f\n", 80.0);
   fprintf(fp, "ionrad=%.1f\n", 1.4);
   fprintf(fp, "salt=%.2f\n", 0);
   fprintf(fp, "bndcon=2\n");
   fprintf(fp, "center(777, 777, 0)\n");
   fprintf(fp, "out(phi,file=\"run01.phi\")\n");
   fprintf(fp, "site(a,c,p)\n");
   fprintf(fp, "energy(g,s)\n");
   fclose(fp);

   QDIFF_(potential);
   
   fp = fopen("fort.10", "w");
   fprintf(fp, "gsize=%d\n", 65);
   fprintf(fp, "scale=%.2f\n", 0.5);
   fprintf(fp, "in(unpdb,file=\"fort.13\")\n");
   fprintf(fp, "indi=%.1f\n", 4.0);
   fprintf(fp, "exdi=%.1f\n", 80.0);
   fprintf(fp, "ionrad=%.1f\n", 1.4);
   fprintf(fp, "salt=%.2f\n", 0);
   fprintf(fp, "bndcon=3\n");
   fprintf(fp, "center(777, 777, 0)\n");
   fprintf(fp, "out(phi,file=\"run01.phi\")\n");
   fprintf(fp, "site(a,c,p)\n");
   fprintf(fp, "energy(g,s)\n");
   fclose(fp);

   QDIFF_(potential);
     
   fp = fopen("fort.10", "w");
   fprintf(fp, "gsize=%d\n", 65);
   fprintf(fp, "scale=%.2f\n", 0.5);
   fprintf(fp, "in(unpdb,file=\"fort.13\")\n");
   fprintf(fp, "indi=%.1f\n", 4.0);
   fprintf(fp, "exdi=%.1f\n", 80.0);
   fprintf(fp, "ionrad=%.1f\n", 1.4);
   fprintf(fp, "salt=%.2f\n", 0);
   fprintf(fp, "bndcon=3\n");
   fprintf(fp, "center(777, 777, 0)\n");
   fprintf(fp, "out(phi,file=\"run01.phi\")\n");
   fprintf(fp, "site(a,c,p)\n");
   fprintf(fp, "energy(g,s)\n");
   fclose(fp);

   
   QDIFF_(potential);

   return 0;
}
