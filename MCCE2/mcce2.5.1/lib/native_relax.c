#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mcce.h"

#define  PRINT_THR        9
#define  contact2        16
#define  CONVERGED_R     1e-6

extern typedef struct RELAX_STRUCT{
    int    i_relax;
    ATOM   *atom_p;
    
    VECTOR r;
    VECTOR r_p;       /* r' */
    //VECTOR r_1step_back;  /* r of 1 step  back */
    VECTOR r_orig;
    
    int tor_atom1;
    int tor_atom2;
    int tor_atom3;
    TORS tors;
    int n_rotate1, n_rotate2;
    int *rotate1_lst, *rotate2_lst;
    
    int    n_ngh;
    struct RELAX_STRUCT **ngh;
    int    n_ngh14;
    struct RELAX_STRUCT **ngh14;
    int    n_constr;
    struct RELAX_STRUCT **constr_list;
    double  *constr_dsq;
    
    VECTOR lj_frc;
    VECTOR elec_frc;
    VECTOR torsion_frc;
    VECTOR constr_frc;
    VECTOR frc;
    int    movable;
    int    moving;
    int    moved;
} RELAX;

void relaxation_setup(PROT prot);
void complete_constr(int i_res, int i_conf, int j_res, int j_conf);
int  in_relax_list(int ia);
void setup_nghlst(PROT prot);
void get_frc(float tors_scale, PROT prot);
void get_rp();
int  shake();
void add_conf2relax(CONF *conf_p, int fix);
void add_atom2relax(int ia, int fix);
int res_in_relax(int i_res);
int pick_sidechain(int i_res, PROT prot);
void collect_ngh(int i_res, PROT prot);
int closer_than(int i_res, int i_conf, int j_res, int j_conf, PROT prot, float crg_thr, float dist_thr);

int native_relaxation() {
    CONSTRAINT2 = env.hv_relax_constraint * env.hv_relax_constraint;
    CONSTRAINT_FRC = env.hv_relax_constraint_frc;

    /* setup relaxation */
    printf("   Start setting up for relaxation.\n"); fflush(stdout);
    relaxation_setup(prot);
    id_conf(prot);
    printf("   Setup for relaxation done.\n"); fflush(stdout);
    
    
