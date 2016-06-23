#ifndef _Q_CLOVER_PHI_H
#define _Q_CLOVER_PHI_H

namespace cvc {

void clover_term_init (double***s, int nmat);
void clover_term_fini (double***s);

void clover_term_eo (double**s, double*gauge_field);
void Q_clover_phi_eo (double *e_new, double *o_new, double *e_old, double *o_old, double *gauge_field, double mass, double *aux, double **cl);
void M_clover_zz (double*s, double*r, double mass, double*cl);

void clover_mzz_matrix (double**mzz, double**cl, double mu, double csw);
void clover_mzz_inv_matrix (double**mzzinv, double**mzz);
void M_clover_zz_matrix (double*s, double*r, double*mzz);
void M_clover_zz_inv_matrix (double*s, double*r, double *mzzinv);
void Q_clover_phi_matrix_eo (double *e_new, double *o_new, double *e_old, double *o_old, double *gauge_field, double *aux, double**mzz);
void C_clover_oo (double*s, double*r, double *gauge_field, double *s_aux, double*mzz, double*mzzinv);
void X_clover_eo (double *even, double *odd, double *gauge_field, double*mzzinv);
void C_clover_from_Xeo (double *t, double *s, double *r, double *gauge_field, double*mzz);

}  /* end of namespace cvc */
#endif