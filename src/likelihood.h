#ifndef LONGPROBITLVM_LIKELIHOOD_H
#define LONGPROBITLVM_LIKELIHOOD_H

void set_cfa_long_options(
    int *, int *, int *, int *, int *, int *, double *, int *);
void set_cfa_long_gamma_time_free(
    int *, int *, int *, int *);
void set_cfa_long_threshold_id(
    int *, int *);
void set_cfa_long_mean_ident(
    int *, int *);
void loglik_long(
    double *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, int *, int *, double *, double *, int *, int *, int *, int *,
    int *, int *, int *, int *, double *);
void loglik_long_agh(
    double *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, int *, int *, double *, double *, int *, int *, int *, int *,
    int *, int *, int *, int *, int *, double *);
void loglik_long_tt_agh(
    double *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, int *, int *, double *, double *, int *, double *, double *, int *,
    int *, int *, int *, int *, int *, int *, int *, int *, double *);
void loglik_long_tt_agh_full_equiv(
    double *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, int *, int *, double *, double *, int *, double *, double *, int *,
    int *, int *, int *, int *, int *, int *, int *, int *, int *,
    double *);
void loglik_long_tt(
    double *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, int *, double *, double *, int *, double *, double *, int *, int *,
    int *, int *, int *, int *, int *, int *, int *, double *);
void set_growth_tt_options(
    int *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, int *, int *, double *, int *, int *);
void loglik_growth_tt(
    double *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, int *, double *, double *, int *, double *, double *, int *, int *,
    int *, int *, int *, int *, int *, int *, int *, int *, double *,
    double *, double *, int *, int *, double *);
void loglik_growth_tt_agh(
    double *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, int *, double *, double *, int *, double *, double *, int *, int *,
    int *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, double *, double *, double *, int *, int *, double *);
void loglik_growth_tt_i_vec(
    double *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, int *, double *, double *, int *, double *, double *, int *, int *,
    int *, int *, int *, int *, int *, int *, int *, int *, double *,
    double *, double *, int *, int *, double *);
void loglik_growth_tt_agh_i_vec(
    double *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, int *, double *, double *, int *, double *, double *, int *, int *,
    int *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, double *, double *, double *, int *, int *, double *);
void loglik_long_tt_agh_full_equiv_wu(
    double *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, int *, int *, double *, double *, int *, double *, double *, int *,
    int *, int *, int *, int *, int *, int *, int *, int *, double *);
void loglik_growth_tt_agh_full_equiv(
    double *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, int *, double *, double *, int *, double *, double *, int *, int *,
    int *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, double *, double *, double *, int *, int *, double *);
void loglik_growth_tt_agh_full_equiv_i_vec(
    double *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, int *, double *, double *, int *, double *, double *, int *, int *,
    int *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, double *, double *, double *, int *, int *, double *);
void loglik_growth_re_tt_agh_nested(
    double *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, int *, double *, double *, int *, double *, double *, int *, double *,
    double *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, int *, int *, int *, double *, double *, double *, int *, int *,
    double *);
void loglik_growth_re_tt_agh_full_equiv(
    double *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, int *, double *, double *, int *, double *, double *, int *, double *,
    double *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, int *, int *, int *, double *, double *, double *, int *, int *,
    double *);
void loglik_growth_re_tt_agh_full_equiv_score_q1(
    double *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, int *, double *, double *, int *, double *, double *, int *, double *,
    double *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, int *, int *, int *, double *, double *, double *, int *, int *,
    double *, double *, int *);
void loglik_growth_tt_agh_full_equiv_score_q1(
    double *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, int *, double *, double *, int *, double *, double *, int *, int *,
    int *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, double *, double *, double *, int *, int *, double *, double *, int *);
void loglik_long_tt_agh_full_equiv_score_q1(
    double *, int *, int *, int *, int *, int *, int *, int *, int *,
    int *, int *, int *, double *, double *, int *, double *, double *, int *,
    int *, int *, int *, int *, int *, int *, int *, int *, int *,
    double *, double *, int *);
void get_agh_tt_grid_info(
    int *, int *, int *, int *, int *);
void get_agh_tt_grid_copy(
    double *, double *, double *, int *, int *, int *, int *, int *);
void set_agh_tt_grid_copy(
    double *, double *, double *, int *, int *, int *, int *, int *);
void get_growth_agh_tt_fe_grid_info(
    int *, int *, int *, int *, int *, int *, int *, int *);
void get_growth_agh_tt_fe_grid_copy(
    double *, double *, double *, double *, double *, int *, int *, int *, int *,
    int *, int *, int *, int *);
void set_growth_agh_tt_fe_grid_copy(
    double *, double *, double *, double *, double *, int *, int *, int *, int *,
    int *, int *, int *, int *);

void get_growth_re_agh_full_equiv_grid_info(
    int *ready0,
    int *fit_id0,
    int *dim_w0,
    int *dim_z0,
    int *nsub0,
    int *nquad_u0,
    int *nquad_zeta0
);

void get_growth_re_agh_full_equiv_grid_copy(
    double *mode_out,
    double *B_out,
    double *logdet_w_out,
    int *ready0,
    int *fit_id0,
    int *dim_w0,
    int *dim_z0,
    int *nsub0,
    int *nquad_u0,
    int *nquad_zeta0
);

void set_growth_re_agh_full_equiv_score_grid(
    double *mode_in,
    double *B_in,
    double *logdet_w_in,
    int *dim_w0,
    int *dim_z0,
    int *nsub0,
    int *nquad_u0,
    int *nquad_zeta0,
    int *fit_id0,
    int *ready0
);

#endif
