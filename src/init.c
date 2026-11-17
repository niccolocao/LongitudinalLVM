#include <stddef.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "likelihood.h"

static const R_CMethodDef methods[] = {
    {"set_cfa_long_options", (DL_FUNC) &set_cfa_long_options, 8, NULL},
    {"set_cfa_long_gamma_time_free", (DL_FUNC) &set_cfa_long_gamma_time_free, 4, NULL},
    {"set_cfa_long_threshold_id", (DL_FUNC) &set_cfa_long_threshold_id, 2, NULL},
    {"set_cfa_long_mean_ident", (DL_FUNC) &set_cfa_long_mean_ident, 2, NULL},
    {"loglik_long", (DL_FUNC) &loglik_long, 23, NULL},
    {"loglik_long_agh", (DL_FUNC) &loglik_long_agh, 24, NULL},
    {"loglik_long_tt_agh", (DL_FUNC) &loglik_long_tt_agh, 27, NULL},
    {"loglik_long_tt_agh_full_equiv", (DL_FUNC) &loglik_long_tt_agh_full_equiv, 28, NULL},
    {"loglik_long_tt", (DL_FUNC) &loglik_long_tt, 26, NULL},
    {"set_growth_tt_options", (DL_FUNC) &set_growth_tt_options, 15, NULL},
    {"loglik_growth_tt", (DL_FUNC) &loglik_growth_tt, 32, NULL},
    {"loglik_growth_tt_agh", (DL_FUNC) &loglik_growth_tt_agh, 34, NULL},
    {"loglik_growth_tt_i_vec", (DL_FUNC) &loglik_growth_tt_i_vec, 32, NULL},
    {"loglik_growth_tt_agh_i_vec", (DL_FUNC) &loglik_growth_tt_agh_i_vec, 34, NULL},
    {"loglik_long_tt_agh_full_equiv_wu", (DL_FUNC) &loglik_long_tt_agh_full_equiv_wu, 27, NULL},
    {"loglik_growth_tt_agh_full_equiv", (DL_FUNC) &loglik_growth_tt_agh_full_equiv, 34, NULL},
    {"loglik_growth_tt_agh_full_equiv_i_vec", (DL_FUNC) &loglik_growth_tt_agh_full_equiv_i_vec, 38, NULL},
    {"loglik_growth_re_tt_agh_nested", (DL_FUNC) &loglik_growth_re_tt_agh_nested, 37, NULL},
    {"loglik_growth_re_tt_agh_full_equiv", (DL_FUNC) &loglik_growth_re_tt_agh_full_equiv, 37, NULL},
    {"loglik_growth_re_tt_agh_full_equiv_score_q1", (DL_FUNC) &loglik_growth_re_tt_agh_full_equiv_score_q1, 39, NULL},
    {"loglik_growth_tt_agh_full_equiv_score_q1", (DL_FUNC) &loglik_growth_tt_agh_full_equiv_score_q1, 36, NULL},
    {"loglik_long_tt_agh_full_equiv_score_q1", (DL_FUNC) &loglik_long_tt_agh_full_equiv_score_q1, 30, NULL},
    {"get_agh_tt_grid_info", (DL_FUNC) &get_agh_tt_grid_info, 5, NULL},
    {"get_agh_tt_grid_copy", (DL_FUNC) &get_agh_tt_grid_copy, 8, NULL},
    {"set_agh_tt_grid_copy", (DL_FUNC) &set_agh_tt_grid_copy, 8, NULL},
    {"get_growth_agh_tt_fe_grid_info", (DL_FUNC) &get_growth_agh_tt_fe_grid_info, 8, NULL},
    {"get_growth_agh_tt_fe_grid_copy", (DL_FUNC) &get_growth_agh_tt_fe_grid_copy, 13, NULL},
    {"set_growth_agh_tt_fe_grid_copy", (DL_FUNC) &set_growth_agh_tt_fe_grid_copy, 13, NULL},
    {"get_growth_re_agh_full_equiv_grid_info", (DL_FUNC) &get_growth_re_agh_full_equiv_grid_info, 7, NULL},
    {"get_growth_re_agh_full_equiv_grid_copy", (DL_FUNC) &get_growth_re_agh_full_equiv_grid_copy, 10, NULL},
    {"set_growth_re_agh_full_equiv_score_grid", (DL_FUNC) &set_growth_re_agh_full_equiv_score_grid, 10, NULL},
    {NULL, NULL, 0, NULL}
};

void attribute_visible R_init_longProbitLVM(DllInfo *dll)
{
    R_registerRoutines(dll, methods, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
