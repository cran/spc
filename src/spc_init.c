#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void ewma_p_arl_be(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ewma_phat_arl_coll(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ewma_phat_crit_coll(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ewma_phat_lambda_coll(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void lns2ewma_arl(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void lns2ewma_crit(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void mewma_ad(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void mewma_arl(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void mewma_arl_f(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void mewma_crit(void *, void *, void *, void *, void *, void *);
extern void mewma_psi(void *, void *, void *, void *, void *, void *, void *);
extern void phat_cdf(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void phat_pdf(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void phat_qf(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void quadrature_nodes_weights(void *, void *, void *, void *, void *);
extern void s_res_ewma_arl(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void scusum_arl(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void scusum_crit(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void scusum_s_arl(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sewma_arl(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sewma_arl_prerun(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sewma_crit(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sewma_crit_prerun(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sewma_q(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sewma_q_crit(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sewma_q_crit_prerun(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sewma_q_prerun(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sewma_sf(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sewma_sf_prerun(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void tol_lim_fac(void *, void *, void *, void *, void *, void *);
extern void tshewhart_ar1_arl(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void x_res_ewma_arl(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xcusum_ad(void *, void *, void *, void *, void *, void *, void *);
extern void xcusum_arl(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xcusum_crit(void *, void *, void *, void *, void *, void *, void *);
extern void xcusum_q(void *, void *, void *, void *, void *, void *, void *, void *);
extern void xcusum_sf(void *, void *, void *, void *, void *, void *, void *, void *);
extern void xDcusum_arl(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xDewma_arl(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xDgrsr_arl(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xewma_ad(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xewma_arl(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xewma_arl_f(void *, void *, void *, void *, void *, void *, void *, void *);
extern void xewma_arl_prerun(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xewma_crit(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xewma_q(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xewma_q_prerun(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xewma_sf(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xewma_sf_prerun(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xgrsr_ad(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xgrsr_arl(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xgrsr_crit(void *, void *, void *, void *, void *, void *, void *, void *);
extern void xsewma_arl(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xsewma_crit(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xsewma_q(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xsewma_q_crit(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xsewma_res_arl(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xsewma_res_pms(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xsewma_sf(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xshewhart_ar1_arl(void *, void *, void *, void *, void *, void *);
extern void xtcusum_arl(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xtewma_ad(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xtewma_arl(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xtewma_q(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void xtewma_sf(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"ewma_p_arl_be",            (DL_FUNC) &ewma_p_arl_be,             9},
    {"ewma_phat_arl_coll",       (DL_FUNC) &ewma_phat_arl_coll,       13},
    {"ewma_phat_crit_coll",      (DL_FUNC) &ewma_phat_crit_coll,      12},
    {"ewma_phat_lambda_coll",    (DL_FUNC) &ewma_phat_lambda_coll,    12},
    {"lns2ewma_arl",             (DL_FUNC) &lns2ewma_arl,              9},
    {"lns2ewma_crit",            (DL_FUNC) &lns2ewma_crit,            11},
    {"mewma_ad",                 (DL_FUNC) &mewma_ad,                 12},
    {"mewma_arl",                (DL_FUNC) &mewma_arl,                10},
    {"mewma_arl_f",              (DL_FUNC) &mewma_arl_f,               9},
    {"mewma_crit",               (DL_FUNC) &mewma_crit,                6},
    {"mewma_psi",                (DL_FUNC) &mewma_psi,                 7},
    {"phat_cdf",                 (DL_FUNC) &phat_cdf,                  9},
    {"phat_pdf",                 (DL_FUNC) &phat_pdf,                  9},
    {"phat_qf",                  (DL_FUNC) &phat_qf,                   9},
    {"quadrature_nodes_weights", (DL_FUNC) &quadrature_nodes_weights,  5},
    {"s_res_ewma_arl",           (DL_FUNC) &s_res_ewma_arl,           11},
    {"scusum_arl",               (DL_FUNC) &scusum_arl,               13},
    {"scusum_crit",              (DL_FUNC) &scusum_crit,              12},
    {"scusum_s_arl",             (DL_FUNC) &scusum_s_arl,             14},
    {"sewma_arl",                (DL_FUNC) &sewma_arl,                11},
    {"sewma_arl_prerun",         (DL_FUNC) &sewma_arl_prerun,         13},
    {"sewma_crit",               (DL_FUNC) &sewma_crit,               14},
    {"sewma_crit_prerun",        (DL_FUNC) &sewma_crit_prerun,        18},
    {"sewma_q",                  (DL_FUNC) &sewma_q,                  11},
    {"sewma_q_crit",             (DL_FUNC) &sewma_q_crit,             16},
    {"sewma_q_crit_prerun",      (DL_FUNC) &sewma_q_crit_prerun,      19},
    {"sewma_q_prerun",           (DL_FUNC) &sewma_q_prerun,           14},
    {"sewma_sf",                 (DL_FUNC) &sewma_sf,                 11},
    {"sewma_sf_prerun",          (DL_FUNC) &sewma_sf_prerun,          14},
    {"tol_lim_fac",              (DL_FUNC) &tol_lim_fac,               6},
    {"tshewhart_ar1_arl",        (DL_FUNC) &tshewhart_ar1_arl,        10},
    {"x_res_ewma_arl",           (DL_FUNC) &x_res_ewma_arl,            9},
    {"xcusum_ad",                (DL_FUNC) &xcusum_ad,                 7},
    {"xcusum_arl",               (DL_FUNC) &xcusum_arl,                9},
    {"xcusum_crit",              (DL_FUNC) &xcusum_crit,               7},
    {"xcusum_q",                 (DL_FUNC) &xcusum_q,                  8},
    {"xcusum_sf",                (DL_FUNC) &xcusum_sf,                 8},
    {"xDcusum_arl",              (DL_FUNC) &xDcusum_arl,              11},
    {"xDewma_arl",               (DL_FUNC) &xDewma_arl,               13},
    {"xDgrsr_arl",               (DL_FUNC) &xDgrsr_arl,               11},
    {"xewma_ad",                 (DL_FUNC) &xewma_ad,                 11},
    {"xewma_arl",                (DL_FUNC) &xewma_arl,                10},
    {"xewma_arl_f",              (DL_FUNC) &xewma_arl_f,               8},
    {"xewma_arl_prerun",         (DL_FUNC) &xewma_arl_prerun,         15},
    {"xewma_crit",               (DL_FUNC) &xewma_crit,               10},
    {"xewma_q",                  (DL_FUNC) &xewma_q,                  11},
    {"xewma_q_prerun",           (DL_FUNC) &xewma_q_prerun,           17},
    {"xewma_sf",                 (DL_FUNC) &xewma_sf,                 11},
    {"xewma_sf_prerun",          (DL_FUNC) &xewma_sf_prerun,          18},
    {"xgrsr_ad",                 (DL_FUNC) &xgrsr_ad,                  9},
    {"xgrsr_arl",                (DL_FUNC) &xgrsr_arl,                10},
    {"xgrsr_crit",               (DL_FUNC) &xgrsr_crit,                8},
    {"xsewma_arl",               (DL_FUNC) &xsewma_arl,               16},
    {"xsewma_crit",              (DL_FUNC) &xsewma_crit,              15},
    {"xsewma_q",                 (DL_FUNC) &xsewma_q,                 16},
    {"xsewma_q_crit",            (DL_FUNC) &xsewma_q_crit,            18},
    {"xsewma_res_arl",           (DL_FUNC) &xsewma_res_arl,           15},
    {"xsewma_res_pms",           (DL_FUNC) &xsewma_res_pms,           16},
    {"xsewma_sf",                (DL_FUNC) &xsewma_sf,                16},
    {"xshewhart_ar1_arl",        (DL_FUNC) &xshewhart_ar1_arl,         6},
    {"xtcusum_arl",              (DL_FUNC) &xtcusum_arl,               9},
    {"xtewma_ad",                (DL_FUNC) &xtewma_ad,                13},
    {"xtewma_arl",               (DL_FUNC) &xtewma_arl,               12},
    {"xtewma_q",                 (DL_FUNC) &xtewma_q,                 13},
    {"xtewma_sf",                (DL_FUNC) &xtewma_sf,                13},
    {NULL, NULL, 0}
};

void R_init_spc(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
