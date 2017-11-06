#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void check_nthread(void *);
extern void check_os(void *);
extern void dgm_bdmcmc_mpl_binary_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_bdmcmc_mpl_binary_ma_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_bdmcmc_mpl_binary_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_bdmcmc_mpl_binary_map_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_bdmcmc_mpl_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_bdmcmc_mpl_ma_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_bdmcmc_mpl_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_bdmcmc_mpl_map_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_rjmcmc_mpl_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dgm_rjmcmc_mpl_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_bdmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_bdmcmc_ma_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_bdmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_bdmcmc_map_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_DMH_bdmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_DMH_bdmcmc_ma_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_DMH_bdmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_DMH_bdmcmc_map_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_DMH_rjmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_DMH_rjmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_rjmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcgm_rjmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_ma_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_map_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_mpl_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_mpl_ma_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_mpl_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_bdmcmc_mpl_map_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_DMH_bdmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_DMH_bdmcmc_ma_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_DMH_bdmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_DMH_bdmcmc_map_multi_update(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_DMH_rjmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_DMH_rjmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_rjmcmc_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_rjmcmc_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_rjmcmc_mpl_ma(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggm_rjmcmc_mpl_map(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void log_exp_mc(void *, void *, void *, void *, void *, void *, void *, void *);
extern void log_mpl_binary_parallel_hc(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void log_mpl_dis(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void omp_set_num_cores(void *);
extern void rgwish_c(void *, void *, void *, void *, void *);
extern void rwish_c(void *, void *, void *, void *);
extern void scale_free(void *, void *);
extern void transfer_data(void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"check_nthread",                          (DL_FUNC) &check_nthread,                           1},
    {"check_os",                               (DL_FUNC) &check_os,                                1},
    {"dgm_bdmcmc_mpl_binary_ma",               (DL_FUNC) &dgm_bdmcmc_mpl_binary_ma,               13},
    {"dgm_bdmcmc_mpl_binary_ma_multi_update",  (DL_FUNC) &dgm_bdmcmc_mpl_binary_ma_multi_update,  14},
    {"dgm_bdmcmc_mpl_binary_map",              (DL_FUNC) &dgm_bdmcmc_mpl_binary_map,              17},
    {"dgm_bdmcmc_mpl_binary_map_multi_update", (DL_FUNC) &dgm_bdmcmc_mpl_binary_map_multi_update, 19},
    {"dgm_bdmcmc_mpl_ma",                      (DL_FUNC) &dgm_bdmcmc_mpl_ma,                      13},
    {"dgm_bdmcmc_mpl_ma_multi_update",         (DL_FUNC) &dgm_bdmcmc_mpl_ma_multi_update,         14},
    {"dgm_bdmcmc_mpl_map",                     (DL_FUNC) &dgm_bdmcmc_mpl_map,                     17},
    {"dgm_bdmcmc_mpl_map_multi_update",        (DL_FUNC) &dgm_bdmcmc_mpl_map_multi_update,        19},
    {"dgm_rjmcmc_mpl_ma",                      (DL_FUNC) &dgm_rjmcmc_mpl_ma,                      13},
    {"dgm_rjmcmc_mpl_map",                     (DL_FUNC) &dgm_rjmcmc_mpl_map,                     17},
    {"gcgm_bdmcmc_ma",                         (DL_FUNC) &gcgm_bdmcmc_ma,                         19},
    {"gcgm_bdmcmc_ma_multi_update",            (DL_FUNC) &gcgm_bdmcmc_ma_multi_update,            20},
    {"gcgm_bdmcmc_map",                        (DL_FUNC) &gcgm_bdmcmc_map,                        23},
    {"gcgm_bdmcmc_map_multi_update",           (DL_FUNC) &gcgm_bdmcmc_map_multi_update,           25},
    {"gcgm_DMH_bdmcmc_ma",                     (DL_FUNC) &gcgm_DMH_bdmcmc_ma,                     20},
    {"gcgm_DMH_bdmcmc_ma_multi_update",        (DL_FUNC) &gcgm_DMH_bdmcmc_ma_multi_update,        21},
    {"gcgm_DMH_bdmcmc_map",                    (DL_FUNC) &gcgm_DMH_bdmcmc_map,                    24},
    {"gcgm_DMH_bdmcmc_map_multi_update",       (DL_FUNC) &gcgm_DMH_bdmcmc_map_multi_update,       26},
    {"gcgm_DMH_rjmcmc_ma",                     (DL_FUNC) &gcgm_DMH_rjmcmc_ma,                     20},
    {"gcgm_DMH_rjmcmc_map",                    (DL_FUNC) &gcgm_DMH_rjmcmc_map,                    24},
    {"gcgm_rjmcmc_ma",                         (DL_FUNC) &gcgm_rjmcmc_ma,                         19},
    {"gcgm_rjmcmc_map",                        (DL_FUNC) &gcgm_rjmcmc_map,                        23},
    {"ggm_bdmcmc_ma",                          (DL_FUNC) &ggm_bdmcmc_ma,                          14},
    {"ggm_bdmcmc_ma_multi_update",             (DL_FUNC) &ggm_bdmcmc_ma_multi_update,             15},
    {"ggm_bdmcmc_map",                         (DL_FUNC) &ggm_bdmcmc_map,                         18},
    {"ggm_bdmcmc_map_multi_update",            (DL_FUNC) &ggm_bdmcmc_map_multi_update,            20},
    {"ggm_bdmcmc_mpl_ma",                      (DL_FUNC) &ggm_bdmcmc_mpl_ma,                       9},
    {"ggm_bdmcmc_mpl_ma_multi_update",         (DL_FUNC) &ggm_bdmcmc_mpl_ma_multi_update,         10},
    {"ggm_bdmcmc_mpl_map",                     (DL_FUNC) &ggm_bdmcmc_mpl_map,                     13},
    {"ggm_bdmcmc_mpl_map_multi_update",        (DL_FUNC) &ggm_bdmcmc_mpl_map_multi_update,        15},
    {"ggm_DMH_bdmcmc_ma",                      (DL_FUNC) &ggm_DMH_bdmcmc_ma,                      16},
    {"ggm_DMH_bdmcmc_ma_multi_update",         (DL_FUNC) &ggm_DMH_bdmcmc_ma_multi_update,         17},
    {"ggm_DMH_bdmcmc_map",                     (DL_FUNC) &ggm_DMH_bdmcmc_map,                     20},
    {"ggm_DMH_bdmcmc_map_multi_update",        (DL_FUNC) &ggm_DMH_bdmcmc_map_multi_update,        22},
    {"ggm_DMH_rjmcmc_ma",                      (DL_FUNC) &ggm_DMH_rjmcmc_ma,                      16},
    {"ggm_DMH_rjmcmc_map",                     (DL_FUNC) &ggm_DMH_rjmcmc_map,                     20},
    {"ggm_rjmcmc_ma",                          (DL_FUNC) &ggm_rjmcmc_ma,                          14},
    {"ggm_rjmcmc_map",                         (DL_FUNC) &ggm_rjmcmc_map,                         18},
    {"ggm_rjmcmc_mpl_ma",                      (DL_FUNC) &ggm_rjmcmc_mpl_ma,                       9},
    {"ggm_rjmcmc_mpl_map",                     (DL_FUNC) &ggm_rjmcmc_mpl_map,                     13},
    {"log_exp_mc",                             (DL_FUNC) &log_exp_mc,                              8},
    {"log_mpl_binary_parallel_hc",             (DL_FUNC) &log_mpl_binary_parallel_hc,              9},
    {"log_mpl_dis",                            (DL_FUNC) &log_mpl_dis,                            11},
    {"omp_set_num_cores",                      (DL_FUNC) &omp_set_num_cores,                       1},
    {"rgwish_c",                               (DL_FUNC) &rgwish_c,                                5},
    {"rwish_c",                                (DL_FUNC) &rwish_c,                                 4},
    {"scale_free",                             (DL_FUNC) &scale_free,                              2},
    {"transfer_data",                          (DL_FUNC) &transfer_data,                           5},
    {NULL, NULL, 0}
};

void R_init_BDgraph(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
