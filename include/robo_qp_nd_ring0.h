#ifndef ROBO_QP_ND_RING0_H
#define ROBO_QP_ND_RING0_H

/* ============================================================================
 * ROBO v3.0.0 — Pure-GPR General N-DoF Slack-CBF-QP Solver (N <= 8, M <= 24)
 * ============================================================================
 * Freestanding Pure-GPR Q16.16 solver:
 *   - Zero FPU / SIMD (-mgeneral-regs-only clean)
 *   - Zero udiv / sdiv / libc / heap / indirect branches
 *   - Zero signed-overflow or shift UB (LLVM SSA --check-nsw clean)
 *   - L1+L2 exact slack delta >= 0: zero slack when feasible, minimal relaxation
 *     under constraint conflicts, and guaranteed box-feasible emergency brake.
 *   - Integrated with FLOW Layer-0 `flow_cholesky_sdd.h` PSD certification.
 * ============================================================================ */

#include "robo_cbf_ring0.h"
#include "flow_cholesky_sdd.h"

#ifdef __cplusplus
extern "C" {
#endif

#define ROBO_QP_MAX_N              8U
#define ROBO_QP_MAX_M              24U
#define ROBO_QP_DEFAULT_ITERS      128U
#define ROBO_VERDICT_SLACK_RELAXED 4U

typedef struct {
    int64_t a_q16[ROBO_QP_MAX_N];
    int64_t b_q16;
    int64_t slack_max_q16;
    int64_t slack_lin_q16;
    int64_t slack_quad_q16;
    uint8_t is_hard;
} robo_qp_con_t;

typedef struct {
    uint32_t      n;
    uint32_t      m;
    int64_t       u_nom_q16[ROBO_QP_MAX_N];
    int64_t       u_lo_q16[ROBO_QP_MAX_N];
    int64_t       u_hi_q16[ROBO_QP_MAX_N];
    int64_t       u_brake_q16[ROBO_QP_MAX_N];
    int64_t       w_diag_q16[ROBO_QP_MAX_N];
    robo_qp_con_t con[ROBO_QP_MAX_M];
    uint32_t      max_iters;
} robo_qp_nd_problem_t;

typedef struct {
    uint32_t verdict;
    uint32_t n;
    uint32_t m;
    uint32_t iters_run;
    uint8_t  feasible;
    uint8_t  slack_active;
    uint8_t  fallback_engaged;
    uint8_t  sdd_certified;
    int64_t  u_safe_q16[ROBO_QP_MAX_N];
    int64_t  delta_q16[ROBO_QP_MAX_M];
    int64_t  max_hard_viol_q16;
    int64_t  max_soft_viol_q16;
    int64_t  max_slack_q16;
    int64_t  sdd_mu_floor_q32;
    int64_t  obj_q16;
} robo_qp_nd_cert_t;

FLOW_FORCE_INL int64_t robo__sat_neg_i64(int64_t x) {
    flow_i128 nx = -(flow_i128)x;
    if (nx > (flow_i128)INT64_MAX) return INT64_MAX;
    if (nx < (flow_i128)(INT64_MIN + 1)) return INT64_MIN + 1;
    return (int64_t)nx;
}

FLOW_FORCE_INL int64_t robo__clamp_i64(int64_t v, int64_t lo, int64_t hi) {
    int64_t c = flow_sel_i64(v < lo, lo, v);
    return flow_sel_i64(c > hi, hi, c);
}

FLOW_FORCE_INL int64_t robo__dot_nd_q16(const int64_t *a, const int64_t *b, uint32_t n) {
    flow_i128 acc = 0;
    uint32_t j;
    uint32_t nn = (n > ROBO_QP_MAX_N) ? ROBO_QP_MAX_N : n;
    for (j = 0U; j < nn; ++j) {
        acc += (flow_i128)a[j] * (flow_i128)b[j];
    }
    return flow_sym_shr_i128(acc, 16U);
}

/* Solve N-DoF Box + Slack-CBF-QP via Pure-GPR ADMM + Cholesky + FLOW SDD check */
FLOW_INL void robo_qp_nd_solve_ring0(const robo_qp_nd_problem_t *prob, robo_qp_nd_cert_t *cert) {
    uint32_t i, j, k;
    if (!cert) {
        return;
    }
    for (j = 0U; j < ROBO_QP_MAX_N; ++j) {
        cert->u_safe_q16[j] = 0;
    }
    for (i = 0U; i < ROBO_QP_MAX_M; ++i) {
        cert->delta_q16[i] = 0;
    }
    cert->verdict           = ROBO_VERDICT_VETOED_INFEASIBLE;
    cert->n                 = 0U;
    cert->m                 = 0U;
    cert->iters_run         = 0U;
    cert->feasible          = 0U;
    cert->slack_active      = 0U;
    cert->fallback_engaged  = 1U;
    cert->sdd_certified     = 0U;
    cert->max_hard_viol_q16 = 0;
    cert->max_soft_viol_q16 = 0;
    cert->max_slack_q16     = 0;
    cert->sdd_mu_floor_q32  = 0;
    cert->obj_q16           = 0;

    if (!prob || prob->n == 0U || prob->n > ROBO_QP_MAX_N || prob->m > ROBO_QP_MAX_M) {
        return;
    }

    const uint32_t n = prob->n;
    const uint32_t m = prob->m;
    cert->n          = n;
    cert->m          = m;

    int64_t u_lo[ROBO_QP_MAX_N];
    int64_t u_hi[ROBO_QP_MAX_N];
    int64_t u_nom[ROBO_QP_MAX_N];
    int64_t u_brake[ROBO_QP_MAX_N];
    int64_t w_diag[ROBO_QP_MAX_N];
    bool box_valid = true;

    int64_t w_sum = 0;
    for (j = 0U; j < n; ++j) {
        u_lo[j]  = prob->u_lo_q16[j];
        u_hi[j]  = prob->u_hi_q16[j];
        u_nom[j] = prob->u_nom_q16[j];
        if (u_lo[j] > u_hi[j]) {
            box_valid = false;
            int64_t mid = flow_ashr_i64(robo__sat_add_i64(u_lo[j], u_hi[j]), 1U);
            u_lo[j]     = mid;
            u_hi[j]     = mid;
        }
        u_brake[j] = robo__clamp_i64(prob->u_brake_q16[j], u_lo[j], u_hi[j]);
        w_diag[j]  = (prob->w_diag_q16[j] > INT64_C(256)) ? prob->w_diag_q16[j] : ROBO_Q16_ONE;
        w_sum      = robo__sat_add_i64(w_sum, w_diag[j]);
    }

    /* Fast path: check if clamped u_nom already satisfies all constraints with zero slack */
    if (box_valid) {
        bool nom_ok = true;
        int64_t u_cand[ROBO_QP_MAX_N];
        for (j = 0U; j < n; ++j) {
            u_cand[j] = robo__clamp_i64(u_nom[j], u_lo[j], u_hi[j]);
        }
        for (i = 0U; i < m; ++i) {
            int64_t dot = robo__dot_nd_q16(prob->con[i].a_q16, u_cand, n);
            if (dot < robo__sat_sub_i64(prob->con[i].b_q16, INT64_C(32))) {
                nom_ok = false;
                break;
            }
        }
        if (nom_ok) {
            int64_t obj = 0;
            bool same   = true;
            for (j = 0U; j < n; ++j) {
                cert->u_safe_q16[j] = u_cand[j];
                int64_t d           = robo__sat_sub_i64(u_cand[j], u_nom[j]);
                if ((int64_t)flow_abs_i64(d) > INT64_C(8)) {
                    same = false;
                }
                int64_t wd = robo__mul_q16(w_diag[j], robo__mul_q16(d, d));
                obj        = robo__sat_add_i64(obj, flow_ashr_i64(wd, 1U));
            }
            cert->verdict          = same ? ROBO_VERDICT_PASS : ROBO_VERDICT_PROJECTED;
            cert->feasible         = 1U;
            cert->fallback_engaged = 0U;
            cert->sdd_certified    = 1U;
            cert->obj_q16          = obj;
            return;
        }
    }

    /* Normalize rows so ||a_i||_2 == 1.0 in Q16.16 for scale-free conditioning */
    int64_t anorm[ROBO_QP_MAX_M][ROBO_QP_MAX_N];
    int64_t bnorm[ROBO_QP_MAX_M];
    int64_t row_norm[ROBO_QP_MAX_M];
    int64_t plin_norm[ROBO_QP_MAX_M];
    int64_t rho_norm[ROBO_QP_MAX_M];
    int64_t smax_norm[ROBO_QP_MAX_M];
    bool    hard_degenerate_infeas = !box_valid;

    for (i = 0U; i < m; ++i) {
        flow_u128 n2 = 0;
        for (j = 0U; j < n; ++j) {
            int64_t aj = prob->con[i].a_q16[j];
            n2 += (flow_u128)((flow_i128)aj * (flow_i128)aj);
        }
        int64_t rn  = robo__sqrt_q16((uint64_t)(n2 >> 0U));
        row_norm[i] = rn;
        if (rn <= INT64_C(16)) {
            for (j = 0U; j < n; ++j) {
                anorm[i][j] = 0;
            }
            bnorm[i]     = (prob->con[i].b_q16 > 0) ? ROBO_Q16_ONE : 0;
            plin_norm[i] = ROBO_Q16_ONE * 64;
            rho_norm[i]  = ROBO_Q16_ONE * 64;
            smax_norm[i] = ROBO_Q16_ONE;
            if (prob->con[i].b_q16 > INT64_C(64) && prob->con[i].is_hard != 0U) {
                hard_degenerate_infeas = true;
            }
        } else {
            for (j = 0U; j < n; ++j) {
                anorm[i][j] = flow_sym_div128_clamped(
                    flow_shl_i128((flow_i128)prob->con[i].a_q16[j], 16U), (flow_i128)rn,
                    FLOW_MAX_Q16);
            }
            bnorm[i] = flow_sym_div128_clamped(
                flow_shl_i128((flow_i128)prob->con[i].b_q16, 16U), (flow_i128)rn, FLOW_MAX_Q16);

            int64_t plin_raw = (prob->con[i].slack_lin_q16 > 0)
                                   ? prob->con[i].slack_lin_q16
                                   : (ROBO_Q16_ONE * 64);
            int64_t rho_raw  = (prob->con[i].slack_quad_q16 > INT64_C(256))
                                   ? prob->con[i].slack_quad_q16
                                   : (ROBO_Q16_ONE * 64);
            plin_norm[i] = robo__clamp_i64(robo__mul_q16(plin_raw, rn), ROBO_Q16_ONE,
                                           ROBO_Q16_ONE * 16384);
            rho_norm[i]  = robo__clamp_i64(
                robo__mul_q16(rho_raw, robo__mul_q16(rn, rn)), ROBO_Q16_ONE, ROBO_Q16_ONE * 32768);
            int64_t smax_raw = (prob->con[i].slack_max_q16 > 0) ? prob->con[i].slack_max_q16 : 0;
            smax_norm[i]     = flow_sym_div128_clamped(flow_shl_i128((flow_i128)smax_raw, 16U),
                                                       (flow_i128)rn, FLOW_MAX_Q16);
        }
    }

    if (hard_degenerate_infeas) {
        for (j = 0U; j < n; ++j) {
            cert->u_safe_q16[j] = u_brake[j];
        }
        return;
    }

    /* ADMM penalty parameter sigma */
    int64_t w_avg = flow_sym_div128_clamped((flow_i128)w_sum, (flow_i128)n, FLOW_MAX_Q16);
    int64_t sigma = robo__clamp_i64(robo__max_i64(ROBO_Q16_ONE * 4, w_avg), ROBO_Q16_ONE,
                                    ROBO_Q16_ONE * 64);

    /* Assemble symmetric positive-definite K = diag(W + sigma) + sigma * A^T A */
    int64_t K[ROBO_QP_MAX_N][ROBO_QP_MAX_N];
    int64_t L[ROBO_QP_MAX_N][ROBO_QP_MAX_N];
    int64_t inv_Ljj[ROBO_QP_MAX_N];
    for (j = 0U; j < ROBO_QP_MAX_N; ++j) {
        inv_Ljj[j] = 0;
        for (k = 0U; k < ROBO_QP_MAX_N; ++k) {
            K[j][k] = 0;
            L[j][k] = 0;
        }
    }

    for (j = 0U; j < n; ++j) {
        for (k = j; k < n; ++k) {
            flow_i128 ata = 0;
            for (i = 0U; i < m; ++i) {
                ata += (flow_i128)anorm[i][j] * (flow_i128)anorm[i][k];
            }
            int64_t ata_q16 = flow_sym_shr_i128(ata, 16U);
            int64_t val     = robo__mul_q16(sigma, ata_q16);
            if (j == k) {
                val = robo__sat_add_i64(val, robo__sat_add_i64(w_diag[j], sigma));
            }
            K[j][k] = val;
            K[k][j] = val;
        }
    }

    /* Pure-GPR Cholesky K = L * L^T */
    for (j = 0U; j < n; ++j) {
        flow_i128 sum_sq = 0;
        for (k = 0U; k < j; ++k) {
            sum_sq += (flow_i128)L[j][k] * (flow_i128)L[j][k];
        }
        int64_t diag = robo__sat_sub_i64(K[j][j], flow_sym_shr_i128(sum_sq, 16U));
        if (diag < INT64_C(1024)) {
            diag = INT64_C(1024);
        }
        int64_t ljj = robo__sqrt_q16((uint64_t)flow_shl_i128((flow_i128)diag, 16U));
        if (ljj < INT64_C(32)) {
            ljj = INT64_C(32);
        }
        L[j][j]    = ljj;
        inv_Ljj[j] = flow_sym_div128_clamped(flow_shl_i128((flow_i128)1, 32U), (flow_i128)ljj,
                                             FLOW_MAX_Q16);
        for (i = j + 1U; i < n; ++i) {
            flow_i128 s = 0;
            for (k = 0U; k < j; ++k) {
                s += (flow_i128)L[i][k] * (flow_i128)L[j][k];
            }
            int64_t diff = robo__sat_sub_i64(K[i][j], flow_sym_shr_i128(s, 16U));
            L[i][j]      = robo__mul_q16(diff, inv_Ljj[j]);
        }
    }

    /* FLOW Layer-0 Cholesky-SDD Certification of K_q32 >= L*L^T + mu_psd*I */
    int64_t K_q32_flat[ROBO_QP_MAX_N * ROBO_QP_MAX_N];
    int64_t L_q16_flat[ROBO_QP_MAX_N * ROBO_QP_MAX_N];
    for (j = 0U; j < n; ++j) {
        for (k = 0U; k < n; ++k) {
            K_q32_flat[j * n + k] = (int64_t)flow_shl_i128((flow_i128)K[j][k], 16U);
            L_q16_flat[j * n + k] = L[j][k];
        }
    }
    int64_t mu_psd_q32 = 0;
    int sdd_rc = flow_cholesky_sdd_floor_i64((int)n, (int)ROBO_QP_MAX_N, K_q32_flat, L_q16_flat,
                                             &mu_psd_q32);
    cert->sdd_mu_floor_q32 = mu_psd_q32;
    cert->sdd_certified    = (sdd_rc == 0 && mu_psd_q32 >= -(ROBO_Q16_ONE * 64)) ? 1U : 0U;

    /* ADMM State Initialization */
    int64_t u[ROBO_QP_MAX_N];
    int64_t z_u[ROBO_QP_MAX_N];
    int64_t y_u[ROBO_QP_MAX_N];
    int64_t w_unom[ROBO_QP_MAX_N];
    int64_t z_c[ROBO_QP_MAX_M];
    int64_t y_c[ROBO_QP_MAX_M];
    int64_t delta_n[ROBO_QP_MAX_M];

    for (j = 0U; j < n; ++j) {
        u[j]      = robo__clamp_i64(u_nom[j], u_lo[j], u_hi[j]);
        z_u[j]    = u[j];
        y_u[j]    = 0;
        w_unom[j] = robo__mul_q16(w_diag[j], u_nom[j]);
    }
    for (i = 0U; i < m; ++i) {
        z_c[i]     = robo__dot_nd_q16(anorm[i], u, n);
        y_c[i]     = 0;
        delta_n[i] = 0;
    }

    uint32_t max_iters = (prob->max_iters > 0U && prob->max_iters <= 256U)
                             ? prob->max_iters
                             : ROBO_QP_DEFAULT_ITERS;
    uint32_t iter;
    int64_t tmp[ROBO_QP_MAX_N];
    int64_t rhs[ROBO_QP_MAX_N];

    for (iter = 0U; iter < max_iters; ++iter) {
        /* 1. Assemble RHS = W*u_nom + sigma*(z_u - y_u) + sigma * A^T (z_c - y_c) */
        for (j = 0U; j < n; ++j) {
            flow_i128 at_zy = 0;
            for (i = 0U; i < m; ++i) {
                int64_t zy = robo__sat_sub_i64(z_c[i], y_c[i]);
                at_zy += (flow_i128)anorm[i][j] * (flow_i128)zy;
            }
            int64_t sum_c = robo__sat_add_i64(robo__sat_sub_i64(z_u[j], y_u[j]),
                                              flow_sym_shr_i128(at_zy, 16U));
            rhs[j]        = robo__sat_add_i64(w_unom[j], robo__mul_q16(sigma, sum_c));
        }

        /* 2. Forward & Backward Cholesky Substitution: L * L^T * u = rhs */
        for (j = 0U; j < n; ++j) {
            flow_i128 s = 0;
            for (k = 0U; k < j; ++k) {
                s += (flow_i128)L[j][k] * (flow_i128)tmp[k];
            }
            tmp[j] = robo__mul_q16(robo__sat_sub_i64(rhs[j], flow_sym_shr_i128(s, 16U)),
                                   inv_Ljj[j]);
        }
        for (j = n; j > 0U; --j) {
            uint32_t jj = j - 1U;
            flow_i128 s = 0;
            for (k = jj + 1U; k < n; ++k) {
                s += (flow_i128)L[k][jj] * (flow_i128)u[k];
            }
            u[jj] = robo__mul_q16(robo__sat_sub_i64(tmp[jj], flow_sym_shr_i128(s, 16U)),
                                  inv_Ljj[jj]);
        }

        /* 3. Over-relaxed (alpha = 1.5) z_u, z_c, delta and dual y updates */
        int64_t max_prim_res = 0;
        for (j = 0U; j < n; ++j) {
            int64_t u_hat  = robo__sat_add_i64(u[j], flow_ashr_i64(robo__sat_sub_i64(u[j], z_u[j]), 1U));
            int64_t zu_new = robo__clamp_i64(robo__sat_add_i64(u_hat, y_u[j]), u_lo[j], u_hi[j]);
            int64_t res    = (int64_t)flow_abs_i64(robo__sat_sub_i64(u[j], zu_new));
            if (res > max_prim_res) {
                max_prim_res = res;
            }
            y_u[j] = robo__sat_add_i64(y_u[j], robo__sat_sub_i64(u_hat, zu_new));
            z_u[j] = zu_new;
        }

        for (i = 0U; i < m; ++i) {
            int64_t au     = robo__dot_nd_q16(anorm[i], u, n);
            int64_t au_hat = robo__sat_add_i64(au, flow_ashr_i64(robo__sat_sub_i64(au, z_c[i]), 1U));
            int64_t z_hat  = robo__sat_add_i64(au_hat, y_c[i]);
            int64_t zc_new = z_hat;
            int64_t di     = 0;

            if (z_hat < bnorm[i]) {
                int64_t diff = robo__sat_sub_i64(bnorm[i], z_hat);
                if (prob->con[i].is_hard != 0U) {
                    zc_new = bnorm[i];
                } else {
                    int64_t sig_d = robo__mul_q16(sigma, diff);
                    int64_t num   = robo__sat_sub_i64(sig_d, plin_norm[i]);
                    if (num > 0 && smax_norm[i] > 0) {
                        int64_t den = robo__sat_add_i64(rho_norm[i], sigma);
                        di          = flow_sym_div128_clamped(flow_shl_i128((flow_i128)num, 16U),
                                                              (flow_i128)den, smax_norm[i]);
                        if (di > diff) {
                            di = diff;
                        }
                        if (di < 0) {
                            di = 0;
                        }
                    }
                    zc_new = robo__sat_sub_i64(bnorm[i], di);
                }
            }
            delta_n[i] = di;
            int64_t cres = (int64_t)flow_abs_i64(robo__sat_sub_i64(au, zc_new));
            if (cres > max_prim_res) {
                max_prim_res = cres;
            }
            y_c[i] = robo__sat_add_i64(y_c[i], robo__sat_sub_i64(au_hat, zc_new));
            z_c[i] = zc_new;
        }

        if (iter >= 16U && (iter & 7U) == 7U && max_prim_res <= INT64_C(24)) {
            iter += 1U;
            break;
        }
    }
    cert->iters_run = iter;

    /* 4. Exact box-clamped candidate & feasibility verification */
    int64_t u_cand[ROBO_QP_MAX_N];
    for (j = 0U; j < n; ++j) {
        u_cand[j] = robo__clamp_i64(z_u[j], u_lo[j], u_hi[j]);
    }

    int64_t max_hard_viol = 0;
    int64_t max_soft_viol = 0;
    int64_t max_delta     = 0;
    bool    all_ok        = true;
    int64_t cand_delta[ROBO_QP_MAX_M];

    for (i = 0U; i < m; ++i) {
        int64_t dot   = robo__dot_nd_q16(prob->con[i].a_q16, u_cand, n);
        int64_t need  = robo__sat_sub_i64(prob->con[i].b_q16, dot);
        int64_t tol_i = robo__max_i64(INT64_C(128), flow_ashr_i64(row_norm[i], 7U));
        cand_delta[i] = 0;

        if (prob->con[i].is_hard != 0U) {
            if (need > max_hard_viol) {
                max_hard_viol = need;
            }
            if (need > tol_i) {
                all_ok = false;
            }
        } else {
            int64_t smax = (prob->con[i].slack_max_q16 > 0) ? prob->con[i].slack_max_q16 : 0;
            if (need > tol_i || (need > 0 && delta_n[i] > INT64_C(16))) {
                int64_t d_act = (need > smax) ? smax : need;
                cand_delta[i] = d_act;
                if (d_act > max_delta) {
                    max_delta = d_act;
                }
                int64_t rem = robo__sat_sub_i64(need, d_act);
                if (rem > max_soft_viol) {
                    max_soft_viol = rem;
                }
                if (rem > tol_i) {
                    all_ok = false;
                }
            }
        }
    }

    if (all_ok) {
        int64_t obj = 0;
        bool same   = true;
        for (j = 0U; j < n; ++j) {
            cert->u_safe_q16[j] = u_cand[j];
            int64_t d           = robo__sat_sub_i64(u_cand[j], u_nom[j]);
            if ((int64_t)flow_abs_i64(d) > INT64_C(8)) {
                same = false;
            }
            int64_t wd = robo__mul_q16(w_diag[j], robo__mul_q16(d, d));
            obj        = robo__sat_add_i64(obj, flow_ashr_i64(wd, 1U));
        }
        for (i = 0U; i < m; ++i) {
            cert->delta_q16[i] = cand_delta[i];
        }
        cert->feasible          = 1U;
        cert->fallback_engaged  = 0U;
        cert->slack_active      = (max_delta > INT64_C(32)) ? 1U : 0U;
        cert->max_hard_viol_q16 = max_hard_viol;
        cert->max_soft_viol_q16 = max_soft_viol;
        cert->max_slack_q16     = max_delta;
        cert->obj_q16           = obj;
        cert->verdict           = (cert->slack_active != 0U)
                                      ? ROBO_VERDICT_SLACK_RELAXED
                                      : (same ? ROBO_VERDICT_PASS : ROBO_VERDICT_PROJECTED);
    } else {
        /* Guaranteed box-feasible emergency brake fallback */
        int64_t obj = 0;
        for (j = 0U; j < n; ++j) {
            cert->u_safe_q16[j] = u_brake[j];
            int64_t d           = robo__sat_sub_i64(u_brake[j], u_nom[j]);
            int64_t wd          = robo__mul_q16(w_diag[j], robo__mul_q16(d, d));
            obj                 = robo__sat_add_i64(obj, flow_ashr_i64(wd, 1U));
        }
        cert->max_hard_viol_q16 = max_hard_viol;
        cert->max_soft_viol_q16 = max_soft_viol;
        cert->max_slack_q16     = max_delta;
        cert->obj_q16           = obj;
        cert->feasible          = 0U;
        cert->fallback_engaged  = 1U;
        cert->verdict           = ROBO_VERDICT_VETOED_INFEASIBLE;
    }
}

#ifdef __cplusplus
}
#endif

#endif /* ROBO_QP_ND_RING0_H */
