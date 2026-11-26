#ifndef ROBO_CBF_RING0_H
#define ROBO_CBF_RING0_H

/* ============================================================================
 * ROBO v2.0.0 — Pure-GPR Exact 2D Active-Set QP & Robust Discrete ECBF Filter
 * ============================================================================
 * Deterministic zero-FPU, zero-udiv/sdiv, zero-heap Ring-0 safety shield:
 *   1. Obstacle Second-Order Robust ECBF Halfspaces (with obstacle accel & dt)
 *   2. Workspace Boundary ECBF Halfspaces (4 walls)
 *   3. Exact Contact Friction Cone + ZMP Support Polygon Box Halfspaces
 *   4. Exact Pure-GPR 2D Active-Set QP Solver (single & pairwise Cramer intersections)
 * ============================================================================ */

#include "robo_dynamics_ring0.h"

#ifdef __cplusplus
extern "C" {
#endif

#define ROBO_VERDICT_PASS              0U
#define ROBO_VERDICT_PROJECTED         1U
#define ROBO_VERDICT_VETOED_CONTACT    2U
#define ROBO_VERDICT_VETOED_INFEASIBLE 3U

#define ROBO_MAX_HALFSPACES            16U /* <= 8 obs + 4 walls + 4 box bounds */

typedef struct {
    int64_t ax_q16;
    int64_t ay_q16;
    int64_t b_q16; /* constraint: ax * u_x + ay * u_y >= b */
} robo_halfspace_t;

typedef struct {
    uint32_t verdict;
    uint32_t n_constraints;
    uint32_t active_pairs_checked;
    int64_t  u_safe_x_q16;
    int64_t  u_safe_y_q16;
    int64_t  min_barrier_h_q16;
    int64_t  min_cbf_margin_q16;
    int64_t  friction_slack_q16;
    int64_t  zmp_slack_q16;
    int64_t  intervention_norm_q16;
} robo_safety_cert_t;

FLOW_INL bool robo__feasible_all(const robo_halfspace_t *hs, uint32_t m, int64_t ux, int64_t uy,
                                 int64_t tol_q16) {
    uint32_t k;
    for (k = 0U; k < m; ++k) {
        int64_t dot = robo__mul_q16(hs[k].ax_q16, ux) + robo__mul_q16(hs[k].ay_q16, uy);
        if (dot < hs[k].b_q16 - tol_q16) {
            return false;
        }
    }
    return true;
}

/* Exact 2D Euclidean-Projection QP Solver over <= 16 halfspaces in Pure-GPR */
FLOW_INL bool robo__solve_qp2d_ring0(const robo_halfspace_t *hs, uint32_t m, int64_t ux_nom,
                                     int64_t uy_nom, int64_t *ux_out, int64_t *uy_out,
                                     uint32_t *checks_out) {
    const int64_t tol   = INT64_C(64); /* ~1e-3 Q16 tolerance for rounding safety */
    uint32_t checks     = 1U;
    if (robo__feasible_all(hs, m, ux_nom, uy_nom, tol)) {
        *ux_out = ux_nom;
        *uy_out = uy_nom;
        if (checks_out) *checks_out = checks;
        return true;
    }

    bool found          = false;
    flow_u128 best_d2   = (flow_u128)UINT64_C(0xFFFFFFFFFFFFFFFF) << 60U;
    int64_t best_ux     = 0;
    int64_t best_uy     = 0;
    uint32_t i, j;

    /* 1-active constraint projections */
    for (i = 0U; i < m; ++i) {
        int64_t dot  = robo__mul_q16(hs[i].ax_q16, ux_nom) + robo__mul_q16(hs[i].ay_q16, uy_nom);
        int64_t viol = (hs[i].b_q16 + INT64_C(16)) - dot;
        if (viol <= 0) {
            continue;
        }
        int64_t n2 = robo__mul_q16(hs[i].ax_q16, hs[i].ax_q16) +
                     robo__mul_q16(hs[i].ay_q16, hs[i].ay_q16);
        if (n2 <= INT64_C(16)) {
            continue;
        }
        int64_t step = robo__div_q16(viol, n2, FLOW_MAX_Q16);
        int64_t cx   = ux_nom + robo__mul_q16(step, hs[i].ax_q16);
        int64_t cy   = uy_nom + robo__mul_q16(step, hs[i].ay_q16);
        checks += 1U;
        if (robo__feasible_all(hs, m, cx, cy, tol)) {
            int64_t dx  = cx - ux_nom;
            int64_t dy  = cy - uy_nom;
            flow_u128 d2 = (flow_u128)((flow_i128)dx * (flow_i128)dx +
                                       (flow_i128)dy * (flow_i128)dy);
            if (!found || d2 < best_d2) {
                found   = true;
                best_d2 = d2;
                best_ux = cx;
                best_uy = cy;
            }
        }
    }

    /* 2-active constraint vertex intersections via exact 128-bit Cramer's rule */
    for (i = 0U; i < m; ++i) {
        for (j = i + 1U; j < m; ++j) {
            flow_i128 det = (flow_i128)hs[i].ax_q16 * (flow_i128)hs[j].ay_q16 -
                            (flow_i128)hs[i].ay_q16 * (flow_i128)hs[j].ax_q16;
            flow_u128 adet = flow_abs_i128(det);
            if (adet < (flow_u128)UINT64_C(4294967) || adet > ((flow_u128)UINT64_C(1) << 62U)) {
                continue;
            }
            flow_i128 nx_q32 = (flow_i128)hs[i].b_q16 * (flow_i128)hs[j].ay_q16 -
                               (flow_i128)hs[i].ay_q16 * (flow_i128)hs[j].b_q16;
            flow_i128 ny_q32 = (flow_i128)hs[i].ax_q16 * (flow_i128)hs[j].b_q16 -
                               (flow_i128)hs[i].b_q16 * (flow_i128)hs[j].ax_q16;
            int64_t cx = flow_sym_div128_clamped(flow_shl_i128(nx_q32, 16U), det, FLOW_MAX_Q16);
            int64_t cy = flow_sym_div128_clamped(flow_shl_i128(ny_q32, 16U), det, FLOW_MAX_Q16);
            checks += 1U;
            if (robo__feasible_all(hs, m, cx, cy, tol)) {
                int64_t dx   = cx - ux_nom;
                int64_t dy   = cy - uy_nom;
                flow_u128 d2 = (flow_u128)((flow_i128)dx * (flow_i128)dx +
                                           (flow_i128)dy * (flow_i128)dy);
                if (!found || d2 < best_d2) {
                    found   = true;
                    best_d2 = d2;
                    best_ux = cx;
                    best_uy = cy;
                }
            }
        }
    }

    if (checks_out) *checks_out = checks;
    if (found) {
        *ux_out = best_ux;
        *uy_out = best_uy;
    }
    return found;
}

/* Compute minimum barrier margin & ECBF halfspaces with obstacle-accel + discrete dt safety buffer */
FLOW_INL uint32_t robo__build_ecbf_halfspaces(const robo_state_t *st,
                                              robo_halfspace_t hs[ROBO_MAX_HALFSPACES],
                                              int64_t *min_h_out) {
    uint32_t m    = 0U;
    int64_t min_h = FLOW_MAX_Q16;
    uint32_t i;
    uint32_t n_obs = (st->n_obs > ROBO_MAX_OBS) ? ROBO_MAX_OBS : st->n_obs;

    for (i = 0U; i < n_obs; ++i) {
        int64_t dx   = st->x_q16 - st->obs[i].ox_q16;
        int64_t dy   = st->y_q16 - st->obs[i].oy_q16;
        int64_t dvx  = st->vx_q16 - st->obs[i].ovx_q16;
        int64_t dvy  = st->vy_q16 - st->obs[i].ovy_q16;
        int64_t dt_b = robo__mul_q16(st->v_max_q16, st->dt_q16);
        int64_t rsum = st->obs[i].radius_q16 + st->r_robot_q16 + st->dist_pos_q16 + dt_b;

        flow_u128 d2_u128 = (flow_u128)((flow_i128)dx * (flow_i128)dx +
                                        (flow_i128)dy * (flow_i128)dy);
        int64_t d2_q16    = flow_sym_shr_i128((flow_i128)d2_u128, 16U);
        int64_t r2_q16    = robo__mul_q16(rsum, rsum);
        int64_t h_q16     = d2_q16 - r2_q16;
        if (h_q16 < min_h) {
            min_h = h_q16;
        }

        int64_t dist_q16  = robo__sqrt_q16((uint64_t)d2_u128);
        int64_t hdot_q16  = flow_sym_shr_i128((flow_i128)dx * (flow_i128)dvx +
                                              (flow_i128)dy * (flow_i128)dvy,
                                              15U) -
                            flow_sym_shr_i128((flow_i128)dist_q16 * (flow_i128)st->dist_vel_q16, 15U);
        int64_t dv2_2_q16 = flow_sym_shr_i128((flow_i128)dvx * (flow_i128)dvx +
                                              (flow_i128)dvy * (flow_i128)dvy,
                                              15U);
        /* Subtract worst-case moving-obstacle accel + discrete step curvature buffer */
        int64_t obs_buf   = flow_sym_shr_i128((flow_i128)dist_q16 * (flow_i128)st->obs_amax_q16, 15U);
        int64_t disc_buf  = robo__mul_q16(flow_shl_i64(st->dt_q16, 3U), robo__mul_q16(st->u_max_q16, st->u_max_q16));

        int64_t rhs_pos   = dv2_2_q16 + robo__mul_q16(st->k1_q16, hdot_q16) +
                            robo__mul_q16(st->k0_q16, h_q16) - obs_buf - disc_buf;
        hs[m].ax_q16      = flow_shl_i64(dx, 1U);
        hs[m].ay_q16      = flow_shl_i64(dy, 1U);
        hs[m].b_q16       = -rhs_pos;
        m += 1U;
    }

    /* 4 Workspace Wall ECBFs */
    {
        int64_t wall_buf = st->dist_pos_q16 + robo__mul_q16(st->v_max_q16, st->dt_q16);
        int64_t dv       = st->dist_vel_q16;

        int64_t hx_lo = (st->x_q16 - st->x_min_q16) - wall_buf;
        int64_t hx_hi = (st->x_max_q16 - st->x_q16) - wall_buf;
        int64_t hy_lo = (st->y_q16 - st->y_min_q16) - wall_buf;
        int64_t hy_hi = (st->y_max_q16 - st->y_q16) - wall_buf;

        min_h = robo__min_i64(min_h, hx_lo);
        min_h = robo__min_i64(min_h, hx_hi);
        min_h = robo__min_i64(min_h, hy_lo);
        min_h = robo__min_i64(min_h, hy_hi);

        hs[m].ax_q16 = ROBO_Q16_ONE;
        hs[m].ay_q16 = 0;
        hs[m].b_q16  = -(robo__mul_q16(st->k1_q16, st->vx_q16 - dv) + robo__mul_q16(st->k0_q16, hx_lo));
        m += 1U;

        hs[m].ax_q16 = -ROBO_Q16_ONE;
        hs[m].ay_q16 = 0;
        hs[m].b_q16  = -(robo__mul_q16(st->k1_q16, -st->vx_q16 - dv) + robo__mul_q16(st->k0_q16, hx_hi));
        m += 1U;

        hs[m].ax_q16 = 0;
        hs[m].ay_q16 = ROBO_Q16_ONE;
        hs[m].b_q16  = -(robo__mul_q16(st->k1_q16, st->vy_q16 - dv) + robo__mul_q16(st->k0_q16, hy_lo));
        m += 1U;

        hs[m].ax_q16 = 0;
        hs[m].ay_q16 = -ROBO_Q16_ONE;
        hs[m].b_q16  = -(robo__mul_q16(st->k1_q16, -st->vy_q16 - dv) + robo__mul_q16(st->k0_q16, hy_hi));
        m += 1U;
    }

    *min_h_out = min_h;
    return m;
}

FLOW_INL robo_safety_cert_t robo_cbf_filter_ring0(const robo_state_t *st, int64_t u_nom_x_q16,
                                                  int64_t u_nom_y_q16) {
    robo_safety_cert_t cert;
    cert.verdict               = ROBO_VERDICT_PASS;
    cert.n_constraints         = 0U;
    cert.active_pairs_checked  = 0U;
    cert.u_safe_x_q16          = 0;
    cert.u_safe_y_q16          = 0;
    cert.min_barrier_h_q16     = 0;
    cert.min_cbf_margin_q16    = 0;
    cert.friction_slack_q16    = 0;
    cert.zmp_slack_q16         = 0;
    cert.intervention_norm_q16 = 0;

    if (!st) {
        cert.verdict = ROBO_VERDICT_VETOED_INFEASIBLE;
        return cert;
    }

    int64_t ux_lo = 0, ux_hi = 0, fric_slack = 0, zmp_slack = 0;
    bool contact_ok = robo_compute_contact_ux_bounds(st, &ux_lo, &ux_hi, &fric_slack, &zmp_slack);
    cert.friction_slack_q16 = fric_slack;
    cert.zmp_slack_q16      = zmp_slack;

    robo_halfspace_t hs[ROBO_MAX_HALFSPACES];
    int64_t min_h          = 0;
    uint32_t m_ecbf        = robo__build_ecbf_halfspaces(st, hs, &min_h);
    cert.min_barrier_h_q16 = min_h;

    if (!contact_ok) {
        cert.n_constraints      = m_ecbf;
        cert.verdict            = ROBO_VERDICT_VETOED_CONTACT;
        cert.u_safe_x_q16       = 0;
        cert.u_safe_y_q16       = 0;
        cert.min_cbf_margin_q16 = -ROBO_Q16_ONE;
        return cert;
    }

    int64_t uy_lo = -st->u_max_q16;
    int64_t uy_hi = st->u_max_q16;

    /* Append exact box bounds as linear halfspaces so Active-Set QP handles all bounds jointly */
    uint32_t m = m_ecbf;
    hs[m++]    = (robo_halfspace_t){ROBO_Q16_ONE, 0, ux_lo};
    hs[m++]    = (robo_halfspace_t){-ROBO_Q16_ONE, 0, -ux_hi};
    hs[m++]    = (robo_halfspace_t){0, ROBO_Q16_ONE, uy_lo};
    hs[m++]    = (robo_halfspace_t){0, -ROBO_Q16_ONE, -uy_hi};
    cert.n_constraints = m;

    int64_t ux = 0, uy = 0;
    uint32_t checks = 0U;
    bool ok         = robo__solve_qp2d_ring0(hs, m, u_nom_x_q16, u_nom_y_q16, &ux, &uy, &checks);
    cert.active_pairs_checked = checks;

    if (!ok) {
        int64_t bx   = flow_shl_i64(-st->vx_q16, 2U);
        int64_t by   = flow_shl_i64(-st->vy_q16, 2U);
        ux           = flow_sel_i64(bx < ux_lo, ux_lo, flow_sel_i64(bx > ux_hi, ux_hi, bx));
        uy           = flow_sel_i64(by < uy_lo, uy_lo, flow_sel_i64(by > uy_hi, uy_hi, by));
        cert.verdict = ROBO_VERDICT_VETOED_INFEASIBLE;
    } else {
        ux           = flow_sel_i64(ux < ux_lo, ux_lo, flow_sel_i64(ux > ux_hi, ux_hi, ux));
        uy           = flow_sel_i64(uy < uy_lo, uy_lo, flow_sel_i64(uy > uy_hi, uy_hi, uy));
        bool same    = (flow_abs_i64(ux - u_nom_x_q16) <= INT64_C(8)) &&
                       (flow_abs_i64(uy - u_nom_y_q16) <= INT64_C(8));
        cert.verdict = same ? ROBO_VERDICT_PASS : ROBO_VERDICT_PROJECTED;
    }

    int64_t min_margin = FLOW_MAX_Q16;
    uint32_t j;
    for (j = 0U; j < m_ecbf; ++j) {
        int64_t dot  = robo__mul_q16(hs[j].ax_q16, ux) + robo__mul_q16(hs[j].ay_q16, uy);
        int64_t marg = dot - hs[j].b_q16;
        if (marg < min_margin) {
            min_margin = marg;
        }
    }

    cert.u_safe_x_q16       = ux;
    cert.u_safe_y_q16       = uy;
    cert.min_cbf_margin_q16 = min_margin;

    int64_t dxu   = ux - u_nom_x_q16;
    int64_t dyu   = uy - u_nom_y_q16;
    flow_u128 du2 = (flow_u128)((flow_i128)dxu * (flow_i128)dxu +
                                (flow_i128)dyu * (flow_i128)dyu);
    cert.intervention_norm_q16 = robo__sqrt_q16((uint64_t)du2);
    return cert;
}

#ifdef __cplusplus
}
#endif

#endif /* ROBO_CBF_RING0_H */
