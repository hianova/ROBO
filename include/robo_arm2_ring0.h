#ifndef ROBO_ARM2_RING0_H
#define ROBO_ARM2_RING0_H

/* ============================================================================
 * ROBO v2.0.0 — Pure-GPR 2-Link Rigid-Body Manipulator Torque CBF Filter
 * ============================================================================
 * Exact analytical M(q) * ddq + C(q, dq) * dq + g(q) = tau mapped through J(q)
 * into torque-space ECBF halfspaces & joint-angle ECBFs, solved via exact 2D QP.
 * Zero FPU, zero udiv/sdiv, zero heap.
 * ============================================================================ */

#include "robo_cbf_ring0.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int64_t q1_q16;
    int64_t q2_q16;
    int64_t dq1_q16;
    int64_t dq2_q16;

    /* Physical rigid-body parameters (Q16.16 SI) */
    int64_t l1_q16;
    int64_t l2_q16;
    int64_t m1_q16;
    int64_t m2_q16;
    int64_t i1_q16;
    int64_t i2_q16;
    int64_t grav_q16; /* 0 for horizontal planar, 9.81 for vertical */

    /* Joint & torque limits + ECBF gains */
    int64_t q1_min_q16, q1_max_q16;
    int64_t q2_min_q16, q2_max_q16;
    int64_t tau_max_q16;
    int64_t r_ee_q16;
    int64_t k1_q16;
    int64_t k0_q16;

    uint32_t   n_obs;
    robo_obs_t obs[ROBO_MAX_OBS];
} robo_arm2_state_t;

typedef struct {
    uint32_t verdict;
    uint32_t n_constraints;
    int64_t  tau1_safe_q16;
    int64_t  tau2_safe_q16;
    int64_t  ee_x_q16;
    int64_t  ee_y_q16;
    int64_t  min_h_q16;
} robo_arm2_cert_t;

/* Exact Pure-GPR 2-Link Rigid-Body Forward Dynamics: ddq = M(q)^-1 (tau - C*dq - g) */
FLOW_INL void robo_arm2_forward_dyn_ring0(const robo_arm2_state_t *st, int64_t tau1_q16,
                                          int64_t tau2_q16, int64_t *ddq1_out, int64_t *ddq2_out) {
    int64_t s1 = 0, c1 = 0, s2 = 0, c2 = 0, s12 = 0, c12 = 0;
    robo_sincos_q16(st->q1_q16, &s1, &c1);
    robo_sincos_q16(st->q2_q16, &s2, &c2);
    robo_sincos_q16(st->q1_q16 + st->q2_q16, &s12, &c12);
    (void)s1; (void)s12;

    int64_t lc1 = flow_ashr_i64(st->l1_q16, 1U);
    int64_t lc2 = flow_ashr_i64(st->l2_q16, 1U);
    int64_t a1  = st->i1_q16 + st->i2_q16 +
                  robo__mul_q16(st->m1_q16, robo__mul_q16(lc1, lc1)) +
                  robo__mul_q16(st->m2_q16, robo__mul_q16(st->l1_q16, st->l1_q16) +
                                                robo__mul_q16(lc2, lc2));
    int64_t a2  = robo__mul_q16(st->m2_q16, robo__mul_q16(st->l1_q16, lc2));
    int64_t a3  = st->i2_q16 + robo__mul_q16(st->m2_q16, robo__mul_q16(lc2, lc2));

    int64_t m11 = a1 + flow_shl_i64(robo__mul_q16(a2, c2), 1U);
    int64_t m12 = a3 + robo__mul_q16(a2, c2);
    int64_t m22 = a3;

    int64_t h_s2 = robo__mul_q16(a2, s2);
    int64_t cvec1 = -robo__mul_q16(h_s2, flow_shl_i64(robo__mul_q16(st->dq1_q16, st->dq2_q16), 1U) +
                                            robo__mul_q16(st->dq2_q16, st->dq2_q16));
    int64_t cvec2 = robo__mul_q16(h_s2, robo__mul_q16(st->dq1_q16, st->dq1_q16));

    int64_t g2 = robo__mul_q16(robo__mul_q16(st->m2_q16, lc2), robo__mul_q16(st->grav_q16, c12));
    int64_t g1 = robo__mul_q16(robo__mul_q16(st->m1_q16, lc1) + robo__mul_q16(st->m2_q16, st->l1_q16),
                               robo__mul_q16(st->grav_q16, c1)) + g2;

    int64_t rhs1 = tau1_q16 - cvec1 - g1;
    int64_t rhs2 = tau2_q16 - cvec2 - g2;

    flow_i128 det = (flow_i128)m11 * (flow_i128)m22 - (flow_i128)m12 * (flow_i128)m12;
    flow_i128 n1  = (flow_i128)m22 * (flow_i128)rhs1 - (flow_i128)m12 * (flow_i128)rhs2;
    flow_i128 n2  = (flow_i128)m11 * (flow_i128)rhs2 - (flow_i128)m12 * (flow_i128)rhs1;

    *ddq1_out = flow_sym_div128_clamped(flow_shl_i128(n1, 16U), det, FLOW_MAX_Q16);
    *ddq2_out = flow_sym_div128_clamped(flow_shl_i128(n2, 16U), det, FLOW_MAX_Q16);
}

FLOW_INL robo_arm2_cert_t robo_arm2_cbf_ring0(const robo_arm2_state_t *st, int64_t tau1_nom_q16,
                                              int64_t tau2_nom_q16) {
    robo_arm2_cert_t cert;
    cert.verdict       = ROBO_VERDICT_PASS;
    cert.n_constraints = 0U;
    cert.tau1_safe_q16 = tau1_nom_q16;
    cert.tau2_safe_q16 = tau2_nom_q16;
    cert.ee_x_q16      = 0;
    cert.ee_y_q16      = 0;
    cert.min_h_q16     = FLOW_MAX_Q16;
    if (!st) {
        cert.verdict = ROBO_VERDICT_VETOED_INFEASIBLE;
        return cert;
    }

    int64_t s1 = 0, c1 = 0, s2 = 0, c2 = 0, s12 = 0, c12 = 0;
    robo_sincos_q16(st->q1_q16, &s1, &c1);
    robo_sincos_q16(st->q2_q16, &s2, &c2);
    robo_sincos_q16(st->q1_q16 + st->q2_q16, &s12, &c12);

    /* Forward kinematics & Jacobian J(q) */
    int64_t l1c1  = robo__mul_q16(st->l1_q16, c1);
    int64_t l1s1  = robo__mul_q16(st->l1_q16, s1);
    int64_t l2c12 = robo__mul_q16(st->l2_q16, c12);
    int64_t l2s12 = robo__mul_q16(st->l2_q16, s12);
    int64_t px    = l1c1 + l2c12;
    int64_t py    = l1s1 + l2s12;
    cert.ee_x_q16 = px;
    cert.ee_y_q16 = py;

    int64_t j00 = -l1s1 - l2s12, j01 = -l2s12;
    int64_t j10 = l1c1 + l2c12,  j11 = l2c12;
    int64_t vx  = robo__mul_q16(j00, st->dq1_q16) + robo__mul_q16(j01, st->dq2_q16);
    int64_t vy  = robo__mul_q16(j10, st->dq1_q16) + robo__mul_q16(j11, st->dq2_q16);

    int64_t dq12  = st->dq1_q16 + st->dq2_q16;
    int64_t dq1sq = robo__mul_q16(st->dq1_q16, st->dq1_q16);
    int64_t dq12sq = robo__mul_q16(dq12, dq12);
    int64_t jdot_x = -robo__mul_q16(l1c1, dq1sq) - robo__mul_q16(l2c12, dq12sq);
    int64_t jdot_y = -robo__mul_q16(l1s1, dq1sq) - robo__mul_q16(l2s12, dq12sq);

    /* Compute M^-1 and M^-1 * h_cg */
    int64_t lc1 = flow_ashr_i64(st->l1_q16, 1U);
    int64_t lc2 = flow_ashr_i64(st->l2_q16, 1U);
    int64_t a1  = st->i1_q16 + st->i2_q16 +
                  robo__mul_q16(st->m1_q16, robo__mul_q16(lc1, lc1)) +
                  robo__mul_q16(st->m2_q16, robo__mul_q16(st->l1_q16, st->l1_q16) +
                                                robo__mul_q16(lc2, lc2));
    int64_t a2  = robo__mul_q16(st->m2_q16, robo__mul_q16(st->l1_q16, lc2));
    int64_t a3  = st->i2_q16 + robo__mul_q16(st->m2_q16, robo__mul_q16(lc2, lc2));

    int64_t m11 = a1 + flow_shl_i64(robo__mul_q16(a2, c2), 1U);
    int64_t m12 = a3 + robo__mul_q16(a2, c2);
    int64_t m22 = a3;
    flow_i128 det = (flow_i128)m11 * (flow_i128)m22 - (flow_i128)m12 * (flow_i128)m12;

    int64_t inv00 = flow_sym_div128_clamped(flow_shl_i128((flow_i128)m22, 32U), det, FLOW_MAX_Q16);
    int64_t inv01 = flow_sym_div128_clamped(flow_shl_i128(-(flow_i128)m12, 32U), det, FLOW_MAX_Q16);
    int64_t inv11 = flow_sym_div128_clamped(flow_shl_i128((flow_i128)m11, 32U), det, FLOW_MAX_Q16);

    int64_t h_s2  = robo__mul_q16(a2, s2);
    int64_t cvec1 = -robo__mul_q16(h_s2, flow_shl_i64(robo__mul_q16(st->dq1_q16, st->dq2_q16), 1U) +
                                             robo__mul_q16(st->dq2_q16, st->dq2_q16));
    int64_t cvec2 = robo__mul_q16(h_s2, dq1sq);
    int64_t g2    = robo__mul_q16(robo__mul_q16(st->m2_q16, lc2), robo__mul_q16(st->grav_q16, c12));
    int64_t g1    = robo__mul_q16(robo__mul_q16(st->m1_q16, lc1) + robo__mul_q16(st->m2_q16, st->l1_q16),
                                  robo__mul_q16(st->grav_q16, c1)) + g2;
    int64_t hcg1  = cvec1 + g1;
    int64_t hcg2  = cvec2 + g2;

    int64_t d0_1 = -(robo__mul_q16(inv00, hcg1) + robo__mul_q16(inv01, hcg2));
    int64_t d0_2 = -(robo__mul_q16(inv01, hcg1) + robo__mul_q16(inv11, hcg2));

    /* G = J * M^-1 ; b0 = J * (-M^-1 * hcg) + Jdot * dq */
    int64_t g00 = robo__mul_q16(j00, inv00) + robo__mul_q16(j01, inv01);
    int64_t g01 = robo__mul_q16(j00, inv01) + robo__mul_q16(j01, inv11);
    int64_t g10 = robo__mul_q16(j10, inv00) + robo__mul_q16(j11, inv01);
    int64_t g11 = robo__mul_q16(j10, inv01) + robo__mul_q16(j11, inv11);

    int64_t b0x = robo__mul_q16(j00, d0_1) + robo__mul_q16(j01, d0_2) + jdot_x;
    int64_t b0y = robo__mul_q16(j10, d0_1) + robo__mul_q16(j11, d0_2) + jdot_y;

    robo_halfspace_t hs[ROBO_MAX_HALFSPACES];
    uint32_t m    = 0U;
    int64_t min_h = FLOW_MAX_Q16;
    uint32_t n_obs = (st->n_obs > ROBO_MAX_OBS) ? ROBO_MAX_OBS : st->n_obs;
    uint32_t i;

    for (i = 0U; i < n_obs && m < 8U; ++i) {
        int64_t dx   = px - st->obs[i].ox_q16;
        int64_t dy   = py - st->obs[i].oy_q16;
        int64_t rsum = st->obs[i].radius_q16 + st->r_ee_q16;
        flow_u128 d2 = (flow_u128)((flow_i128)dx * (flow_i128)dx + (flow_i128)dy * (flow_i128)dy);
        int64_t h    = flow_sym_shr_i128((flow_i128)d2, 16U) - robo__mul_q16(rsum, rsum);
        if (h < min_h) min_h = h;

        int64_t hdot  = flow_sym_shr_i128((flow_i128)dx * (flow_i128)vx + (flow_i128)dy * (flow_i128)vy, 15U);
        int64_t v2_2  = flow_sym_shr_i128((flow_i128)vx * (flow_i128)vx + (flow_i128)vy * (flow_i128)vy, 15U);
        int64_t ax_c  = flow_shl_i64(dx, 1U);
        int64_t ay_c  = flow_shl_i64(dy, 1U);
        int64_t b_c   = -(v2_2 + robo__mul_q16(st->k1_q16, hdot) + robo__mul_q16(st->k0_q16, h) - INT64_C(32768));

        /* Map ax_c*ddp_x + ay_c*ddp_y >= b_c  ->  (a^T G) tau >= b_c - a^T b0 */
        hs[m].ax_q16 = robo__mul_q16(ax_c, g00) + robo__mul_q16(ay_c, g10);
        hs[m].ay_q16 = robo__mul_q16(ax_c, g01) + robo__mul_q16(ay_c, g11);
        hs[m].b_q16  = b_c - (robo__mul_q16(ax_c, b0x) + robo__mul_q16(ay_c, b0y));
        m += 1U;
    }

    /* Joint angle limits via ECBF mapped through M^-1 */
    {
        int64_t hq1_lo = st->q1_q16 - st->q1_min_q16;
        int64_t hq1_hi = st->q1_max_q16 - st->q1_q16;
        int64_t hq2_lo = st->q2_q16 - st->q2_min_q16;
        int64_t hq2_hi = st->q2_max_q16 - st->q2_q16;
        min_h          = robo__min_i64(min_h, robo__min_i64(hq1_lo, hq1_hi));
        min_h          = robo__min_i64(min_h, robo__min_i64(hq2_lo, hq2_hi));

        int64_t bq1_lo = -(robo__mul_q16(st->k1_q16, st->dq1_q16) + robo__mul_q16(st->k0_q16, hq1_lo));
        int64_t bq1_hi = -(robo__mul_q16(st->k1_q16, -st->dq1_q16) + robo__mul_q16(st->k0_q16, hq1_hi));
        int64_t bq2_lo = -(robo__mul_q16(st->k1_q16, st->dq2_q16) + robo__mul_q16(st->k0_q16, hq2_lo));
        int64_t bq2_hi = -(robo__mul_q16(st->k1_q16, -st->dq2_q16) + robo__mul_q16(st->k0_q16, hq2_hi));

        hs[m++] = (robo_halfspace_t){inv00, inv01, bq1_lo - d0_1};
        hs[m++] = (robo_halfspace_t){-inv00, -inv01, bq1_hi + d0_1};
        hs[m++] = (robo_halfspace_t){inv01, inv11, bq2_lo - d0_2};
        hs[m++] = (robo_halfspace_t){-inv01, -inv11, bq2_hi + d0_2};
    }

    /* Torque box bounds [-tau_max, tau_max]^2 */
    hs[m++] = (robo_halfspace_t){ROBO_Q16_ONE, 0, -st->tau_max_q16};
    hs[m++] = (robo_halfspace_t){-ROBO_Q16_ONE, 0, -st->tau_max_q16};
    hs[m++] = (robo_halfspace_t){0, ROBO_Q16_ONE, -st->tau_max_q16};
    hs[m++] = (robo_halfspace_t){0, -ROBO_Q16_ONE, -st->tau_max_q16};

    cert.n_constraints = m;
    cert.min_h_q16     = min_h;

    int64_t t1 = 0, t2 = 0;
    uint32_t checks = 0U;
    bool ok         = robo__solve_qp2d_ring0(hs, m, tau1_nom_q16, tau2_nom_q16, &t1, &t2, &checks);
    if (!ok) {
        t1           = flow_clamp_sym(-flow_shl_i64(st->dq1_q16, 4U), st->tau_max_q16);
        t2           = flow_clamp_sym(-flow_shl_i64(st->dq2_q16, 4U), st->tau_max_q16);
        cert.verdict = ROBO_VERDICT_VETOED_INFEASIBLE;
    } else {
        t1           = flow_clamp_sym(t1, st->tau_max_q16);
        t2           = flow_clamp_sym(t2, st->tau_max_q16);
        bool same    = (flow_abs_i64(t1 - tau1_nom_q16) <= INT64_C(8)) &&
                       (flow_abs_i64(t2 - tau2_nom_q16) <= INT64_C(8));
        cert.verdict = same ? ROBO_VERDICT_PASS : ROBO_VERDICT_PROJECTED;
    }
    cert.tau1_safe_q16 = t1;
    cert.tau2_safe_q16 = t2;
    return cert;
}

#ifdef __cplusplus
}
#endif

#endif /* ROBO_ARM2_RING0_H */
