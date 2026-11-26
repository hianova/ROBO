#ifndef ROBO_DYNAMICS_RING0_H
#define ROBO_DYNAMICS_RING0_H

/* ============================================================================
 * ROBO v2.0.0 — Pure-GPR Ring-0 Q16.16 Kinematics, Trig & Contact Limits
 * ============================================================================
 * 100% freestanding Pure-GPR compatible with C17 & C++20:
 *   - Zero FPU / SIMD (-mgeneral-regs-only clean)
 *   - Zero udiv / sdiv / libc calls / heap / indirect jump tables
 * ============================================================================ */

#include "flow.h"

#ifdef __cplusplus
extern "C" {
#endif

#define ROBO_Q16_ONE      FLOW_Q16_ONE
#define ROBO_MAX_OBS      8U
#define ROBO_PI_HALF_Q16  INT64_C(102944)     /* round(pi/2 * 65536) */
#define ROBO_INV_HP_Q32   INT64_C(2734260856) /* round((2/pi) * 2^32) */

typedef struct {
    int64_t ox_q16;
    int64_t oy_q16;
    int64_t radius_q16;
    int64_t ovx_q16;
    int64_t ovy_q16;
} robo_obs_t;

typedef struct {
    /* Planar state (Q16.16 meters, m/s) */
    int64_t x_q16;
    int64_t y_q16;
    int64_t vx_q16;
    int64_t vy_q16;

    /* Workspace safety box (Q16.16 meters) */
    int64_t x_min_q16;
    int64_t x_max_q16;
    int64_t y_min_q16;
    int64_t y_max_q16;

    /* Actuator, robot body & robust sensor/discrete bounds (Q16.16) */
    int64_t u_max_q16;
    int64_t v_max_q16;
    int64_t r_robot_q16;
    int64_t dist_pos_q16;
    int64_t dist_vel_q16;
    int64_t obs_amax_q16; /* worst-case obstacle acceleration bound (m/s^2) */
    int64_t dt_q16;       /* discrete plant sample period (s) */

    /* Contact friction cone & ZMP stability envelope (Q16.16) */
    int64_t fn_q16;         /* normal force N; <= 0 means contact loss */
    int64_t ft_q16;         /* baseline tangential load N */
    int64_t mu_q16;         /* Coulomb friction coeff */
    int64_t mass_q16;       /* effective mass kg */
    int64_t zmp_x_q16;      /* baseline ZMP offset m */
    int64_t zmp_half_w_q16; /* support polygon half-width m */
    int64_t h_over_g_q16;   /* CoM height / gravity (s^2) */

    /* ECBF gains (Q16.16): k1 = 2*lambda, k0 = lambda^2 */
    int64_t k1_q16;
    int64_t k0_q16;

    uint32_t   n_obs;
    robo_obs_t obs[ROBO_MAX_OBS];
} robo_state_t;

FLOW_FORCE_INL int64_t robo__mul_q16(int64_t a, int64_t b) {
    return flow_sym_shr_i128((flow_i128)a * (flow_i128)b, 16U);
}

FLOW_FORCE_INL int64_t robo__div_q16(int64_t num, int64_t den, int64_t cap) {
    return flow_sym_div128_clamped(flow_shl_i128((flow_i128)num, 16U), (flow_i128)den, cap);
}

FLOW_FORCE_INL int64_t robo__sqrt_q16(uint64_t x_q32) {
    return (int64_t)flow_isqrt_round_u128((flow_u128)x_q32);
}

FLOW_FORCE_INL int64_t robo__min_i64(int64_t a, int64_t b) {
    return flow_sel_i64(a < b, a, b);
}

FLOW_FORCE_INL int64_t robo__max_i64(int64_t a, int64_t b) {
    return flow_sel_i64(a > b, a, b);
}

/* Pure-GPR Q16.16 sin(theta) & cos(theta): zero FPU, zero division, Q24 inner precision */
FLOW_INL void robo_sincos_q16(int64_t theta_q16, int64_t *sin_out, int64_t *cos_out) {
    flow_i128 scaled = (flow_i128)theta_q16 * (flow_i128)ROBO_INV_HP_Q32;
    int64_t q        = flow_sym_shr_i128(scaled, 48U);
    flow_i128 rem    = scaled - flow_shl_i128((flow_i128)q, 48U);
    int64_t rem_q24  = flow_sym_shr_i128(rem, 24U);
    int64_t x        = flow_sym_shr_i128((flow_i128)rem_q24 * INT64_C(26353589), 24U);

    int64_t x2 = flow_sym_shr_i128((flow_i128)x * (flow_i128)x, 24U);
    int64_t x3 = flow_sym_shr_i128((flow_i128)x2 * (flow_i128)x, 24U);
    int64_t x4 = flow_sym_shr_i128((flow_i128)x2 * (flow_i128)x2, 24U);
    int64_t x5 = flow_sym_shr_i128((flow_i128)x3 * (flow_i128)x2, 24U);
    int64_t x6 = flow_sym_shr_i128((flow_i128)x4 * (flow_i128)x2, 24U);

    int64_t s0 = flow_sym_shr_i64(
        x - flow_sym_shr_i128((flow_i128)x3 * INT64_C(2796203), 24U) +
            flow_sym_shr_i128((flow_i128)x5 * INT64_C(139810), 24U),
        8U);
    int64_t c0 = flow_sym_shr_i64(
        INT64_C(16777216) - flow_ashr_i64(x2, 1U) +
            flow_sym_shr_i128((flow_i128)x4 * INT64_C(699051), 24U) -
            flow_sym_shr_i128((flow_i128)x6 * INT64_C(23302), 24U),
        8U);

    uint32_t quad = (uint32_t)q & 3U;
    int64_t s_res = s0;
    int64_t c_res = c0;
    if (quad == 1U) {
        s_res = c0;
        c_res = -s0;
    } else if (quad == 2U) {
        s_res = -s0;
        c_res = -c0;
    } else if (quad == 3U) {
        s_res = -c0;
        c_res = s0;
    }
    if (sin_out) {
        *sin_out = s_res;
    }
    if (cos_out) {
        *cos_out = c_res;
    }
}

/* 2-link planar forward kinematics check (l1, l2, theta1, theta2 in Q16.16) */
FLOW_INL void robo_fk2_ring0(int64_t l1_q16, int64_t l2_q16, int64_t th1_q16, int64_t th2_q16,
                             int64_t *ex_q16, int64_t *ey_q16) {
    int64_t s1 = 0, c1 = 0, s12 = 0, c12 = 0;
    robo_sincos_q16(th1_q16, &s1, &c1);
    robo_sincos_q16(th1_q16 + th2_q16, &s12, &c12);
    *ex_q16 = robo__mul_q16(l1_q16, c1) + robo__mul_q16(l2_q16, c12);
    *ey_q16 = robo__mul_q16(l1_q16, s1) + robo__mul_q16(l2_q16, s12);
}

/* Compute exact Q16.16 box [ux_lo, ux_hi] respecting actuation, Coulomb friction, and ZMP */
FLOW_INL bool robo_compute_contact_ux_bounds(const robo_state_t *st, int64_t *ux_lo_out,
                                             int64_t *ux_hi_out, int64_t *fric_slack_out,
                                             int64_t *zmp_slack_out) {
    if (st->fn_q16 <= 0 || st->mass_q16 <= 0) {
        *ux_lo_out      = 0;
        *ux_hi_out      = 0;
        *fric_slack_out = -ROBO_Q16_ONE;
        *zmp_slack_out  = -ROBO_Q16_ONE;
        return false;
    }
    int64_t lo = -st->u_max_q16;
    int64_t hi = st->u_max_q16;

    /* |ft + m * u_x| <= mu * fn */
    int64_t f_lim      = robo__mul_q16(st->mu_q16, st->fn_q16);
    int64_t fric_slack = f_lim - (int64_t)flow_abs_i64(st->ft_q16);
    int64_t u_f_lo     = robo__div_q16(-f_lim - st->ft_q16, st->mass_q16, FLOW_MAX_Q16);
    int64_t u_f_hi     = robo__div_q16(f_lim - st->ft_q16, st->mass_q16, FLOW_MAX_Q16);
    lo                 = robo__max_i64(lo, u_f_lo);
    hi                 = robo__min_i64(hi, u_f_hi);

    /* |zmp_x + (h/g) * u_x| <= zmp_half_w (1-LSB inward margin for exact float64 signoff) */
    int64_t zmp_eff_w = st->zmp_half_w_q16 - INT64_C(64);
    int64_t zmp_slack = zmp_eff_w - (int64_t)flow_abs_i64(st->zmp_x_q16);
    if (st->h_over_g_q16 > 0 && zmp_eff_w > 0) {
        int64_t u_z_lo = robo__div_q16(-zmp_eff_w - st->zmp_x_q16, st->h_over_g_q16, FLOW_MAX_Q16);
        int64_t u_z_hi = robo__div_q16(zmp_eff_w - st->zmp_x_q16, st->h_over_g_q16, FLOW_MAX_Q16);
        lo             = robo__max_i64(lo, u_z_lo);
        hi             = robo__min_i64(hi, u_z_hi);
    }

    *ux_lo_out      = lo;
    *ux_hi_out      = hi;
    *fric_slack_out = fric_slack;
    *zmp_slack_out  = zmp_slack;
    return (lo <= hi) && (fric_slack >= -ROBO_Q16_ONE) && (zmp_slack >= -ROBO_Q16_ONE);
}

#ifdef __cplusplus
}
#endif

#endif /* ROBO_DYNAMICS_RING0_H */
