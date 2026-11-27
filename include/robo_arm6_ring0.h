#ifndef ROBO_ARM6_RING0_H
#define ROBO_ARM6_RING0_H

/* ============================================================================
 * ROBO v3.0.0 — Pure-GPR 3D 6-DoF Spatial Manipulator + Latency Tube-CBF
 * ============================================================================
 * Freestanding Pure-GPR Q16.16 spatial kinematics, SO(3)/quaternion algebra,
 * analytical 6-DoF Jacobian & Jdot*dq, exact 3D capsule-to-sphere &
 * capsule-to-capsule self-collision distances/gradients, rigid-body M(q)/g(q),
 * table contact-force barrier, and explicit 30-50 ms VLA delay Tube-CBF.
 *   - Zero FPU / SIMD (-mgeneral-regs-only clean)
 *   - Zero udiv / sdiv / libc / heap / indirect branches
 *   - Zero signed-overflow / shift UB
 * ============================================================================ */

#include "robo_qp_nd_ring0.h"

#ifdef __cplusplus
extern "C" {
#endif

#define ROBO_ARM6_DOF        6U
#define ROBO_ARM6_MAX_OBS    8U
#define ROBO_ARM6_AXIS_X     0U
#define ROBO_ARM6_AXIS_Y     1U
#define ROBO_ARM6_AXIS_Z     2U

typedef struct {
    int64_t x;
    int64_t y;
    int64_t z;
} robo_vec3_t;

typedef struct {
    int64_t w;
    int64_t x;
    int64_t y;
    int64_t z;
} robo_quat_t;

typedef struct {
    int64_t m[3][3];
} robo_mat3_t;

typedef struct {
    int64_t ox_q16;
    int64_t oy_q16;
    int64_t oz_q16;
    int64_t radius_q16;
    int64_t ovx_q16;
    int64_t ovy_q16;
    int64_t ovz_q16;
} robo_obs3_t;

typedef struct {
    robo_vec3_t offset_q16;
    uint32_t    axis; /* 0=X, 1=Y, 2=Z */
    robo_vec3_t cap_a_local_q16;
    robo_vec3_t cap_b_local_q16;
    int64_t     cap_radius_q16;
    int64_t     mass_q16;
    int64_t     inertia_q16;
    robo_vec3_t com_local_q16;
} robo_arm6_link_t;

typedef struct {
    int64_t          q_q16[ROBO_ARM6_DOF];
    int64_t          dq_q16[ROBO_ARM6_DOF];
    robo_arm6_link_t links[ROBO_ARM6_DOF];
    robo_vec3_t      ee_offset_q16;
    int64_t          ee_radius_q16;
    int64_t          payload_mass_q16;
    int64_t          grav_q16;

    /* Joint, velocity, acceleration, torque & contact limits */
    int64_t q_min_q16[ROBO_ARM6_DOF];
    int64_t q_max_q16[ROBO_ARM6_DOF];
    int64_t dq_max_q16[ROBO_ARM6_DOF];
    int64_t ddq_max_q16[ROBO_ARM6_DOF];
    int64_t tau_max_q16[ROBO_ARM6_DOF];

    int64_t table_z_q16;
    int64_t contact_k_q16;
    int64_t contact_b_q16;
    int64_t contact_fmax_q16;

    /* Explicit VLA action-delay Tube-CBF & uncertainty bounds */
    int64_t tau_delay_q16; /* 0.030 - 0.050 s */
    int64_t dt_q16;        /* 0.001 s (1 kHz shield) */
    int64_t dist_pos_q16;
    int64_t dist_vel_q16;
    int64_t dist_acc_q16;
    int64_t obs_amax_q16;

    /* ECBF gains */
    int64_t k1_q16;
    int64_t k0_q16;
    int64_t kv_q16;

    uint32_t    n_obs;
    robo_obs3_t obs[ROBO_ARM6_MAX_OBS];
} robo_arm6_state_t;

typedef struct {
    robo_vec3_t joint_pos[ROBO_ARM6_DOF];
    robo_vec3_t joint_vel[ROBO_ARM6_DOF];
    robo_vec3_t joint_axis[ROBO_ARM6_DOF];
    robo_vec3_t omega[ROBO_ARM6_DOF];
    robo_quat_t link_quat[ROBO_ARM6_DOF];
    robo_mat3_t link_rot[ROBO_ARM6_DOF];
    robo_vec3_t cap_a[ROBO_ARM6_DOF];
    robo_vec3_t cap_b[ROBO_ARM6_DOF];
    robo_vec3_t com_world[ROBO_ARM6_DOF];
    robo_vec3_t ee_pos;
    robo_vec3_t ee_vel;
    robo_quat_t ee_quat;
    int64_t     J_ee[3][ROBO_ARM6_DOF];
    int64_t     J_rot[3][ROBO_ARM6_DOF];
} robo_arm6_fk_t;

typedef struct {
    uint32_t verdict;
    uint32_t n_constraints;
    uint32_t qp_iters;
    uint8_t  feasible;
    uint8_t  slack_active;
    uint8_t  fallback_engaged;
    uint8_t  mass_sdd_certified;
    int64_t  ddq_safe_q16[ROBO_ARM6_DOF];
    int64_t  tau_safe_q16[ROBO_ARM6_DOF];
    robo_vec3_t ee_pos;
    robo_quat_t ee_quat;
    int64_t  min_obs_h_q16;
    int64_t  min_self_h_q16;
    int64_t  min_contact_h_q16;
    int64_t  min_joint_h_q16;
    int64_t  min_stl_h_q16;
    int64_t  tube_margin_q16;
    int64_t  max_slack_q16;
    int64_t  mass_sdd_mu_q32;
} robo_arm6_cert_t;

/* --- 3D Vector, SO(3) & Unit-Quaternion Pure-GPR Q16.16 Algebra --- */

FLOW_FORCE_INL robo_vec3_t robo_v3_make(int64_t x, int64_t y, int64_t z) {
    robo_vec3_t v;
    v.x = x;
    v.y = y;
    v.z = z;
    return v;
}

FLOW_FORCE_INL robo_vec3_t robo_v3_add(robo_vec3_t a, robo_vec3_t b) {
    return robo_v3_make(robo__sat_add_i64(a.x, b.x), robo__sat_add_i64(a.y, b.y),
                        robo__sat_add_i64(a.z, b.z));
}

FLOW_FORCE_INL robo_vec3_t robo_v3_sub(robo_vec3_t a, robo_vec3_t b) {
    return robo_v3_make(robo__sat_sub_i64(a.x, b.x), robo__sat_sub_i64(a.y, b.y),
                        robo__sat_sub_i64(a.z, b.z));
}

FLOW_FORCE_INL robo_vec3_t robo_v3_scale(robo_vec3_t a, int64_t s) {
    return robo_v3_make(robo__mul_q16(a.x, s), robo__mul_q16(a.y, s), robo__mul_q16(a.z, s));
}

FLOW_FORCE_INL int64_t robo_v3_dot(robo_vec3_t a, robo_vec3_t b) {
    flow_i128 d = (flow_i128)a.x * (flow_i128)b.x + (flow_i128)a.y * (flow_i128)b.y +
                  (flow_i128)a.z * (flow_i128)b.z;
    return flow_sym_shr_i128(d, 16U);
}

FLOW_FORCE_INL flow_u128 robo_v3_norm2_q32(robo_vec3_t a) {
    return (flow_u128)((flow_i128)a.x * (flow_i128)a.x + (flow_i128)a.y * (flow_i128)a.y +
                       (flow_i128)a.z * (flow_i128)a.z);
}

FLOW_FORCE_INL int64_t robo_v3_norm(robo_vec3_t a) {
    return robo__sqrt_q16((uint64_t)robo_v3_norm2_q32(a));
}

FLOW_FORCE_INL robo_vec3_t robo_v3_cross(robo_vec3_t a, robo_vec3_t b) {
    flow_i128 cx = (flow_i128)a.y * (flow_i128)b.z - (flow_i128)a.z * (flow_i128)b.y;
    flow_i128 cy = (flow_i128)a.z * (flow_i128)b.x - (flow_i128)a.x * (flow_i128)b.z;
    flow_i128 cz = (flow_i128)a.x * (flow_i128)b.y - (flow_i128)a.y * (flow_i128)b.x;
    return robo_v3_make(flow_sym_shr_i128(cx, 16U), flow_sym_shr_i128(cy, 16U),
                        flow_sym_shr_i128(cz, 16U));
}

FLOW_INL robo_vec3_t robo_v3_unit(robo_vec3_t a, int64_t *norm_out) {
    int64_t n = robo_v3_norm(a);
    if (norm_out) {
        *norm_out = n;
    }
    if (n <= INT64_C(16)) {
        return robo_v3_make(0, 0, ROBO_Q16_ONE);
    }
    return robo_v3_make(
        flow_sym_div128_clamped(flow_shl_i128((flow_i128)a.x, 16U), (flow_i128)n, ROBO_Q16_ONE),
        flow_sym_div128_clamped(flow_shl_i128((flow_i128)a.y, 16U), (flow_i128)n, ROBO_Q16_ONE),
        flow_sym_div128_clamped(flow_shl_i128((flow_i128)a.z, 16U), (flow_i128)n, ROBO_Q16_ONE));
}

FLOW_FORCE_INL robo_quat_t robo_quat_identity(void) {
    robo_quat_t q;
    q.w = ROBO_Q16_ONE;
    q.x = 0;
    q.y = 0;
    q.z = 0;
    return q;
}

FLOW_INL robo_quat_t robo_quat_normalize(robo_quat_t q) {
    flow_u128 n2 = (flow_u128)((flow_i128)q.w * (flow_i128)q.w + (flow_i128)q.x * (flow_i128)q.x +
                               (flow_i128)q.y * (flow_i128)q.y + (flow_i128)q.z * (flow_i128)q.z);
    int64_t n    = robo__sqrt_q16((uint64_t)n2);
    if (n <= INT64_C(16)) {
        return robo_quat_identity();
    }
    robo_quat_t out;
    out.w = flow_sym_div128_clamped(flow_shl_i128((flow_i128)q.w, 16U), (flow_i128)n, ROBO_Q16_ONE);
    out.x = flow_sym_div128_clamped(flow_shl_i128((flow_i128)q.x, 16U), (flow_i128)n, ROBO_Q16_ONE);
    out.y = flow_sym_div128_clamped(flow_shl_i128((flow_i128)q.y, 16U), (flow_i128)n, ROBO_Q16_ONE);
    out.z = flow_sym_div128_clamped(flow_shl_i128((flow_i128)q.z, 16U), (flow_i128)n, ROBO_Q16_ONE);
    return out;
}

FLOW_INL robo_quat_t robo_quat_mul(robo_quat_t a, robo_quat_t b) {
    robo_quat_t r;
    flow_i128 w = (flow_i128)a.w * (flow_i128)b.w - (flow_i128)a.x * (flow_i128)b.x -
                  (flow_i128)a.y * (flow_i128)b.y - (flow_i128)a.z * (flow_i128)b.z;
    flow_i128 x = (flow_i128)a.w * (flow_i128)b.x + (flow_i128)a.x * (flow_i128)b.w +
                  (flow_i128)a.y * (flow_i128)b.z - (flow_i128)a.z * (flow_i128)b.y;
    flow_i128 y = (flow_i128)a.w * (flow_i128)b.y - (flow_i128)a.x * (flow_i128)b.z +
                  (flow_i128)a.y * (flow_i128)b.w + (flow_i128)a.z * (flow_i128)b.x;
    flow_i128 z = (flow_i128)a.w * (flow_i128)b.z + (flow_i128)a.x * (flow_i128)b.y -
                  (flow_i128)a.y * (flow_i128)b.x + (flow_i128)a.z * (flow_i128)b.w;
    r.w = flow_sym_shr_i128(w, 16U);
    r.x = flow_sym_shr_i128(x, 16U);
    r.y = flow_sym_shr_i128(y, 16U);
    r.z = flow_sym_shr_i128(z, 16U);
    return robo_quat_normalize(r);
}

FLOW_INL robo_quat_t robo_quat_from_axis_angle(uint32_t axis, int64_t theta_q16) {
    int64_t s = 0, c = 0;
    robo_sincos_q16(flow_ashr_i64(theta_q16, 1U), &s, &c);
    robo_quat_t q;
    q.w = c;
    q.x = (axis == ROBO_ARM6_AXIS_X) ? s : 0;
    q.y = (axis == ROBO_ARM6_AXIS_Y) ? s : 0;
    q.z = (axis == ROBO_ARM6_AXIS_Z) ? s : 0;
    return q;
}

FLOW_INL robo_mat3_t robo_quat_to_mat3(robo_quat_t q) {
    int64_t xx = robo__mul_q16(q.x, q.x);
    int64_t yy = robo__mul_q16(q.y, q.y);
    int64_t zz = robo__mul_q16(q.z, q.z);
    int64_t xy = robo__mul_q16(q.x, q.y);
    int64_t xz = robo__mul_q16(q.x, q.z);
    int64_t yz = robo__mul_q16(q.y, q.z);
    int64_t wx = robo__mul_q16(q.w, q.x);
    int64_t wy = robo__mul_q16(q.w, q.y);
    int64_t wz = robo__mul_q16(q.w, q.z);

    robo_mat3_t R;
    R.m[0][0] = robo__sat_sub_i64(ROBO_Q16_ONE, flow_shl_i64(robo__sat_add_i64(yy, zz), 1U));
    R.m[0][1] = flow_shl_i64(robo__sat_sub_i64(xy, wz), 1U);
    R.m[0][2] = flow_shl_i64(robo__sat_add_i64(xz, wy), 1U);

    R.m[1][0] = flow_shl_i64(robo__sat_add_i64(xy, wz), 1U);
    R.m[1][1] = robo__sat_sub_i64(ROBO_Q16_ONE, flow_shl_i64(robo__sat_add_i64(xx, zz), 1U));
    R.m[1][2] = flow_shl_i64(robo__sat_sub_i64(yz, wx), 1U);

    R.m[2][0] = flow_shl_i64(robo__sat_sub_i64(xz, wy), 1U);
    R.m[2][1] = flow_shl_i64(robo__sat_add_i64(yz, wx), 1U);
    R.m[2][2] = robo__sat_sub_i64(ROBO_Q16_ONE, flow_shl_i64(robo__sat_add_i64(xx, yy), 1U));
    return R;
}

FLOW_FORCE_INL robo_vec3_t robo_mat3_mul_v3(const robo_mat3_t *R, robo_vec3_t v) {
    flow_i128 x = (flow_i128)R->m[0][0] * (flow_i128)v.x + (flow_i128)R->m[0][1] * (flow_i128)v.y +
                  (flow_i128)R->m[0][2] * (flow_i128)v.z;
    flow_i128 y = (flow_i128)R->m[1][0] * (flow_i128)v.x + (flow_i128)R->m[1][1] * (flow_i128)v.y +
                  (flow_i128)R->m[1][2] * (flow_i128)v.z;
    flow_i128 z = (flow_i128)R->m[2][0] * (flow_i128)v.x + (flow_i128)R->m[2][1] * (flow_i128)v.y +
                  (flow_i128)R->m[2][2] * (flow_i128)v.z;
    return robo_v3_make(flow_sym_shr_i128(x, 16U), flow_sym_shr_i128(y, 16U),
                        flow_sym_shr_i128(z, 16U));
}

/* --- Default 6-DoF Spatial Cobot Configuration (Pure-GPR integer Q16.16) --- */

FLOW_INL void robo_arm6_init_default(robo_arm6_state_t *st) {
    uint32_t i;
    if (!st) {
        return;
    }
    for (i = 0U; i < ROBO_ARM6_DOF; ++i) {
        st->q_q16[i]  = 0;
        st->dq_q16[i] = 0;
    }

    /* Link 0: Base Yaw (Z) at (0, 0, 0.16) */
    st->links[0].offset_q16      = robo_v3_make(0, 0, INT64_C(10486)); /* 0.16 m */
    st->links[0].axis            = ROBO_ARM6_AXIS_Z;
    st->links[0].cap_a_local_q16 = robo_v3_make(0, 0, INT64_C(-7864)); /* -0.12 m */
    st->links[0].cap_b_local_q16 = robo_v3_make(0, 0, INT64_C(5243));  /* +0.08 m */
    st->links[0].cap_radius_q16  = INT64_C(4260);                      /* 0.065 m */
    st->links[0].mass_q16        = INT64_C(229376);                    /* 3.5 kg */
    st->links[0].inertia_q16     = INT64_C(1966);                      /* 0.03 kg*m^2 */
    st->links[0].com_local_q16   = robo_v3_make(0, 0, INT64_C(-2621));

    /* Link 1: Shoulder Pitch (Y), upper arm 0.42 m along +Z */
    st->links[1].offset_q16      = robo_v3_make(0, 0, 0);
    st->links[1].axis            = ROBO_ARM6_AXIS_Y;
    st->links[1].cap_a_local_q16 = robo_v3_make(0, 0, INT64_C(2621));  /* 0.04 m */
    st->links[1].cap_b_local_q16 = robo_v3_make(0, 0, INT64_C(27525)); /* 0.42 m */
    st->links[1].cap_radius_q16  = INT64_C(3604);                      /* 0.055 m */
    st->links[1].mass_q16        = INT64_C(196608);                    /* 3.0 kg */
    st->links[1].inertia_q16     = INT64_C(2949);                      /* 0.045 */
    st->links[1].com_local_q16   = robo_v3_make(0, 0, INT64_C(13763)); /* 0.21 m */

    /* Link 2: Elbow Pitch (Y), forearm 0.39 m along +X */
    st->links[2].offset_q16      = robo_v3_make(0, 0, INT64_C(27525)); /* 0.42 m */
    st->links[2].axis            = ROBO_ARM6_AXIS_Y;
    st->links[2].cap_a_local_q16 = robo_v3_make(INT64_C(1966), 0, 0);  /* 0.03 m */
    st->links[2].cap_b_local_q16 = robo_v3_make(INT64_C(25559), 0, 0); /* 0.39 m */
    st->links[2].cap_radius_q16  = INT64_C(3146);                      /* 0.048 m */
    st->links[2].mass_q16        = INT64_C(144179);                    /* 2.2 kg */
    st->links[2].inertia_q16     = INT64_C(1966);                      /* 0.030 */
    st->links[2].com_local_q16   = robo_v3_make(INT64_C(12780), 0, 0); /* 0.195 m */

    /* Link 3: Wrist Roll (X), 0.12 m along +X */
    st->links[3].offset_q16      = robo_v3_make(INT64_C(25559), 0, 0); /* 0.39 m */
    st->links[3].axis            = ROBO_ARM6_AXIS_X;
    st->links[3].cap_a_local_q16 = robo_v3_make(0, 0, 0);
    st->links[3].cap_b_local_q16 = robo_v3_make(INT64_C(7864), 0, 0);  /* 0.12 m */
    st->links[3].cap_radius_q16  = INT64_C(2621);                      /* 0.040 m */
    st->links[3].mass_q16        = INT64_C(78643);                     /* 1.2 kg */
    st->links[3].inertia_q16     = INT64_C(655);                       /* 0.010 */
    st->links[3].com_local_q16   = robo_v3_make(INT64_C(3932), 0, 0);  /* 0.06 m */

    /* Link 4: Wrist Pitch (Y), 0.11 m along -Z */
    st->links[4].offset_q16      = robo_v3_make(INT64_C(7864), 0, 0);  /* 0.12 m */
    st->links[4].axis            = ROBO_ARM6_AXIS_Y;
    st->links[4].cap_a_local_q16 = robo_v3_make(0, 0, 0);
    st->links[4].cap_b_local_q16 = robo_v3_make(0, 0, INT64_C(-7209)); /* -0.11 m */
    st->links[4].cap_radius_q16  = INT64_C(2359);                      /* 0.036 m */
    st->links[4].mass_q16        = INT64_C(58982);                     /* 0.9 kg */
    st->links[4].inertia_q16     = INT64_C(459);                       /* 0.007 */
    st->links[4].com_local_q16   = robo_v3_make(0, 0, INT64_C(-3604)); /* -0.055 m */

    /* Link 5: Wrist Yaw / Flange (Z), 0.09 m along -Z */
    st->links[5].offset_q16      = robo_v3_make(0, 0, INT64_C(-7209)); /* -0.11 m */
    st->links[5].axis            = ROBO_ARM6_AXIS_Z;
    st->links[5].cap_a_local_q16 = robo_v3_make(0, 0, 0);
    st->links[5].cap_b_local_q16 = robo_v3_make(0, 0, INT64_C(-5898)); /* -0.09 m */
    st->links[5].cap_radius_q16  = INT64_C(2294);                      /* 0.035 m */
    st->links[5].mass_q16        = INT64_C(45875);                     /* 0.7 kg */
    st->links[5].inertia_q16     = INT64_C(328);                       /* 0.005 */
    st->links[5].com_local_q16   = robo_v3_make(0, 0, INT64_C(-2949)); /* -0.045 m */

    st->ee_offset_q16    = robo_v3_make(0, 0, INT64_C(-5898)); /* -0.09 m */
    st->ee_radius_q16    = INT64_C(2621);                      /* 0.040 m */
    st->payload_mass_q16 = INT64_C(65536);                     /* 1.0 kg */
    st->grav_q16         = INT64_C(642908);                    /* 9.81 m/s^2 */

    /* Joint limits (rad in Q16.16) */
    st->q_min_q16[0] = INT64_C(-183501); st->q_max_q16[0] = INT64_C(183501); /* +-2.8 rad */
    st->q_min_q16[1] = INT64_C(-117965); st->q_max_q16[1] = INT64_C(117965); /* +-1.8 rad */
    st->q_min_q16[2] = INT64_C(-131072); st->q_max_q16[2] = INT64_C(131072); /* +-2.0 rad */
    st->q_min_q16[3] = INT64_C(-183501); st->q_max_q16[3] = INT64_C(183501); /* +-2.8 rad */
    st->q_min_q16[4] = INT64_C(-137626); st->q_max_q16[4] = INT64_C(137626); /* +-2.1 rad */
    st->q_min_q16[5] = INT64_C(-196608); st->q_max_q16[5] = INT64_C(196608); /* +-3.0 rad */

    for (i = 0U; i < ROBO_ARM6_DOF; ++i) {
        st->dq_max_q16[i]  = INT64_C(196608);  /* 3.0 rad/s */
        st->ddq_max_q16[i] = INT64_C(1638400); /* 25.0 rad/s^2 */
    }
    st->tau_max_q16[0] = INT64_C(7864320); /* 120 Nm */
    st->tau_max_q16[1] = INT64_C(7864320); /* 120 Nm */
    st->tau_max_q16[2] = INT64_C(5242880); /* 80 Nm */
    st->tau_max_q16[3] = INT64_C(2621440); /* 40 Nm */
    st->tau_max_q16[4] = INT64_C(2621440); /* 40 Nm */
    st->tau_max_q16[5] = INT64_C(1966080); /* 30 Nm */

    st->table_z_q16      = 0;
    st->contact_k_q16    = INT64_C(131072000); /* 2000 N/m */
    st->contact_b_q16    = INT64_C(2621440);   /* 40 N*s/m */
    st->contact_fmax_q16 = INT64_C(1638400);   /* 25 N */

    st->tau_delay_q16 = INT64_C(2621);   /* 0.040 s (40 ms) */
    st->dt_q16        = INT64_C(66);     /* 0.001 s (1 kHz) */
    st->dist_pos_q16  = INT64_C(655);    /* 0.010 m/rad */
    st->dist_vel_q16  = INT64_C(3932);   /* 0.060 m/s */
    st->dist_acc_q16  = INT64_C(229376); /* 3.5 m/s^2 (+/-20% payload + discretization) */
    st->obs_amax_q16  = INT64_C(32768);  /* 0.5 m/s^2 */

    st->k1_q16 = INT64_C(1179648); /* 18.0 */
    st->k0_q16 = INT64_C(5308416); /* 81.0 */
    st->kv_q16 = INT64_C(1048576); /* 16.0 */

    st->n_obs = 0U;
    for (i = 0U; i < ROBO_ARM6_MAX_OBS; ++i) {
        st->obs[i].ox_q16     = 0;
        st->obs[i].oy_q16     = 0;
        st->obs[i].oz_q16     = 0;
        st->obs[i].radius_q16 = 0;
        st->obs[i].ovx_q16    = 0;
        st->obs[i].ovy_q16    = 0;
        st->obs[i].ovz_q16    = 0;
    }
}

/* --- 6-DoF Forward Kinematics, Link Capsules & Analytical Spatial Jacobian --- */

FLOW_INL void robo_arm6_point_jacobian(const robo_arm6_fk_t *fk, uint32_t link_idx,
                                       robo_vec3_t world_pt, int64_t Jv[3][ROBO_ARM6_DOF]) {
    uint32_t k;
    uint32_t lim = (link_idx < ROBO_ARM6_DOF) ? link_idx : (ROBO_ARM6_DOF - 1U);
    for (k = 0U; k < ROBO_ARM6_DOF; ++k) {
        if (k <= lim) {
            robo_vec3_t arm = robo_v3_sub(world_pt, fk->joint_pos[k]);
            robo_vec3_t col = robo_v3_cross(fk->joint_axis[k], arm);
            Jv[0][k]        = col.x;
            Jv[1][k]        = col.y;
            Jv[2][k]        = col.z;
        } else {
            Jv[0][k] = 0;
            Jv[1][k] = 0;
            Jv[2][k] = 0;
        }
    }
}

FLOW_INL robo_vec3_t robo_arm6_point_vel(const robo_arm6_fk_t *fk, const int64_t dq[ROBO_ARM6_DOF],
                                         uint32_t link_idx, robo_vec3_t world_pt) {
    int64_t Jv[3][ROBO_ARM6_DOF];
    robo_arm6_point_jacobian(fk, link_idx, world_pt, Jv);
    return robo_v3_make(robo__dot_nd_q16(Jv[0], dq, ROBO_ARM6_DOF),
                        robo__dot_nd_q16(Jv[1], dq, ROBO_ARM6_DOF),
                        robo__dot_nd_q16(Jv[2], dq, ROBO_ARM6_DOF));
}

/* Exact spatial bias acceleration Jdot(q, dq)*dq for world_pt rigidly attached to link_idx */
FLOW_INL robo_vec3_t robo_arm6_point_jdot_dq(const robo_arm6_fk_t *fk,
                                             const int64_t dq[ROBO_ARM6_DOF], uint32_t link_idx,
                                             robo_vec3_t world_pt) {
    robo_vec3_t p_vel = robo_arm6_point_vel(fk, dq, link_idx, world_pt);
    robo_vec3_t acc   = robo_v3_make(0, 0, 0);
    uint32_t lim      = (link_idx < ROBO_ARM6_DOF) ? link_idx : (ROBO_ARM6_DOF - 1U);
    uint32_t k;

    for (k = 0U; k <= lim; ++k) {
        robo_vec3_t w_prev = (k == 0U) ? robo_v3_make(0, 0, 0) : fk->omega[k - 1U];
        robo_vec3_t zdot   = robo_v3_cross(w_prev, fk->joint_axis[k]);
        robo_vec3_t arm    = robo_v3_sub(world_pt, fk->joint_pos[k]);
        robo_vec3_t vdiff  = robo_v3_sub(p_vel, fk->joint_vel[k]);
        robo_vec3_t term   = robo_v3_add(robo_v3_cross(zdot, arm),
                                         robo_v3_cross(fk->joint_axis[k], vdiff));
        acc                = robo_v3_add(acc, robo_v3_scale(term, dq[k]));
    }
    return acc;
}

FLOW_INL void robo_arm6_fk_ring0(const robo_arm6_state_t *st, robo_arm6_fk_t *fk) {
    if (!st || !fk) {
        return;
    }
    robo_vec3_t p_cur = robo_v3_make(0, 0, 0);
    robo_quat_t q_cur = robo_quat_identity();
    robo_mat3_t R_cur = robo_quat_to_mat3(q_cur);
    robo_vec3_t w_cur = robo_v3_make(0, 0, 0);
    uint32_t k, r;

    for (k = 0U; k < ROBO_ARM6_DOF; ++k) {
        robo_vec3_t step = robo_mat3_mul_v3(&R_cur, st->links[k].offset_q16);
        p_cur            = robo_v3_add(p_cur, step);
        fk->joint_pos[k] = p_cur;

        robo_vec3_t local_ax = robo_v3_make(
            (st->links[k].axis == ROBO_ARM6_AXIS_X) ? ROBO_Q16_ONE : 0,
            (st->links[k].axis == ROBO_ARM6_AXIS_Y) ? ROBO_Q16_ONE : 0,
            (st->links[k].axis == ROBO_ARM6_AXIS_Z) ? ROBO_Q16_ONE : 0);
        robo_vec3_t z_world = robo_mat3_mul_v3(&R_cur, local_ax);
        fk->joint_axis[k]   = z_world;

        robo_quat_t q_j   = robo_quat_from_axis_angle(st->links[k].axis, st->q_q16[k]);
        q_cur             = robo_quat_mul(q_cur, q_j);
        R_cur             = robo_quat_to_mat3(q_cur);
        fk->link_quat[k]  = q_cur;
        fk->link_rot[k]   = R_cur;

        w_cur         = robo_v3_add(w_cur, robo_v3_scale(z_world, st->dq_q16[k]));
        fk->omega[k]  = w_cur;

        /* Velocity of joint_pos[k] due to preceding joints 0..k-1 */
        robo_vec3_t jv = robo_v3_make(0, 0, 0);
        for (r = 0U; r < k; ++r) {
            robo_vec3_t arm = robo_v3_sub(p_cur, fk->joint_pos[r]);
            jv              = robo_v3_add(jv, robo_v3_scale(robo_v3_cross(fk->joint_axis[r], arm),
                                                            st->dq_q16[r]));
        }
        fk->joint_vel[k] = jv;

        fk->cap_a[k] = robo_v3_add(p_cur, robo_mat3_mul_v3(&R_cur, st->links[k].cap_a_local_q16));
        fk->cap_b[k] = robo_v3_add(p_cur, robo_mat3_mul_v3(&R_cur, st->links[k].cap_b_local_q16));
        fk->com_world[k] = robo_v3_add(p_cur, robo_mat3_mul_v3(&R_cur, st->links[k].com_local_q16));
    }

    fk->ee_quat = q_cur;
    fk->ee_pos  = robo_v3_add(p_cur, robo_mat3_mul_v3(&R_cur, st->ee_offset_q16));
    robo_arm6_point_jacobian(fk, ROBO_ARM6_DOF - 1U, fk->ee_pos, fk->J_ee);
    for (k = 0U; k < ROBO_ARM6_DOF; ++k) {
        fk->J_rot[0][k] = fk->joint_axis[k].x;
        fk->J_rot[1][k] = fk->joint_axis[k].y;
        fk->J_rot[2][k] = fk->joint_axis[k].z;
    }
    fk->ee_vel = robo_v3_make(robo__dot_nd_q16(fk->J_ee[0], st->dq_q16, ROBO_ARM6_DOF),
                              robo__dot_nd_q16(fk->J_ee[1], st->dq_q16, ROBO_ARM6_DOF),
                              robo__dot_nd_q16(fk->J_ee[2], st->dq_q16, ROBO_ARM6_DOF));
}

/* --- Exact 3D Capsule-to-Sphere & Capsule-to-Capsule Distances & Normals --- */

FLOW_INL int64_t robo_capsule_sphere_dist_ring0(robo_vec3_t a, robo_vec3_t b, int64_t r_cap,
                                                robo_vec3_t center, int64_t r_sph,
                                                robo_vec3_t *cp_out, robo_vec3_t *normal_out) {
    robo_vec3_t d   = robo_v3_sub(b, a);
    flow_u128   l2  = robo_v3_norm2_q32(d);
    int64_t     l2q = flow_sym_shr_i128((flow_i128)l2, 16U);
    int64_t     t   = 0;

    if (l2q > INT64_C(16)) {
        int64_t ca_dot_d = robo_v3_dot(robo_v3_sub(center, a), d);
        t = flow_sym_div128_clamped(flow_shl_i128((flow_i128)ca_dot_d, 16U), (flow_i128)l2q,
                                    ROBO_Q16_ONE);
        t = robo__clamp_i64(t, 0, ROBO_Q16_ONE);
    }

    robo_vec3_t cp = robo_v3_add(a, robo_v3_scale(d, t));
    robo_vec3_t v  = robo_v3_sub(cp, center);
    int64_t     dist = 0;
    robo_vec3_t n    = robo_v3_unit(v, &dist);

    if (cp_out) {
        *cp_out = cp;
    }
    if (normal_out) {
        *normal_out = n;
    }
    return robo__sat_sub_i64(dist, robo__sat_add_i64(r_cap, r_sph));
}

FLOW_INL int64_t robo_capsule_capsule_dist_ring0(robo_vec3_t a1, robo_vec3_t b1, int64_t r1,
                                                 robo_vec3_t a2, robo_vec3_t b2, int64_t r2,
                                                 robo_vec3_t *p1_out, robo_vec3_t *p2_out,
                                                 robo_vec3_t *normal_out) {
    robo_vec3_t d1 = robo_v3_sub(b1, a1);
    robo_vec3_t d2 = robo_v3_sub(b2, a2);
    robo_vec3_t r  = robo_v3_sub(a1, a2);

    int64_t a = flow_sym_shr_i128((flow_i128)robo_v3_norm2_q32(d1), 16U);
    int64_t e = flow_sym_shr_i128((flow_i128)robo_v3_norm2_q32(d2), 16U);
    int64_t f = robo_v3_dot(d2, r);
    int64_t s = 0;
    int64_t t = 0;

    if (a <= INT64_C(16) && e <= INT64_C(16)) {
        s = 0;
        t = 0;
    } else if (a <= INT64_C(16)) {
        s = 0;
        t = robo__clamp_i64(
            flow_sym_div128_clamped(flow_shl_i128((flow_i128)robo__sat_neg_i64(f), 16U),
                                    (flow_i128)e, ROBO_Q16_ONE),
            0, ROBO_Q16_ONE);
    } else {
        int64_t c = robo_v3_dot(d1, r);
        if (e <= INT64_C(16)) {
            t = 0;
            s = robo__clamp_i64(
                flow_sym_div128_clamped(flow_shl_i128((flow_i128)robo__sat_neg_i64(c), 16U),
                                        (flow_i128)a, ROBO_Q16_ONE),
                0, ROBO_Q16_ONE);
        } else {
            int64_t   b     = robo_v3_dot(d1, d2);
            flow_i128 denom = (flow_i128)a * (flow_i128)e - (flow_i128)b * (flow_i128)b;
            if (denom > (flow_i128)INT64_C(65536)) {
                flow_i128 s_num = (flow_i128)b * (flow_i128)f - (flow_i128)c * (flow_i128)e;
                s = robo__clamp_i64(flow_sym_div128_clamped(flow_shl_i128(s_num, 16U), denom,
                                                            ROBO_Q16_ONE),
                                    0, ROBO_Q16_ONE);
            } else {
                s = 0;
            }
            flow_i128 t_num = (flow_i128)b * (flow_i128)s + flow_shl_i128((flow_i128)f, 16U);
            if (t_num < 0) {
                t = 0;
                s = robo__clamp_i64(
                    flow_sym_div128_clamped(flow_shl_i128((flow_i128)robo__sat_neg_i64(c), 16U),
                                            (flow_i128)a, ROBO_Q16_ONE),
                    0, ROBO_Q16_ONE);
            } else {
                flow_i128 e_q32 = flow_shl_i128((flow_i128)e, 16U);
                if (t_num > e_q32) {
                    t = ROBO_Q16_ONE;
                    s = robo__clamp_i64(
                        flow_sym_div128_clamped(
                            flow_shl_i128((flow_i128)robo__sat_sub_i64(b, c), 16U),
                            (flow_i128)a, ROBO_Q16_ONE),
                        0, ROBO_Q16_ONE);
                } else {
                    t = flow_sym_div128_clamped(t_num, (flow_i128)e, ROBO_Q16_ONE);
                }
            }
        }
    }

    robo_vec3_t p1   = robo_v3_add(a1, robo_v3_scale(d1, s));
    robo_vec3_t p2   = robo_v3_add(a2, robo_v3_scale(d2, t));
    robo_vec3_t v    = robo_v3_sub(p1, p2);
    int64_t     dist = 0;
    robo_vec3_t n    = robo_v3_unit(v, &dist);

    if (p1_out) {
        *p1_out = p1;
    }
    if (p2_out) {
        *p2_out = p2;
    }
    if (normal_out) {
        *normal_out = n;
    }
    return robo__sat_sub_i64(dist, robo__sat_add_i64(r1, r2));
}

/* --- 6-DoF Rigid-Body Mass Matrix M(q), Gravity g(q) & SDD Certification --- */

FLOW_INL void robo_arm6_mass_gravity_ring0(const robo_arm6_state_t *st, const robo_arm6_fk_t *fk,
                                           int64_t M[ROBO_ARM6_DOF][ROBO_ARM6_DOF],
                                           int64_t g[ROBO_ARM6_DOF], uint8_t *sdd_ok_out,
                                           int64_t *sdd_mu_out) {
    uint32_t i, j, k, c;
    int64_t Jv[3][ROBO_ARM6_DOF];
    for (i = 0U; i < ROBO_ARM6_DOF; ++i) {
        g[i] = 0;
        for (j = 0U; j < ROBO_ARM6_DOF; ++j) {
            M[i][j] = (i == j) ? INT64_C(6554) : 0; /* 0.10 kg*m^2 rotor inertia */
        }
    }

    for (k = 0U; k < ROBO_ARM6_DOF; ++k) {
        int64_t mk = st->links[k].mass_q16;
        int64_t ik = st->links[k].inertia_q16;
        robo_arm6_point_jacobian(fk, k, fk->com_world[k], Jv);
        for (i = 0U; i < ROBO_ARM6_DOF; ++i) {
            g[i] = robo__sat_add_i64(
                g[i], robo__mul_q16(mk, robo__mul_q16(st->grav_q16, Jv[2][i])));
            for (j = i; j < ROBO_ARM6_DOF; ++j) {
                flow_i128 jtj = 0;
                for (c = 0U; c < 3U; ++c) {
                    jtj += (flow_i128)Jv[c][i] * (flow_i128)Jv[c][j];
                }
                int64_t add_m = robo__mul_q16(mk, flow_sym_shr_i128(jtj, 16U));
                if (i <= k && j <= k) {
                    int64_t rdot = robo_v3_dot(fk->joint_axis[i], fk->joint_axis[j]);
                    add_m        = robo__sat_add_i64(add_m, robo__mul_q16(ik, rdot));
                }
                M[i][j] = robo__sat_add_i64(M[i][j], add_m);
                M[j][i] = M[i][j];
            }
        }
    }

    /* Payload contribution at end-effector */
    int64_t mp = st->payload_mass_q16;
    if (mp > 0) {
        for (i = 0U; i < ROBO_ARM6_DOF; ++i) {
            g[i] = robo__sat_add_i64(
                g[i], robo__mul_q16(mp, robo__mul_q16(st->grav_q16, fk->J_ee[2][i])));
            for (j = i; j < ROBO_ARM6_DOF; ++j) {
                flow_i128 jtj = 0;
                for (c = 0U; c < 3U; ++c) {
                    jtj += (flow_i128)fk->J_ee[c][i] * (flow_i128)fk->J_ee[c][j];
                }
                int64_t add_m = robo__mul_q16(mp, flow_sym_shr_i128(jtj, 16U));
                M[i][j]       = robo__sat_add_i64(M[i][j], add_m);
                M[j][i]       = M[i][j];
            }
        }
    }

    /* Certify SPD of M(q) via Cholesky + FLOW flow_cholesky_sdd_floor_i64 */
    if (sdd_ok_out || sdd_mu_out) {
        int64_t L[ROBO_ARM6_DOF][ROBO_ARM6_DOF];
        int64_t M_q32[ROBO_ARM6_DOF * ROBO_ARM6_DOF];
        int64_t L_q16[ROBO_ARM6_DOF * ROBO_ARM6_DOF];
        for (i = 0U; i < ROBO_ARM6_DOF; ++i) {
            for (j = 0U; j < ROBO_ARM6_DOF; ++j) {
                L[i][j]                        = 0;
                M_q32[i * ROBO_ARM6_DOF + j]   = (int64_t)flow_shl_i128((flow_i128)M[i][j], 16U);
            }
        }
        for (j = 0U; j < ROBO_ARM6_DOF; ++j) {
            flow_i128 ssq = 0;
            for (k = 0U; k < j; ++k) {
                ssq += (flow_i128)L[j][k] * (flow_i128)L[j][k];
            }
            int64_t d = robo__sat_sub_i64(M[j][j], flow_sym_shr_i128(ssq, 16U));
            if (d < INT64_C(256)) d = INT64_C(256);
            L[j][j] = robo__sqrt_q16((uint64_t)flow_shl_i128((flow_i128)d, 16U));
            for (i = j + 1U; i < ROBO_ARM6_DOF; ++i) {
                flow_i128 s = 0;
                for (k = 0U; k < j; ++k) {
                    s += (flow_i128)L[i][k] * (flow_i128)L[j][k];
                }
                int64_t diff = robo__sat_sub_i64(M[i][j], flow_sym_shr_i128(s, 16U));
                L[i][j]      = flow_sym_div128_clamped(flow_shl_i128((flow_i128)diff, 16U),
                                                       (flow_i128)L[j][j], FLOW_MAX_Q16);
            }
        }
        for (i = 0U; i < ROBO_ARM6_DOF; ++i) {
            for (j = 0U; j < ROBO_ARM6_DOF; ++j) {
                L_q16[i * ROBO_ARM6_DOF + j] = L[i][j];
            }
        }
        int64_t mu_q32 = 0;
        int rc = flow_cholesky_sdd_floor_i64((int)ROBO_ARM6_DOF, (int)ROBO_ARM6_DOF, M_q32, L_q16,
                                             &mu_q32);
        if (sdd_mu_out) *sdd_mu_out = mu_q32;
        if (sdd_ok_out) *sdd_ok_out = (rc == 0 && mu_q32 >= -(ROBO_Q16_ONE * 16)) ? 1U : 0U;
    }
}

/* --- 6-DoF Spatial Tube-CBF Safety Filter (Acceleration Input) --- */

typedef struct {
    int64_t h_raw_q16;
    int64_t h_tube_q16;
    int64_t hdot_q16;
    int64_t tube_pad_q16;
    int64_t gq_q16[ROBO_ARM6_DOF];
    int64_t rhs_bias_q16;
} robo__arm6_cand_con_t;

FLOW_INL void robo__insert_top_con(robo__arm6_cand_con_t *arr, uint32_t *count, uint32_t cap,
                                   const robo__arm6_cand_con_t *cand) {
    uint32_t n = *count;
    uint32_t pos = n;
    uint32_t k;
    for (k = 0U; k < n; ++k) {
        if (cand->h_tube_q16 < arr[k].h_tube_q16) {
            pos = k;
            break;
        }
    }
    if (pos >= cap) {
        return;
    }
    uint32_t limit = (n < cap) ? n : (cap - 1U);
    for (k = limit; k > pos; --k) {
        arr[k] = arr[k - 1U];
    }
    arr[pos] = *cand;
    if (n < cap) {
        *count = n + 1U;
    }
}

FLOW_INL void robo_arm6_cbf_accel_ring0(const robo_arm6_state_t *st,
                                        const int64_t ddq_nom_q16[ROBO_ARM6_DOF],
                                        robo_arm6_cert_t *cert) {
    uint32_t i, j;
    if (!cert) {
        return;
    }
    for (j = 0U; j < ROBO_ARM6_DOF; ++j) {
        cert->ddq_safe_q16[j] = 0;
        cert->tau_safe_q16[j] = 0;
    }
    cert->verdict            = ROBO_VERDICT_VETOED_INFEASIBLE;
    cert->n_constraints      = 0U;
    cert->qp_iters           = 0U;
    cert->feasible           = 0U;
    cert->slack_active       = 0U;
    cert->fallback_engaged   = 1U;
    cert->mass_sdd_certified = 0U;
    cert->ee_pos             = robo_v3_make(0, 0, 0);
    cert->ee_quat            = robo_quat_identity();
    cert->min_obs_h_q16      = FLOW_MAX_Q16;
    cert->min_self_h_q16     = FLOW_MAX_Q16;
    cert->min_contact_h_q16  = FLOW_MAX_Q16;
    cert->min_joint_h_q16    = FLOW_MAX_Q16;
    cert->min_stl_h_q16      = FLOW_MAX_Q16;
    cert->tube_margin_q16    = 0;
    cert->max_slack_q16      = 0;
    cert->mass_sdd_mu_q32    = 0;

    if (!st || !ddq_nom_q16) {
        return;
    }

    robo_arm6_fk_t fk;
    robo_arm6_fk_ring0(st, &fk);
    cert->ee_pos  = fk.ee_pos;
    cert->ee_quat = fk.ee_quat;

    int64_t M[ROBO_ARM6_DOF][ROBO_ARM6_DOF];
    int64_t grav[ROBO_ARM6_DOF];
    robo_arm6_mass_gravity_ring0(st, &fk, M, grav, &cert->mass_sdd_certified,
                                 &cert->mass_sdd_mu_q32);

    /* Explicit VLA action-delay Tube-CBF parameters */
    int64_t tau_eff  = robo__sat_add_i64(st->tau_delay_q16, st->dt_q16);
    int64_t tau_sq_h = flow_ashr_i64(robo__mul_q16(tau_eff, tau_eff), 1U);
    int64_t a_pert   = robo__sat_add_i64(st->dist_acc_q16,
                                         robo__sat_add_i64(st->obs_amax_q16, INT64_C(262144)));
    int64_t base_pad = robo__sat_add_i64(st->dist_pos_q16, robo__mul_q16(a_pert, tau_sq_h));

    robo_qp_nd_problem_t prob;
    prob.n         = ROBO_ARM6_DOF;
    prob.m         = 0U;
    prob.max_iters = 64U;
    for (j = 0U; j < ROBO_QP_MAX_N; ++j) {
        prob.u_nom_q16[j]   = 0;
        prob.u_lo_q16[j]    = 0;
        prob.u_hi_q16[j]    = 0;
        prob.u_brake_q16[j] = 0;
        prob.w_diag_q16[j]  = ROBO_Q16_ONE;
    }
    for (i = 0U; i < ROBO_QP_MAX_M; ++i) {
        for (j = 0U; j < ROBO_QP_MAX_N; ++j) {
            prob.con[i].a_q16[j] = 0;
        }
        prob.con[i].b_q16          = 0;
        prob.con[i].slack_max_q16  = 0;
        prob.con[i].slack_lin_q16  = 0;
        prob.con[i].slack_quad_q16 = 0;
        prob.con[i].is_hard        = 0U;
    }

    /* 1. Joint Position + Velocity + Acceleration Box Bounds (exact 1D coordinate bounds) */
    int64_t min_j_h = FLOW_MAX_Q16;
    for (j = 0U; j < ROBO_ARM6_DOF; ++j) {
        prob.u_nom_q16[j]  = ddq_nom_q16[j];
        prob.w_diag_q16[j] = ROBO_Q16_ONE;

        int64_t dq_abs   = (int64_t)flow_abs_i64(st->dq_q16[j]);
        int64_t j_pad    = robo__sat_add_i64(st->dist_pos_q16,
                                             robo__sat_add_i64(robo__mul_q16(dq_abs, tau_eff),
                                                               robo__mul_q16(st->ddq_max_q16[j],
                                                                             tau_sq_h)));
        int64_t h_jlo    = robo__sat_sub_i64(st->q_q16[j], st->q_min_q16[j]);
        int64_t h_jhi    = robo__sat_sub_i64(st->q_max_q16[j], st->q_q16[j]);
        int64_t h_vlo    = robo__sat_add_i64(st->dq_max_q16[j], st->dq_q16[j]);
        int64_t h_vhi    = robo__sat_sub_i64(st->dq_max_q16[j], st->dq_q16[j]);
        min_j_h          = robo__min_i64(min_j_h, robo__min_i64(h_jlo, h_jhi));
        min_j_h          = robo__min_i64(min_j_h, robo__min_i64(h_vlo, h_vhi));

        int64_t h_jlo_t  = robo__sat_sub_i64(h_jlo, j_pad);
        int64_t h_jhi_t  = robo__sat_sub_i64(h_jhi, j_pad);
        int64_t dq_lo_w  = robo__sat_sub_i64(st->dq_q16[j], st->dist_vel_q16);
        int64_t dq_hi_w  = robo__sat_sub_i64(robo__sat_neg_i64(st->dq_q16[j]), st->dist_vel_q16);

        int64_t lo_pos = robo__sat_add_i64(
            st->dist_acc_q16,
            robo__sat_neg_i64(robo__sat_add_i64(robo__mul_q16(st->k1_q16, dq_lo_w),
                                                robo__mul_q16(st->k0_q16, h_jlo_t))));
        int64_t hi_pos = robo__sat_sub_i64(
            robo__sat_add_i64(robo__mul_q16(st->k1_q16, dq_hi_w),
                              robo__mul_q16(st->k0_q16, h_jhi_t)),
            st->dist_acc_q16);

        int64_t vel_pad = robo__sat_add_i64(st->dist_vel_q16,
                                            robo__mul_q16(st->ddq_max_q16[j], tau_eff));
        int64_t lo_vel  = robo__sat_add_i64(
            st->dist_acc_q16,
            robo__sat_neg_i64(robo__mul_q16(st->kv_q16, robo__sat_sub_i64(h_vlo, vel_pad))));
        int64_t hi_vel  = robo__sat_sub_i64(
            robo__mul_q16(st->kv_q16, robo__sat_sub_i64(h_vhi, vel_pad)), st->dist_acc_q16);

        int64_t neg_ddq_max = robo__sat_neg_i64(st->ddq_max_q16[j]);
        int64_t lo = robo__max_i64(neg_ddq_max, robo__max_i64(lo_pos, lo_vel));
        int64_t hi = robo__min_i64(st->ddq_max_q16[j], robo__min_i64(hi_pos, hi_vel));

        int64_t brake_j = robo__clamp_i64(
            robo__sat_neg_i64(flow_shl_i64(st->dq_q16[j], 4U)), neg_ddq_max, st->ddq_max_q16[j]);
        if (lo > hi) {
            /* If near a limit at high speed, bias box to maximum braking away from limit */
            int64_t clamped_b = (h_jlo < h_jhi) ? st->ddq_max_q16[j] : neg_ddq_max;
            lo = clamped_b;
            hi = clamped_b;
        }
        prob.u_lo_q16[j]    = lo;
        prob.u_hi_q16[j]    = hi;
        prob.u_brake_q16[j] = robo__clamp_i64(brake_j, lo, hi);
    }
    cert->min_joint_h_q16 = min_j_h;

    /* 2. 3D Capsule-to-Sphere Obstacle Tube-ECBFs (select top 6 closest across links 1..5 + EE) */
    robo__arm6_cand_con_t obs_cands[6];
    uint32_t n_obs_cands = 0U;
    int64_t  min_obs_h   = FLOW_MAX_Q16;
    int64_t  max_tube    = base_pad;
    uint32_t n_obs       = (st->n_obs > ROBO_ARM6_MAX_OBS) ? ROBO_ARM6_MAX_OBS : st->n_obs;
    uint32_t l_idx, o_idx;

    for (l_idx = 1U; l_idx < ROBO_ARM6_DOF; ++l_idx) {
        for (o_idx = 0U; o_idx < n_obs; ++o_idx) {
            robo_vec3_t obs_c = robo_v3_make(st->obs[o_idx].ox_q16, st->obs[o_idx].oy_q16,
                                             st->obs[o_idx].oz_q16);
            robo_vec3_t obs_v = robo_v3_make(st->obs[o_idx].ovx_q16, st->obs[o_idx].ovy_q16,
                                             st->obs[o_idx].ovz_q16);
            robo_vec3_t cp, n_out;
            int64_t h_raw = robo_capsule_sphere_dist_ring0(
                fk.cap_a[l_idx], fk.cap_b[l_idx], st->links[l_idx].cap_radius_q16, obs_c,
                st->obs[o_idx].radius_q16, &cp, &n_out);
            if (h_raw < min_obs_h) {
                min_obs_h = h_raw;
            }

            int64_t Jv[3][ROBO_ARM6_DOF];
            robo_arm6_point_jacobian(&fk, l_idx, cp, Jv);
            robo__arm6_cand_con_t cand;
            for (j = 0U; j < ROBO_ARM6_DOF; ++j) {
                flow_i128 g = (flow_i128)n_out.x * (flow_i128)Jv[0][j] +
                              (flow_i128)n_out.y * (flow_i128)Jv[1][j] +
                              (flow_i128)n_out.z * (flow_i128)Jv[2][j];
                cand.gq_q16[j] = flow_sym_shr_i128(g, 16U);
            }
            robo_vec3_t p_vel   = robo_arm6_point_vel(&fk, st->dq_q16, l_idx, cp);
            robo_vec3_t rel_vel = robo_v3_sub(p_vel, obs_v);
            int64_t     hdot    = robo_v3_dot(n_out, rel_vel);
            int64_t     v_app   = (hdot < 0) ? robo__sat_neg_i64(hdot) : 0;
            int64_t     t_pad   = robo__sat_add_i64(base_pad, robo__mul_q16(v_app, tau_eff));
            if (t_pad > max_tube) {
                max_tube = t_pad;
            }

            robo_vec3_t jdot_dq = robo_arm6_point_jdot_dq(&fk, st->dq_q16, l_idx, cp);
            int64_t     bias_a  = robo_v3_dot(n_out, jdot_dq);
            int64_t     h_tube  = robo__sat_sub_i64(h_raw, t_pad);
            int64_t     hd_w    = robo__sat_sub_i64(hdot, st->dist_vel_q16);
            int64_t     rhs     = robo__sat_sub_i64(
                robo__sat_sub_i64(
                    robo__sat_add_i64(st->dist_acc_q16, st->obs_amax_q16), bias_a),
                robo__sat_add_i64(robo__mul_q16(st->k1_q16, hd_w),
                                  robo__mul_q16(st->k0_q16, h_tube)));

            cand.h_raw_q16    = h_raw;
            cand.h_tube_q16   = h_tube;
            cand.hdot_q16     = hdot;
            cand.tube_pad_q16 = t_pad;
            cand.rhs_bias_q16 = rhs;
            robo__insert_top_con(obs_cands, &n_obs_cands, 6U, &cand);
        }
    }
    cert->min_obs_h_q16 = min_obs_h;

    for (i = 0U; i < n_obs_cands; ++i) {
        uint32_t m = prob.m;
        for (j = 0U; j < ROBO_ARM6_DOF; ++j) {
            prob.con[m].a_q16[j] = obs_cands[i].gq_q16[j];
        }
        prob.con[m].b_q16          = obs_cands[i].rhs_bias_q16;
        prob.con[m].slack_max_q16  = INT64_C(1310720);   /* 20.0 */
        prob.con[m].slack_lin_q16  = INT64_C(134217728); /* 2048.0 exact L1 barrier */
        prob.con[m].slack_quad_q16 = INT64_C(134217728); /* 2048.0 */
        prob.con[m].is_hard        = 0U;
        prob.m                     = m + 1U;
    }

    /* 3. 3D Capsule-to-Capsule Self-Collision Tube-ECBFs (non-adjacent pairs, top 4 closest) */
    const uint8_t self_pairs[9][2] = {
        {0U, 2U}, {0U, 3U}, {0U, 4U}, {0U, 5U}, {1U, 3U}, {1U, 4U}, {1U, 5U}, {2U, 4U}, {2U, 5U}
    };
    robo__arm6_cand_con_t self_cands[4];
    uint32_t n_self_cands = 0U;
    int64_t  min_self_h   = FLOW_MAX_Q16;

    for (i = 0U; i < 9U; ++i) {
        uint32_t li = (uint32_t)self_pairs[i][0];
        uint32_t lj = (uint32_t)self_pairs[i][1];
        robo_vec3_t p1, p2, n_out;
        int64_t h_raw = robo_capsule_capsule_dist_ring0(
            fk.cap_a[li], fk.cap_b[li], st->links[li].cap_radius_q16, fk.cap_a[lj], fk.cap_b[lj],
            st->links[lj].cap_radius_q16, &p1, &p2, &n_out);
        if (h_raw < min_self_h) {
            min_self_h = h_raw;
        }

        int64_t J1[3][ROBO_ARM6_DOF], J2[3][ROBO_ARM6_DOF];
        robo_arm6_point_jacobian(&fk, li, p1, J1);
        robo_arm6_point_jacobian(&fk, lj, p2, J2);

        robo__arm6_cand_con_t cand;
        for (j = 0U; j < ROBO_ARM6_DOF; ++j) {
            int64_t djx = robo__sat_sub_i64(J1[0][j], J2[0][j]);
            int64_t djy = robo__sat_sub_i64(J1[1][j], J2[1][j]);
            int64_t djz = robo__sat_sub_i64(J1[2][j], J2[2][j]);
            flow_i128 g = (flow_i128)n_out.x * (flow_i128)djx +
                          (flow_i128)n_out.y * (flow_i128)djy +
                          (flow_i128)n_out.z * (flow_i128)djz;
            cand.gq_q16[j] = flow_sym_shr_i128(g, 16U);
        }

        robo_vec3_t v1      = robo_arm6_point_vel(&fk, st->dq_q16, li, p1);
        robo_vec3_t v2      = robo_arm6_point_vel(&fk, st->dq_q16, lj, p2);
        int64_t     hdot    = robo_v3_dot(n_out, robo_v3_sub(v1, v2));
        int64_t     v_app   = (hdot < 0) ? robo__sat_neg_i64(hdot) : 0;
        int64_t     t_pad   = robo__sat_add_i64(base_pad, robo__mul_q16(v_app, tau_eff));
        if (t_pad > max_tube) {
            max_tube = t_pad;
        }

        robo_vec3_t jd1     = robo_arm6_point_jdot_dq(&fk, st->dq_q16, li, p1);
        robo_vec3_t jd2     = robo_arm6_point_jdot_dq(&fk, st->dq_q16, lj, p2);
        int64_t     bias_a  = robo_v3_dot(n_out, robo_v3_sub(jd1, jd2));
        int64_t     h_tube  = robo__sat_sub_i64(h_raw, t_pad);
        int64_t     hd_w    = robo__sat_sub_i64(hdot, st->dist_vel_q16);
        int64_t     rhs     = robo__sat_sub_i64(
            robo__sat_sub_i64(st->dist_acc_q16, bias_a),
            robo__sat_add_i64(robo__mul_q16(st->k1_q16, hd_w), robo__mul_q16(st->k0_q16, h_tube)));

        cand.h_raw_q16    = h_raw;
        cand.h_tube_q16   = h_tube;
        cand.hdot_q16     = hdot;
        cand.tube_pad_q16 = t_pad;
        cand.rhs_bias_q16 = rhs;
        robo__insert_top_con(self_cands, &n_self_cands, 4U, &cand);
    }
    cert->min_self_h_q16 = min_self_h;

    for (i = 0U; i < n_self_cands; ++i) {
        uint32_t m = prob.m;
        for (j = 0U; j < ROBO_ARM6_DOF; ++j) {
            prob.con[m].a_q16[j] = self_cands[i].gq_q16[j];
        }
        prob.con[m].b_q16          = self_cands[i].rhs_bias_q16;
        prob.con[m].slack_max_q16  = INT64_C(1310720);
        prob.con[m].slack_lin_q16  = INT64_C(134217728);
        prob.con[m].slack_quad_q16 = INT64_C(134217728);
        prob.con[m].is_hard        = 0U;
        prob.m                     = m + 1U;
    }

    /* 4. Table Contact-Force Envelope Tube-ECBF (EE + Wrist Link 4) */
    {
        int64_t pen_max = (st->contact_k_q16 > INT64_C(65536))
                              ? flow_sym_div128_clamped(
                                    flow_shl_i128((flow_i128)st->contact_fmax_q16, 16U),
                                    (flow_i128)st->contact_k_q16, FLOW_MAX_Q16)
                              : 0;
        /* Leave 15% safety margin on F_max for damping + mass uncertainty */
        int64_t pen_safe = robo__mul_q16(pen_max, INT64_C(52428)); /* 0.80 * pen_max */
        int64_t z_floor_ee = robo__sat_sub_i64(
            robo__sat_add_i64(st->table_z_q16, st->ee_radius_q16), pen_safe);
        int64_t h_ee = robo__sat_sub_i64(fk.ee_pos.z, z_floor_ee);
        int64_t v_app_ee = (fk.ee_vel.z < 0) ? robo__sat_neg_i64(fk.ee_vel.z) : 0;
        int64_t pad_ee   = robo__sat_add_i64(base_pad, robo__mul_q16(v_app_ee, tau_eff));
        if (pad_ee > max_tube) max_tube = pad_ee;
        int64_t h_ee_t   = robo__sat_sub_i64(h_ee, pad_ee);
        robo_vec3_t jd_ee = robo_arm6_point_jdot_dq(&fk, st->dq_q16, ROBO_ARM6_DOF - 1U, fk.ee_pos);
        int64_t hd_ee_w  = robo__sat_sub_i64(fk.ee_vel.z, st->dist_vel_q16);
        int64_t rhs_ee   = robo__sat_sub_i64(
            robo__sat_sub_i64(st->dist_acc_q16, jd_ee.z),
            robo__sat_add_i64(robo__mul_q16(st->k1_q16, hd_ee_w),
                              robo__mul_q16(st->k0_q16, h_ee_t)));

        uint32_t m = prob.m;
        for (j = 0U; j < ROBO_ARM6_DOF; ++j) {
            prob.con[m].a_q16[j] = fk.J_ee[2][j];
        }
        prob.con[m].b_q16          = rhs_ee;
        prob.con[m].slack_max_q16  = INT64_C(1310720);
        prob.con[m].slack_lin_q16  = INT64_C(134217728);
        prob.con[m].slack_quad_q16 = INT64_C(134217728);
        prob.con[m].is_hard        = 0U;
        prob.m                     = m + 1U;

        /* Wrist joint 4 height guard above table */
        int64_t z_floor_w = robo__sat_add_i64(st->table_z_q16, st->links[4].cap_radius_q16);
        int64_t h_w       = robo__sat_sub_i64(fk.joint_pos[4].z, z_floor_w);
        int64_t Jw[3][ROBO_ARM6_DOF];
        robo_arm6_point_jacobian(&fk, 4U, fk.joint_pos[4], Jw);
        robo_vec3_t vw     = robo_arm6_point_vel(&fk, st->dq_q16, 4U, fk.joint_pos[4]);
        robo_vec3_t jdw    = robo_arm6_point_jdot_dq(&fk, st->dq_q16, 4U, fk.joint_pos[4]);
        int64_t     v_appw = (vw.z < 0) ? robo__sat_neg_i64(vw.z) : 0;
        int64_t     pad_w  = robo__sat_add_i64(base_pad, robo__mul_q16(v_appw, tau_eff));
        int64_t     h_w_t  = robo__sat_sub_i64(h_w, pad_w);
        int64_t     hd_ww  = robo__sat_sub_i64(vw.z, st->dist_vel_q16);
        int64_t     rhs_w  = robo__sat_sub_i64(
            robo__sat_sub_i64(st->dist_acc_q16, jdw.z),
            robo__sat_add_i64(robo__mul_q16(st->k1_q16, hd_ww), robo__mul_q16(st->k0_q16, h_w_t)));

        m = prob.m;
        for (j = 0U; j < ROBO_ARM6_DOF; ++j) {
            prob.con[m].a_q16[j] = Jw[2][j];
        }
        prob.con[m].b_q16          = rhs_w;
        prob.con[m].slack_max_q16  = INT64_C(1310720);
        prob.con[m].slack_lin_q16  = INT64_C(134217728);
        prob.con[m].slack_quad_q16 = INT64_C(134217728);
        prob.con[m].is_hard        = 0U;
        prob.m                     = m + 1U;

        cert->min_contact_h_q16 = robo__min_i64(h_ee, h_w);
    }

    /* 5. Joint Torque Envelope Constraints: -tau_max <= M(q)*ddq + g(q) <= tau_max (12 rows, M <= 24) */
    for (i = 0U; i < ROBO_ARM6_DOF && prob.m + 1U < ROBO_QP_MAX_M; ++i) {
        /* 5% headroom for +/-20% payload mass mismatch on torque saturation */
        int64_t tau_eff_lim = robo__mul_q16(st->tau_max_q16[i], INT64_C(62259)); /* 0.95 * tau_max */
        int64_t neg_tau_lim = robo__sat_neg_i64(tau_eff_lim);
        uint32_t m_lo       = prob.m;
        uint32_t m_hi       = prob.m + 1U;
        for (j = 0U; j < ROBO_ARM6_DOF; ++j) {
            prob.con[m_lo].a_q16[j] = M[i][j];
            prob.con[m_hi].a_q16[j] = robo__sat_neg_i64(M[i][j]);
        }
        prob.con[m_lo].b_q16          = robo__sat_sub_i64(neg_tau_lim, grav[i]);
        prob.con[m_lo].slack_max_q16  = INT64_C(1310720);
        prob.con[m_lo].slack_lin_q16  = INT64_C(4194304);
        prob.con[m_lo].slack_quad_q16 = INT64_C(4194304);
        prob.con[m_lo].is_hard        = 0U;

        prob.con[m_hi].b_q16          = robo__sat_add_i64(neg_tau_lim, grav[i]);
        prob.con[m_hi].slack_max_q16  = INT64_C(1310720);
        prob.con[m_hi].slack_lin_q16  = INT64_C(4194304);
        prob.con[m_hi].slack_quad_q16 = INT64_C(4194304);
        prob.con[m_hi].is_hard        = 0U;
        prob.m                        = m_hi + 1U;
    }

    robo_qp_nd_cert_t qp_cert;
    robo_qp_nd_solve_ring0(&prob, &qp_cert);

    cert->verdict          = qp_cert.verdict;
    cert->n_constraints    = prob.m;
    cert->qp_iters         = qp_cert.iters_run;
    cert->feasible         = qp_cert.feasible;
    cert->slack_active     = qp_cert.slack_active;
    cert->fallback_engaged = qp_cert.fallback_engaged;
    cert->tube_margin_q16  = max_tube;
    cert->max_slack_q16    = qp_cert.max_slack_q16;
    cert->min_stl_h_q16    = robo__min_i64(
        robo__min_i64(cert->min_obs_h_q16, cert->min_self_h_q16),
        robo__min_i64(cert->min_contact_h_q16, cert->min_joint_h_q16));

    for (i = 0U; i < ROBO_ARM6_DOF; ++i) {
        cert->ddq_safe_q16[i] = qp_cert.u_safe_q16[i];
        int64_t mddq          = robo__dot_nd_q16(M[i], cert->ddq_safe_q16, ROBO_ARM6_DOF);
        int64_t tau_i         = robo__sat_add_i64(mddq, grav[i]);
        cert->tau_safe_q16[i] = robo__clamp_i64(tau_i, robo__sat_neg_i64(st->tau_max_q16[i]),
                                                st->tau_max_q16[i]);
    }
}

#ifdef __cplusplus
}
#endif

#endif /* ROBO_ARM6_RING0_H */
