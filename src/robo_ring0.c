#include "robo_dynamics_ring0.h"
#include "robo_cbf_ring0.h"
#include "robo_arm2_ring0.h"
#include "robo_qp_nd_ring0.h"
#include "robo_arm6_ring0.h"

robo_safety_cert_t robo_cbf_filter_ring0_export(const robo_state_t *st, int64_t u_nom_x_q16,
                                                int64_t u_nom_y_q16) {
    return robo_cbf_filter_ring0(st, u_nom_x_q16, u_nom_y_q16);
}

void robo_kinematics_check_ring0_export(int64_t l1_q16, int64_t l2_q16, int64_t th1_q16,
                                        int64_t th2_q16, int64_t *ex_q16, int64_t *ey_q16) {
    robo_fk2_ring0(l1_q16, l2_q16, th1_q16, th2_q16, ex_q16, ey_q16);
}

robo_arm2_cert_t robo_arm2_cbf_ring0_export(const robo_arm2_state_t *st, int64_t tau1_nom_q16,
                                            int64_t tau2_nom_q16) {
    return robo_arm2_cbf_ring0(st, tau1_nom_q16, tau2_nom_q16);
}

void robo_qp_nd_solve_ring0_export(const robo_qp_nd_problem_t *prob, robo_qp_nd_cert_t *cert) {
    robo_qp_nd_solve_ring0(prob, cert);
}

void robo_arm6_init_default_export(robo_arm6_state_t *st) {
    robo_arm6_init_default(st);
}

void robo_arm6_fk_ring0_export(const robo_arm6_state_t *st, robo_arm6_fk_t *fk) {
    robo_arm6_fk_ring0(st, fk);
}

void robo_arm6_cbf_accel_ring0_export(const robo_arm6_state_t *st,
                                      const int64_t ddq_nom_q16[ROBO_ARM6_DOF],
                                      robo_arm6_cert_t *cert) {
    robo_arm6_cbf_accel_ring0(st, ddq_nom_q16, cert);
}
