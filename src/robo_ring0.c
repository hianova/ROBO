#include "robo_dynamics_ring0.h"
#include "robo_cbf_ring0.h"
#include "robo_arm2_ring0.h"

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
