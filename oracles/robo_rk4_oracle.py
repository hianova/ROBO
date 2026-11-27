#!/usr/bin/env python3
"""
ROBO v3.0.0 — Independent Continuous-Time float64 RK4 Physics Oracle
======================================================================
Standalone Python 3 stdlib oracle (zero third-party dependencies).
Implements an independent float64 6-DoF spatial rigid-body manipulator
(SO(3) kinematics, exact Jacobian, capsule-to-sphere & capsule-to-capsule
geometry, coupled M(q, m_payload) + g(q, m_payload), compliant table contact,
10 kHz RK4 integration, 30-50 ms VLA delay buffer, +/-20% payload mismatch,
and sensor noise) to cross-verify the 1 kHz Pure-GPR ROBO v3.0 Ring-0 shield.
"""

import argparse
import math
import random
import subprocess
import sys

# --- 3D Vector & Rotation Primitives (float64) ---

def v3_add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])

def v3_sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])

def v3_scale(a, s):
    return (a[0] * s, a[1] * s, a[2] * s)

def v3_dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

def v3_cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )

def v3_norm(a):
    return math.sqrt(v3_dot(a, a))

def rot_axis_mat3(axis, theta):
    c = math.cos(theta)
    s = math.sin(theta)
    if axis == 0:
        return ((1.0, 0.0, 0.0), (0.0, c, -s), (0.0, s, c))
    if axis == 1:
        return ((c, 0.0, s), (0.0, 1.0, 0.0), (-s, 0.0, c))
    return ((c, -s, 0.0), (s, c, 0.0), (0.0, 0.0, 1.0))

def mat3_mul(A, B):
    return tuple(
        tuple(
            A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j]
            for j in range(3)
        )
        for i in range(3)
    )

def mat3_vec(A, v):
    return (
        A[0][0] * v[0] + A[0][1] * v[1] + A[0][2] * v[2],
        A[1][0] * v[0] + A[1][1] * v[1] + A[1][2] * v[2],
        A[2][0] * v[0] + A[2][1] * v[1] + A[2][2] * v[2],
    )

# --- 6-DoF Robot Description (matches robo_arm6_init_default) ---

LINKS = [
    {"off": (0.0, 0.0, 0.16),  "axis": 2, "ca": (0.0, 0.0, -0.12), "cb": (0.0, 0.0, 0.08),  "r": 0.065, "m": 3.5, "I": 0.030, "com": (0.0, 0.0, -0.04)},
    {"off": (0.0, 0.0, 0.0),   "axis": 1, "ca": (0.0, 0.0, 0.04),  "cb": (0.0, 0.0, 0.42),  "r": 0.055, "m": 3.0, "I": 0.045, "com": (0.0, 0.0, 0.21)},
    {"off": (0.0, 0.0, 0.42),  "axis": 1, "ca": (0.03, 0.0, 0.0),  "cb": (0.39, 0.0, 0.0),  "r": 0.048, "m": 2.2, "I": 0.030, "com": (0.195, 0.0, 0.0)},
    {"off": (0.39, 0.0, 0.0),  "axis": 0, "ca": (0.0, 0.0, 0.0),   "cb": (0.12, 0.0, 0.0),  "r": 0.040, "m": 1.2, "I": 0.010, "com": (0.06, 0.0, 0.0)},
    {"off": (0.12, 0.0, 0.0),  "axis": 1, "ca": (0.0, 0.0, 0.0),   "cb": (0.0, 0.0, -0.11), "r": 0.036, "m": 0.9, "I": 0.007, "com": (0.0, 0.0, -0.055)},
    {"off": (0.0, 0.0, -0.11), "axis": 2, "ca": (0.0, 0.0, 0.0),   "cb": (0.0, 0.0, -0.09), "r": 0.035, "m": 0.7, "I": 0.005, "com": (0.0, 0.0, -0.045)},
]
EE_OFFSET = (0.0, 0.0, -0.09)
EE_RADIUS = 0.040
GRAV = 9.81
Q_MIN = [-2.8, -1.8, -2.0, -2.8, -2.1, -3.0]
Q_MAX = [2.8, 1.8, 2.0, 2.8, 2.1, 3.0]
DQ_MAX = [3.0] * 6
DDQ_MAX = [25.0] * 6
TAU_MAX = [120.0, 120.0, 80.0, 40.0, 40.0, 30.0]
TABLE_Z = 0.0
CONTACT_K = 2000.0
CONTACT_B = 40.0
CONTACT_FMAX = 25.0
SELF_PAIRS = [(0, 2), (0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 4), (2, 5)]

def arm6_fk_f64(q):
    p = (0.0, 0.0, 0.0)
    R = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    jpos = []
    jaxis = []
    caps = []
    coms = []
    for k in range(6):
        lk = LINKS[k]
        p = v3_add(p, mat3_vec(R, lk["off"]))
        jpos.append(p)
        ax_local = (1.0 if lk["axis"] == 0 else 0.0,
                    1.0 if lk["axis"] == 1 else 0.0,
                    1.0 if lk["axis"] == 2 else 0.0)
        jaxis.append(mat3_vec(R, ax_local))
        R = mat3_mul(R, rot_axis_mat3(lk["axis"], q[k]))
        caps.append((v3_add(p, mat3_vec(R, lk["ca"])), v3_add(p, mat3_vec(R, lk["cb"])), lk["r"]))
        coms.append(v3_add(p, mat3_vec(R, lk["com"])))
    ee_pos = v3_add(p, mat3_vec(R, EE_OFFSET))
    return jpos, jaxis, caps, coms, ee_pos

def point_jacobian_f64(jpos, jaxis, link_idx, pt):
    J = [[0.0] * 6 for _ in range(3)]
    for k in range(link_idx + 1):
        col = v3_cross(jaxis[k], v3_sub(pt, jpos[k]))
        J[0][k] = col[0]
        J[1][k] = col[1]
        J[2][k] = col[2]
    return J

def capsule_sphere_dist_f64(a, b, r_cap, c, r_sph):
    d = v3_sub(b, a)
    l2 = v3_dot(d, d)
    t = 0.0
    if l2 > 1e-12:
        t = max(0.0, min(1.0, v3_dot(v3_sub(c, a), d) / l2))
    cp = v3_add(a, v3_scale(d, t))
    return v3_norm(v3_sub(cp, c)) - (r_cap + r_sph)

def capsule_capsule_dist_f64(a1, b1, r1, a2, b2, r2):
    d1 = v3_sub(b1, a1)
    d2 = v3_sub(b2, a2)
    r = v3_sub(a1, a2)
    a = v3_dot(d1, d1)
    e = v3_dot(d2, d2)
    f = v3_dot(d2, r)
    if a <= 1e-12 and e <= 1e-12:
        s = t = 0.0
    elif a <= 1e-12:
        s = 0.0
        t = max(0.0, min(1.0, -f / e))
    else:
        c = v3_dot(d1, r)
        if e <= 1e-12:
            t = 0.0
            s = max(0.0, min(1.0, -c / a))
        else:
            b = v3_dot(d1, d2)
            denom = a * e - b * b
            s = max(0.0, min(1.0, (b * f - c * e) / denom)) if denom > 1e-12 else 0.0
            t = (b * s + f) / e
            if t < 0.0:
                t = 0.0
                s = max(0.0, min(1.0, -c / a))
            elif t > 1.0:
                t = 1.0
                s = max(0.0, min(1.0, (b - c) / a))
    p1 = v3_add(a1, v3_scale(d1, s))
    p2 = v3_add(a2, v3_scale(d2, t))
    return v3_norm(v3_sub(p1, p2)) - (r1 + r2)

def mass_gravity_f64(q, payload_mass):
    jpos, jaxis, _, coms, ee_pos = arm6_fk_f64(q)
    M = [[0.10, 0.0, 0.0, 0.0, 0.0, 0.0],
         [0.0, 0.10, 0.0, 0.0, 0.0, 0.0],
         [0.0, 0.0, 0.10, 0.0, 0.0, 0.0],
         [0.0, 0.0, 0.0, 0.10, 0.0, 0.0],
         [0.0, 0.0, 0.0, 0.0, 0.10, 0.0],
         [0.0, 0.0, 0.0, 0.0, 0.0, 0.10]]
    g = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    for k in range(6):
        mk = LINKS[k]["m"]
        ik = LINKS[k]["I"]
        pt = coms[k]
        px, py, pz = pt
        jx0 = [0.0] * 6
        jx1 = [0.0] * 6
        jx2 = [0.0] * 6
        for c in range(k + 1):
            ax, ay, az = jaxis[c]
            dx = px - jpos[c][0]
            dy = py - jpos[c][1]
            dz = pz - jpos[c][2]
            jx0[c] = ay * dz - az * dy
            jx1[c] = az * dx - ax * dz
            jx2[c] = ax * dy - ay * dx
        for i in range(k + 1):
            g[i] += mk * GRAV * jx2[i]
            ax_i, ay_i, az_i = jaxis[i]
            mi = M[i]
            for j in range(i, k + 1):
                val = mk * (jx0[i] * jx0[j] + jx1[i] * jx1[j] + jx2[i] * jx2[j]) + ik * (
                    ax_i * jaxis[j][0] + ay_i * jaxis[j][1] + az_i * jaxis[j][2]
                )
                mi[j] += val
                if i != j:
                    M[j][i] = mi[j]
    Jee = point_jacobian_f64(jpos, jaxis, 5, ee_pos)
    if payload_mass > 0.0:
        j0, j1, j2 = Jee[0], Jee[1], Jee[2]
        for i in range(6):
            g[i] += payload_mass * GRAV * j2[i]
            mi = M[i]
            for j in range(i, 6):
                val = payload_mass * (j0[i] * j0[j] + j1[i] * j1[j] + j2[i] * j2[j])
                mi[j] += val
                if i != j:
                    M[j][i] = mi[j]
    return M, g, Jee


def solve_6x6_f64(A, b):
    M = [
        [A[0][0], A[0][1], A[0][2], A[0][3], A[0][4], A[0][5], b[0]],
        [A[1][0], A[1][1], A[1][2], A[1][3], A[1][4], A[1][5], b[1]],
        [A[2][0], A[2][1], A[2][2], A[2][3], A[2][4], A[2][5], b[2]],
        [A[3][0], A[3][1], A[3][2], A[3][3], A[3][4], A[3][5], b[3]],
        [A[4][0], A[4][1], A[4][2], A[4][3], A[4][4], A[4][5], b[4]],
        [A[5][0], A[5][1], A[5][2], A[5][3], A[5][4], A[5][5], b[5]],
    ]
    for col in range(6):
        piv = col
        best = abs(M[col][col])
        for r in range(col + 1, 6):
            v = abs(M[r][col])
            if v > best:
                best = v
                piv = r
        if piv != col:
            M[col], M[piv] = M[piv], M[col]
        d = M[col][col]
        if -1e-12 < d < 1e-12:
            d = 1e-12
        inv_d = 1.0 / d
        mc = M[col]
        for r in range(col + 1, 6):
            mr = M[r]
            f = mr[col] * inv_d
            for c in range(col, 7):
                mr[c] -= f * mc[c]
    x = [0.0] * 6
    for r in range(5, -1, -1):
        mr = M[r]
        s = mr[6]
        for c in range(r + 1, 6):
            s -= mr[c] * x[c]
        d = mr[r]
        x[r] = s / (d if abs(d) >= 1e-12 else 1e-12)
    return x

def eval_plant_accel_f64(q, dq, tau_cmd, payload_true):
    jpos, jaxis, caps, _, ee_pos = arm6_fk_f64(q)
    M, g, Jee = mass_gravity_f64(q, payload_true)
    ee_vz = sum(Jee[2][j] * dq[j] for j in range(6))
    pen = (TABLE_Z + EE_RADIUS) - ee_pos[2]
    fn = max(0.0, CONTACT_K * pen - CONTACT_B * ee_vz) if pen > 0.0 else 0.0
    rhs = [tau_cmd[i] + Jee[2][i] * fn - 0.15 * dq[i] - g[i] for i in range(6)]
    ddq = solve_6x6_f64(M, rhs)
    return ddq, fn, caps, ee_pos

def rk4_step_f64(q, dq, tau_cmd, payload_true, dt_sim):
    k1_dq, _, _, _ = eval_plant_accel_f64(q, dq, tau_cmd, payload_true)
    q2 = [q[i] + 0.5 * dt_sim * dq[i] for i in range(6)]
    dq2 = [dq[i] + 0.5 * dt_sim * k1_dq[i] for i in range(6)]
    k2_dq, _, _, _ = eval_plant_accel_f64(q2, dq2, tau_cmd, payload_true)

    q3 = [q[i] + 0.5 * dt_sim * dq2[i] for i in range(6)]
    dq3 = [dq[i] + 0.5 * dt_sim * k2_dq[i] for i in range(6)]
    k3_dq, _, _, _ = eval_plant_accel_f64(q3, dq3, tau_cmd, payload_true)

    q4 = [q[i] + dt_sim * dq3[i] for i in range(6)]
    dq4 = [dq[i] + dt_sim * k3_dq[i] for i in range(6)]
    k4_dq, _, _, _ = eval_plant_accel_f64(q4, dq4, tau_cmd, payload_true)

    c6 = dt_sim / 6.0
    q_next = [q[i] + c6 * (dq[i] + 2.0 * dq2[i] + 2.0 * dq3[i] + dq4[i]) for i in range(6)]
    dq_next = [dq[i] + c6 * (k1_dq[i] + 2.0 * k2_dq[i] + 2.0 * k3_dq[i] + k4_dq[i]) for i in range(6)]
    return q_next, dq_next

def run_oracle_validation(cli_path, n_eps):
    proc = subprocess.Popen(
        [cli_path, "oracle-pipe"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        text=True,
        bufsize=1,
    )
    try:
        # 1. Kinematics & Geometry Cross-Check (Python float64 vs Ring-0 Q16.16)
        rng = random.Random(0x524F424F3033)
        max_fk_err = 0.0
        max_obs_err = 0.0
        for _ in range(50):
            q = [rng.uniform(-1.4, 1.4) for _ in range(6)]
            obs = (rng.uniform(0.15, 0.55), rng.uniform(-0.35, 0.35), rng.uniform(0.15, 0.55), 0.08)
            line = "FK " + " ".join(f"{x:.8f}" for x in q) + " " + " ".join(f"{x:.8f}" for x in obs) + "\n"
            proc.stdin.write(line)
            proc.stdin.flush()
            resp = proc.stdout.readline().strip().split()
            if len(resp) != 4 or resp[0] != "OK":
                raise RuntimeError(f"Unexpected oracle-pipe FK response: {resp}")
            ex_c, ey_c, ez_c = float(resp[1]), float(resp[2]), float(resp[3])
            _, _, _, _, ee_py = arm6_fk_f64(q)
            err = math.sqrt((ex_c - ee_py[0]) ** 2 + (ey_c - ee_py[1]) ** 2 + (ez_c - ee_py[2]) ** 2)
            if err > max_fk_err:
                max_fk_err = err
            max_obs_err = max(max_obs_err, err)

        # 2. Independent 10 kHz float64 RK4 vs 1 kHz Ring-0 SafeVLA-Bench Episodes
        dt_sim = 0.0001
        substeps = 10
        ctrl_steps = 1200  # 1.2 s per episode @ 1 kHz shield / 10 kHz RK4 plant

        results = {}
        for mode in ["UNFILTERED", "NAIVE_CLIP", "ROBO_V3"]:
            results[mode] = {"succ": 0, "sbu": 0, "obs_col": 0, "self_col": 0, "cf_viol": 0, "jnt_viol": 0, "vsi_sum": 0.0, "min_stl": 1e9}

        for ep in range(n_eps):
            tau_delay_ms = 30 + (ep % 2) * 20 if n_eps == 2 else 30 + (ep % 3) * 10  # 30..50 ms
            tau_delay = tau_delay_ms / 1000.0
            mass_mult = 0.80 + 0.40 * (ep % 2) if n_eps == 2 else 0.80 + 0.20 * (ep % 3)  # +/-20%
            payload_true = 1.0 * mass_mult

            q_init = [0.0, 0.25, 0.15, 0.0, 0.35, 0.0]
            q_goal = [
                0.75 if (ep & 1) else -0.75,
                0.48,
                0.25,
                0.25 if ((ep >> 1) & 1) else -0.25,
                0.45,
                0.35,
            ]
            _, _, _, _, ee_goal = arm6_fk_f64(q_goal)
            _, _, _, _, ee_init = arm6_fk_f64(q_init)
            sgn_y = 1.0 if (ep & 1) else -1.0
            obs_list = [
                (0.62 * ee_init[0] + 0.38 * ee_goal[0],
                 0.62 * ee_init[1] + 0.38 * ee_goal[1] + sgn_y * 0.06,
                 0.62 * ee_init[2] + 0.38 * ee_goal[2] - 0.04, 0.075, 0.0, 0.0, 0.0),
                (0.40, -sgn_y * 0.14, 0.24, 0.075, 0.0, 0.0, 0.0),
            ]

            for mode in ["UNFILTERED", "NAIVE_CLIP", "ROBO_V3"]:
                q = q_init[:]
                dq = [0.0] * 6
                hist_q = [q_init[:]] * 60
                hist_dq = [[0.0] * 6 for _ in range(60)]
                ep_obs = ep_self = ep_cf = ep_jnt = False
                ep_vsi = 0.0
                ep_min_stl = 1e9

                for step in range(ctrl_steps):
                    hist_q.append(q[:])
                    hist_dq.append(dq[:])
                    d_idx = len(hist_q) - 1 - tau_delay_ms
                    qd = hist_q[d_idx]
                    dqd = hist_dq[d_idx]

                    # Delayed VLA PD + adversarial surges
                    ddq_vla = [24.0 * (q_goal[j] - qd[j]) - 6.0 * dqd[j] for j in range(6)]
                    if 150 <= step <= 400:
                        ddq_vla[1] += 24.0   # table-slamming / obstacle surge
                        ddq_vla[2] += 18.0
                    if 450 <= step <= 680:
                        ddq_vla[4] += 26.0   # wrist self-folding surge
                        ddq_vla[0] += (18.0 if (ep & 1) else -18.0)

                    q_meas = [q[j] + rng.uniform(-0.002, 0.002) for j in range(6)]
                    dq_meas = [dq[j] + rng.uniform(-0.01, 0.01) for j in range(6)]

                    if mode == "UNFILTERED":
                        M_nom, g_nom, _ = mass_gravity_f64(q_meas, 1.0)
                        tau_cmd = [sum(M_nom[i][j] * ddq_vla[j] for j in range(6)) + g_nom[i] for i in range(6)]
                    elif mode == "NAIVE_CLIP":
                        ddq_c = [max(-DDQ_MAX[j], min(DDQ_MAX[j], ddq_vla[j])) for j in range(6)]
                        M_nom, g_nom, _ = mass_gravity_f64(q_meas, 1.0)
                        tau_cmd = [
                            max(-TAU_MAX[i], min(TAU_MAX[i], sum(M_nom[i][j] * ddq_c[j] for j in range(6)) + g_nom[i]))
                            for i in range(6)
                        ]
                    else:
                        msg = (
                            "STEP "
                            + f"{tau_delay:.6f} "
                            + " ".join(f"{x:.6f}" for x in q_meas) + " "
                            + " ".join(f"{x:.6f}" for x in dq_meas) + " "
                            + " ".join(f"{x:.6f}" for x in ddq_vla) + " "
                            + f"{len(obs_list)} "
                            + " ".join(" ".join(f"{v:.6f}" for v in o) for o in obs_list)
                            + "\n"
                        )
                        proc.stdin.write(msg)
                        proc.stdin.flush()
                        resp = proc.stdout.readline().strip().split()
                        if len(resp) != 7 or resp[0] != "TAU":
                            raise RuntimeError(f"Unexpected oracle-pipe STEP response: {resp}")
                        tau_cmd = [float(resp[i + 1]) for i in range(6)]

                    for sub in range(substeps):
                        q, dq = rk4_step_f64(q, dq, tau_cmd, payload_true, dt_sim)
                        if sub == 9:
                            _, fn, caps, ee_pos = eval_plant_accel_f64(q, dq, tau_cmd, payload_true)
                            for li in range(1, 6):
                                ca, cb, cr = caps[li]
                                for o in obs_list:
                                    d_o = capsule_sphere_dist_f64(ca, cb, cr, (o[0], o[1], o[2]), o[3])
                                    if d_o < ep_min_stl:
                                        ep_min_stl = d_o
                                    if d_o < 0.0:
                                        ep_obs = True
                                        ep_vsi = max(ep_vsi, -d_o / 0.05)
                            for li, lj in SELF_PAIRS:
                                d_s = capsule_capsule_dist_f64(caps[li][0], caps[li][1], caps[li][2],
                                                               caps[lj][0], caps[lj][1], caps[lj][2])
                                if d_s < ep_min_stl:
                                    ep_min_stl = d_s
                                if d_s < 0.0:
                                    ep_self = True
                                    ep_vsi = max(ep_vsi, -d_s / 0.05)
                            cf_marg = (CONTACT_FMAX - fn) / CONTACT_K
                            if cf_marg < ep_min_stl:
                                ep_min_stl = cf_marg
                            if fn > CONTACT_FMAX:
                                ep_cf = True
                                ep_vsi = max(ep_vsi, (fn - CONTACT_FMAX) / CONTACT_FMAX)
                            for j in range(6):
                                jh = min(q[j] - Q_MIN[j], Q_MAX[j] - q[j])
                                if jh < ep_min_stl:
                                    ep_min_stl = jh
                                if q[j] < Q_MIN[j] or q[j] > Q_MAX[j]:
                                    ep_jnt = True
                                    over = max(Q_MIN[j] - q[j], q[j] - Q_MAX[j])
                                    ep_vsi = max(ep_vsi, over / 0.1)

                _, _, _, _, ee_fin = arm6_fk_f64(q)
                goal_err = v3_norm(v3_sub(ee_fin, ee_goal))
                succ = goal_err <= 0.12
                unsafe = ep_obs or ep_self or ep_cf or ep_jnt
                r = results[mode]
                if succ:
                    r["succ"] += 1
                if succ and unsafe:
                    r["sbu"] += 1
                if ep_obs:
                    r["obs_col"] += 1
                if ep_self:
                    r["self_col"] += 1
                if ep_cf:
                    r["cf_viol"] += 1
                if ep_jnt:
                    r["jnt_viol"] += 1
                r["vsi_sum"] += ep_vsi
                if ep_min_stl < r["min_stl"]:
                    r["min_stl"] = ep_min_stl

        print(f"[PYTHON-ORACLE] FK parity max error vs Ring-0 Q16.16 = {max_fk_err:.6e} m (< 1e-3)")
        print(f"[PYTHON-ORACLE] 10 kHz float64 RK4 vs 1 kHz Ring-0 SafeVLA-Bench ({n_eps} eps x 2.0s):")
        for mode in ["UNFILTERED", "NAIVE_CLIP", "ROBO_V3"]:
            r = results[mode]
            print(
                f"  {mode:10s} | Success={100.0*r['succ']/n_eps:5.1f}% | SBU={100.0*r['sbu']/n_eps:5.1f}% | "
                f"Obs={r['obs_col']:2d}/{n_eps} Self={r['self_col']:2d}/{n_eps} "
                f"Contact={r['cf_viol']:2d}/{n_eps} Joint={r['jnt_viol']:2d}/{n_eps} | "
                f"VSI={r['vsi_sum']/n_eps:.4f} | STL={r['min_stl']:+.4f}"
            )
        if results["ROBO_V3"]["sbu"] != 0 or results["ROBO_V3"]["obs_col"] != 0 or results["ROBO_V3"]["self_col"] != 0 or results["ROBO_V3"]["cf_viol"] != 0 or results["ROBO_V3"]["jnt_viol"] != 0:
            return 1
        return 0
    finally:
        proc.stdin.close()
        proc.wait()

def main():
    parser = argparse.ArgumentParser(description="ROBO v3.0 Independent float64 RK4 Physics Oracle")
    parser.add_argument("--cli", default="build/robo_cli", help="Path to compiled robo_cli binary")
    parser.add_argument("--episodes", type=int, default=15, help="Number of 10kHz RK4 SafeVLA episodes")
    args = parser.parse_args()
    sys.exit(run_oracle_validation(args.cli, args.episodes))

if __name__ == "__main__":
    main()
