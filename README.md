# ROBO v3.0.0 — Ring-0 Embodied AI Physical Safety Filter

Deterministic Pure-GPR physical safety shield for Vision-Language-Action (VLA), contact locomotion, and 6-DoF spatial rigid-body manipulators built on [`../FLOW`](../FLOW/README.md) (`core/flow.h`, `core/flow_cholesky_sdd.h`), aligned with SafeVLA-Bench & MuJoCo Playground/SIMPLE safety standards:

- [`include/robo_dynamics_ring0.h`](include/robo_dynamics_ring0.h): Freestanding Q16.16 kinematics, `sincos`, Coulomb friction cone, and ZMP support-polygon limits (`0` FPU/SIMD, `0` `udiv`/`sdiv`, `0` heap, `0` `blr`/`br`).
- [`include/robo_cbf_ring0.h`](include/robo_cbf_ring0.h): Exact Pure-GPR 2D Active-Set QP solver (`robo__solve_qp2d_ring0`) + robust discrete ECBF with obstacle-accel and `dt` braking margins.
- [`include/robo_arm2_ring0.h`](include/robo_arm2_ring0.h): Analytical 2-Link rigid-body dynamics `M(q)ddq + C(q,dq)dq + g(q) = tau` mapped through `J(q)` into torque-space ECBF & joint-limit halfspaces.
- [`include/robo_qp_nd_ring0.h`](include/robo_qp_nd_ring0.h): General N-DoF ($N \le 8, M \le 24$) Pure-GPR ADMM Slack-CBF-QP solver with `flow_cholesky_sdd_floor_i64` PSD certification, L1+L2 slack $\delta \ge 0$, and guaranteed box-feasible emergency-brake fallback.
- [`include/robo_arm6_ring0.h`](include/robo_arm6_ring0.h): Pure-GPR 3D unit-quaternion / SO(3), 6-DoF spatial kinematics & analytical Jacobian $J(q)$ + $\dot{J}\dot{q}$, exact capsule-to-sphere & capsule-to-capsule self-collision distances/gradients, table contact-force envelope, and explicit 30–50 ms VLA action-delay Tube-CBF inflation.
- [`lean/RoboSafety.lean`](lean/RoboSafety.lean), [`ir/robo_wire.fspec`](ir/robo_wire.fspec) & [`proofs/robo_falsify.smt2`](proofs/robo_falsify.smt2): Pure-`Init` Lean 4 theorems (`lean` + `leanchecker`, `0 sorryAx`, `0 native_decide`), `FLOW v10` 64-byte ArenaSchema, 8 Z3 soundness/falsification checks, and `0` `--check-nsw --strict` findings.
- [`src/robo_cli.c`](src/robo_cli.c) & [`oracles/robo_rk4_oracle.py`](oracles/robo_rk4_oracle.py): 10 kHz continuous-time `float64` RK4 vs 1 kHz Ring-0 SafeVLA-Bench adversarial evaluation (30–50 ms action delay, $\pm 20\%$ payload mass mismatch, sensor noise) generating [`ROBO_VERIFICATION_REPORT.md`](ROBO_VERIFICATION_REPORT.md).

```sh
make clean all test demo
```
