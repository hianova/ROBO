# ROBO v3.0.0 — Embodied AI Ring-0 Physical Safety Verification Report

- **Filter Kernels**: Pure-GPR 2D Active-Set QP + General N-DoF ($N \le 8, M \le 24$) ADMM Slack-CBF-QP + 3D 6-DoF Spatial Tube-CBF (`0` FPU, `0` `udiv`/`sdiv`, `0` heap, `0` `blr`/`br`, `0` NSW UB)
- **Ring-0 Latency**: Planar **`3615.8 ns/step`** | 2-Link Rigid Arm **`1085.5 ns/step`** | 6-DoF Spatial Tube-CBF **`133375.5 ns/step`** (1 kHz deadline = `1,000,000 ns`)

## 1. SafeVLA-Bench 6-DoF Spatial Manipulator (10 kHz `float64` RK4 vs 1 kHz Ring-0)

- **Adversarial Stress**: `25` episodes × `2.0 s` (`500000` RK4 steps/ctrl), `30–50 ms` VLA action-delay buffer, `±20%` payload mass mismatch, sensor noise, capsule-to-sphere obstacles, capsule-to-capsule self-collision, and compliant table contact force (`F_max = 25 N`).

| Controller | Task Success Rate | SBU Rate (Target 0%) | Obs Collisions | Self-Collisions | Contact-Force Viol | Joint Viol | Mean VSI | STL Robustness |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Unfiltered VLA (Adversarial) | 100.0% (25/25) | **100.0% (25/25)** | 25 / 25 | 25 / 25 | 25 / 25 | 0 / 25 | 1.8381 | -0.1059 |
| Naive Box Torque Clipping | 100.0% (25/25) | **100.0% (25/25)** | 25 / 25 | 25 / 25 | 25 / 25 | 0 / 25 | 1.8382 | -0.1059 |
| ROBO v3.0 Exact-QP Ring-0 CBF | 52.0% (13/25) | **0.0% (0/25)** | 0 / 25 | 0 / 25 | 0 / 25 | 0 / 25 | 0.0000 | +0.0125 |

## 2. Moving-Obstacle Planar VLA Contrast (`500` eps × `200` steps = `100,000` steps)

| Controller | Obstacle Collisions (eps) | Wall Breaches (eps) | Friction Cone Viol (steps) | ZMP Viol (steps) | Mean Final Goal Err (m) |
| :--- | ---: | ---: | ---: | ---: | ---: |
| Unfiltered VLA (Adversarial) | 474 / 500 | 0 / 500 | 2131 | 27098 | 0.3867 |
| Naive Box Torque Clipping | 474 / 500 | 73 / 500 | 0 | 26738 | 0.3786 |
| ROBO v3.0 Exact-QP Ring-0 CBF | 0 / 500 | 0 / 500 | 0 | 0 | 0.3636 |

## 3. 2-Link Rigid-Body Manipulator Torque-Space CBF (`200` eps × `200` steps)

| Controller | Cartesian Obstacle Collisions (eps) | Joint Limit Breaches (eps) |
| :--- | ---: | ---: |
| Unfiltered VLA Torques | 100 / 200 | 0 / 200 |
| **ROBO v3.0 Torque-Space CBF** | **0 / 200** | **0 / 200** |

## 4. Formal Lean 4 Proofs, Z3 Falsification & FLOW/Auditor Verification
- **Pure-`Init` Lean 4 Proofs ([`lean/RoboSafety.lean`](lean/RoboSafety.lean))**: Verified via `lean -DwarningAsError=true` + `leanchecker RoboSafety` (`0 sorryAx`, `0 native_decide`, `0` custom axioms) covering Tube-CBF inductive & multi-step horizon invariance, Slack-QP guaranteed feasibility, emergency-brake box feasibility, Coulomb friction, ZMP envelope, compliant contact-force ceiling, and zero-SBU STL theorem.
- **ArenaSchema ([`ir/robo_wire.fspec`](ir/robo_wire.fspec))**: 64-byte `flow_cell` Z3 `QF_LIA` `unsat` (`0` overlap / `0` gap).
- **Z3 Soundness & Falsification ([`proofs/robo_falsify.smt2`](proofs/robo_falsify.smt2))**: 8 checks (`sat` anti-vacuity, `unsat` ECBF invariance, `sat` shock falsification, `unsat` halfspace separation, `unsat` Tube-CBF delay invariance, `sat` unmodeled delay falsification, `unsat` Slack-QP feasibility, `unsat` contact-force bound).
- **Clang LLVM SSA (`auditor bug --check-nsw --strict`)**: **`0` findings** across `src/robo_ring0.c`; `auditor scan`: **`VIOLATION=0`**.
