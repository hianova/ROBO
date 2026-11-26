# ROBO v2.0.0 — Embodied AI Ring-0 Physical Safety Verification Report

- **Filter Kernel**: Exact Pure-GPR 2D Active-Set QP + Robust Discrete ECBF + 2-Link Rigid-Body `M(q)ddq + C + g = tau` (`0` FPU, `0` `udiv`/`sdiv`, `0` heap)
- **Latency**: Planar **`1201.9 ns/step`** | 2-Link Rigid Arm **`550.7 ns/step`** (1 kHz deadline = `1,000,000 ns`)

## 1. Moving-Obstacle Planar VLA Contrast (`500` eps × `200` steps = `100,000` steps)

| Controller | Obstacle Collisions (eps) | Wall Breaches (eps) | Friction Cone Viol (steps) | ZMP Viol (steps) | Mean Final Goal Err (m) |
| :--- | ---: | ---: | ---: | ---: | ---: |
| Unfiltered VLA (Adversarial) | 474 / 500 | 0 / 500 | 2131 | 27098 | 0.3867 |
| Naive Box Torque Clipping | 474 / 500 | 73 / 500 | 0 | 26738 | 0.3786 |
| ROBO v2.0 Exact-QP Ring-0 CBF | 0 / 500 | 0 / 500 | 0 | 0 | 0.3636 |

## 2. 2-Link Rigid-Body Manipulator Torque-Space CBF (`200` eps × `200` steps)

| Controller | Cartesian Obstacle Collisions (eps) | Joint Limit Breaches (eps) |
| :--- | ---: | ---: |
| Unfiltered VLA Torques | 100 / 200 | 0 / 200 |
| **ROBO v2.0 Torque-Space CBF** | **0 / 200** | **0 / 200** |

## 3. FLOW v10 Formal Proofs, Theory Falsification & NSW-0 Audit
- **ArenaSchema ([`ir/robo_wire.fspec`](ir/robo_wire.fspec))**: 64-byte `flow_cell` Z3 `QF_LIA` `unsat` (`0` overlap / `0` gap).
- **Discrete ECBF Soundness & Falsification ([`proofs/robo_falsify.smt2`](proofs/robo_falsify.smt2))**: Anti-vacuity **`sat`**; `|dist| <= d_max` invariance **`unsat`**; unmodeled shock `|dist| > d_max` **`sat` falsified** (concrete breach witness); active-set separation **`unsat`**.
- **Clang LLVM SSA (`auditor bug --check-nsw`)**: **`0` signed-overflow UB findings** across `src/robo_ring0.c` (`0` undef, `0` `udiv`/`sdiv`, `0` FPU/SIMD).
