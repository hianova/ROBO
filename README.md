# ROBO v2.0.0 — Ring-0 Embodied AI Physical Safety Filter

Deterministic Pure-GPR physical safety shield for Vision-Language-Action (VLA), contact locomotion, and rigid-body manipulators built on [`../FLOW`](../FLOW/README.md) (`core/flow.h`):

- [`include/robo_dynamics_ring0.h`](include/robo_dynamics_ring0.h): Freestanding Q16.16 kinematics, `sincos`, Coulomb friction cone, and ZMP support-polygon limits (`0` FPU/SIMD, `0` `udiv`/`sdiv`, `0` heap).
- [`include/robo_cbf_ring0.h`](include/robo_cbf_ring0.h): Exact Pure-GPR 2D Active-Set QP solver (`robo__solve_qp2d_ring0`) + robust discrete ECBF with obstacle-accel and `dt` braking margins.
- [`include/robo_arm2_ring0.h`](include/robo_arm2_ring0.h): Analytical 2-Link rigid-body dynamics `M(q)ddq + C(q,dq)dq + g(q) = tau` mapped through `J(q)` into torque-space ECBF & joint-limit halfspaces.
- [`ir/robo_wire.fspec`](ir/robo_wire.fspec) & [`proofs/robo_falsify.smt2`](proofs/robo_falsify.smt2): `FLOW v10` 64-byte ArenaSchema + Z3 discrete ECBF soundness & shock-falsification suite (`0` `--check-nsw` findings).
- [`src/robo_cli.c`](src/robo_cli.c): Exact QP optimality tests + planar & 2-link arm adversarial benchmarks generating [`ROBO_VERIFICATION_REPORT.md`](ROBO_VERIFICATION_REPORT.md).

```sh
make clean all test demo
```
