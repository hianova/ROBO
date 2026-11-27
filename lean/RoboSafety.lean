/-!
ROBO v3.0.0 — Pure `Init` Formal Safety Proofs for Discrete Tube-CBF,
Slack-QP Feasibility, Contact Force, Friction/ZMP, and Joint-Box Invariants.

Zero external packages, kernel-checked only (no unverified axioms or holes).
Verified via `lean` + `leanchecker`.
-/

namespace RoboSafety

/-- Clamp an integer `x` into `[lo, hi]`. -/
def clampInt (x lo hi : Int) : Int :=
  if x < lo then lo else if x > hi then hi else x

/-- Discrete-time second-order barrier step with action/model disturbance. -/
def barrierNext (h dh u dist : Int) : Int :=
  h + dh + u + dist

/-- Discrete delay + acceleration tube margin (scaled integer units). -/
def tubePad (vApp tauDelay aPertTerm distPos : Int) : Int :=
  vApp * tauDelay + aPertTerm + distPos

/-- Emergency-brake clamping is always box-feasible whenever `u_lo <= u_hi`. -/
theorem emergency_brake_box_feasible (u_brake u_lo u_hi : Int) (h_box : u_lo ≤ u_hi) :
    u_lo ≤ clampInt u_brake u_lo u_hi ∧ clampInt u_brake u_lo u_hi ≤ u_hi := by
  unfold clampInt
  split
  · exact ⟨Int.le_refl u_lo, h_box⟩
  · rename_i h_not_lt
    split
    · exact ⟨h_box, Int.le_refl u_hi⟩
    · rename_i h_not_gt
      have h1 : u_lo ≤ u_brake := Int.not_lt.mp h_not_lt
      have h2 : u_brake ≤ u_hi := Int.not_lt.mp h_not_gt
      exact ⟨h1, h2⟩

/-- Slack-QP feasibility: for any box-clamped `u_safe` and any linear constraint
    `dot >= b`, choosing `delta = max (b - dot) 0` makes `delta >= 0` and
    `dot + delta >= b`. -/
theorem slack_qp_always_feasible (dot b : Int) :
    ∃ delta : Int, delta ≥ 0 ∧ dot + delta ≥ b := by
  by_cases h : dot ≥ b
  · exact ⟨0, Int.le_refl 0, by omega⟩
  · refine ⟨b - dot, ?_, ?_⟩
    · omega
    · omega

/-- Bounded-slack dichotomy: either the soft constraint is satisfied within
    `[0, s_max]`, or the required slack strictly exceeds `s_max` (triggering
    the box-feasible emergency-brake fallback). -/
theorem bounded_slack_dichotomy (dot b s_max : Int) (hs : s_max ≥ 0) :
    (∃ delta : Int, 0 ≤ delta ∧ delta ≤ s_max ∧ dot + delta ≥ b) ∨ (b - dot > s_max) := by
  by_cases h1 : b - dot ≤ s_max
  · left
    by_cases h2 : dot ≥ b
    · exact ⟨0, Int.le_refl 0, hs, by omega⟩
    · refine ⟨b - dot, by omega, h1, by omega⟩
  · right
    omega

/-- Non-negativity of the latency Tube-CBF margin under non-negative speeds and bounds. -/
theorem tube_pad_nonneg (vApp tauDelay aPertTerm distPos : Int)
    (hv : vApp ≥ 0) (ht : tauDelay ≥ 0) (ha : aPertTerm ≥ 0) (hp : distPos ≥ 0) :
    tubePad vApp tauDelay aPertTerm distPos ≥ 0 := by
  unfold tubePad
  have hmul : vApp * tauDelay ≥ 0 := Int.mul_nonneg hv ht
  omega

/-- Single-step Tube-CBF inductive forward invariance under bounded delay/disturbance.
    If the barrier starts inside the tube-inflated safe set (`h_k >= pad`),
    the discrete CBF condition holds with decay `decay <= h_k - pad` and margin `d_max`,
    and the actual disturbance/delay error satisfies `dist >= -d_max`,
    then `h_{k+1} >= pad`. -/
theorem tube_cbf_inductive_step (h_k dh_k u_cmd dist d_max decay pad : Int)
    (h_decay_le : decay ≤ h_k - pad)
    (h_cbf : h_k + dh_k + u_cmd ≥ h_k - decay + d_max)
    (h_dist : dist ≥ -d_max) :
    barrierNext h_k dh_k u_cmd dist ≥ pad := by
  unfold barrierNext
  omega

/-- Multi-step inductive invariance over arbitrary finite horizon `N : Nat`. -/
theorem tube_cbf_horizon_invariant (pad : Int) (h_pad : pad ≥ 0)
    (h_seq dh_seq u_seq dist_seq d_max_seq decay_seq : Nat → Int)
    (h_step_def : ∀ k, h_seq (k + 1) = barrierNext (h_seq k) (dh_seq k) (u_seq k) (dist_seq k))
    (h_init : h_seq 0 ≥ pad)
    (h_decay : ∀ k, decay_seq k ≤ h_seq k - pad)
    (h_cbf : ∀ k, h_seq k + dh_seq k + u_seq k ≥ h_seq k - decay_seq k + d_max_seq k)
    (h_dist : ∀ k, dist_seq k ≥ -(d_max_seq k)) :
    ∀ k : Nat, h_seq k ≥ pad ∧ h_seq k ≥ 0 := by
  intro k
  induction k with
  | zero =>
    exact ⟨h_init, Int.le_trans h_pad h_init⟩
  | succ m ih =>
    have h_next : h_seq (m + 1) ≥ pad := by
      rw [h_step_def m]
      exact tube_cbf_inductive_step (h_seq m) (dh_seq m) (u_seq m) (dist_seq m)
        (d_max_seq m) (decay_seq m) pad (h_decay m) (h_cbf m) (h_dist m)
    exact ⟨h_next, Int.le_trans h_pad h_next⟩

/-- Anti-vacuity witness: concrete values satisfy all hypotheses of
    `tube_cbf_inductive_step` and yield `h_{k+1} >= pad > 0`. -/
theorem tube_cbf_witness_nonvacuous :
    let h_k : Int := 65536
    let pad : Int := 6554
    let decay : Int := 14745
    let dh_k : Int := -3276
    let d_max : Int := 3276
    let u_cmd : Int := 0
    let dist : Int := -3276
    decay ≤ h_k - pad ∧
    h_k + dh_k + u_cmd ≥ h_k - decay + d_max ∧
    dist ≥ -d_max ∧
    barrierNext h_k dh_k u_cmd dist ≥ pad := by
  decide

/-- Coulomb friction cone invariance from Ring-0 tangential force bounds. -/
theorem friction_cone_invariant (ft m u_x f_lim : Int)
    (h_lo : -f_lim - ft ≤ m * u_x) (h_hi : m * u_x ≤ f_lim - ft) :
    -f_lim ≤ ft + m * u_x ∧ ft + m * u_x ≤ f_lim := by
  omega

/-- ZMP support-polygon invariance from Ring-0 acceleration bounds. -/
theorem zmp_envelope_invariant (zmp0 hg u_x w_half : Int)
    (h_lo : -w_half - zmp0 ≤ hg * u_x) (h_hi : hg * u_x ≤ w_half - zmp0) :
    -w_half ≤ zmp0 + hg * u_x ∧ zmp0 + hg * u_x ≤ w_half := by
  omega

/-- Compliant table contact-force bound: when penetration `pen` is kept within
    `k_env * pen <= f_max` and upward velocity damping `b_damp * v_up >= 0`,
    total contact force never exceeds `f_max`. -/
theorem contact_force_bounded (k_env pen b_damp v_up f_max : Int)
    (h_pen_bound : k_env * pen ≤ f_max)
    (hb : b_damp ≥ 0) (hv : v_up ≥ 0) :
    k_env * pen - b_damp * v_up ≤ f_max := by
  have hd : b_damp * v_up ≥ 0 := Int.mul_nonneg hb hv
  omega

/-- Joint-angle box robustness: if `q_min <= q <= q_max`, both boundary margins
    `q - q_min` and `q_max - q` are non-negative. -/
theorem joint_box_stl_nonneg (q q_min q_max : Int)
    (h_lo : q_min ≤ q) (h_hi : q ≤ q_max) :
    q - q_min ≥ 0 ∧ q_max - q ≥ 0 := by
  omega

/-- Zero SBU (Success-But-Unsafe) theorem: if every safety barrier across an
    episode has non-negative STL robustness (`stl_min >= 0`), the episode cannot
    be counted as Success-But-Unsafe (`sbu_flag = 0`). -/
theorem stl_nonneg_implies_zero_sbu (task_ok : Bool) (stl_min : Int)
    (h_safe : stl_min ≥ 0) :
    (if task_ok && decide (stl_min < 0) then (1 : Nat) else 0) = 0 := by
  have hnot : ¬ (stl_min < 0) := Int.not_lt.mpr h_safe
  simp [hnot]

end RoboSafety
