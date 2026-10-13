import Mathlib

/-!
# Lazy periodic rules keep a valid lower bound

Machine-checked core of the lazy-rule argument for the factor-2 shortest common
superstring approximation (openai/math result 128).

Words form an abstract type `W` with a length function; `V` is the finite set of
distinct substrings. `N r` is the number of occurrences of `r` in an arbitrary
common superstring. The facts about `N` that the paper proves for true occurrence
counts enter as explicit hypotheses:

* `R1`, `R2`, `R3`: the two cases of the paper's Lemma 2.1 for a periodic rule (2.6)
  and its equation (2.4);
* extension inequalities `N s ≥ ∑ N (c s)` and `N s ≥ ∑ N (s c)`, and `N s ≥ 1` for
  required strings.

The results proved here:

* `rule_sound`: a periodic rule evaluated on *any* vector that is at most `N` on the
  relevant words is sound (the paper proves this only for its own count vector);
* `certificate`: if every count is justified by sound bounds that read only longer
  words, every count is at most `N`. This is exactly what the implementation's
  independent checker verifies.
-/

open Finset

set_option linter.unusedSectionVars false

namespace LazySCS

variable {W : Type*} [DecidableEq W]

/-- Blocking functional `B(·; x) = ∑ b r * x r` over the finite word set `V`. -/
def B (V : Finset W) (b : W → ℕ) (x : W → ℕ) : ℕ := ∑ r ∈ V, b r * x r

/-- Monotonicity of `B` in the count vector, needing the order only where `b ≠ 0`. -/
lemma B_mono (V : Finset W) (b : W → ℕ) {x y : W → ℕ}
    (h : ∀ r ∈ V, b r ≠ 0 → x r ≤ y r) : B V b x ≤ B V b y := by
  unfold B
  apply Finset.sum_le_sum
  intro r hr
  by_cases h0 : b r = 0
  · simp [h0]
  · exact Nat.mul_le_mul_left _ (h r hr h0)

/-- **Rule soundness for any lower vector.** Let `x ≤ N` on every word that the rule
reads (`v` and every word with a nonzero blocking coefficient). If the rule's trigger
holds for `x`, its conclusion holds for the true counts `N`. -/
theorem rule_sound (V : Finset W) (N x : W → ℕ) (v w : W) (k : ℕ) (bv bw : W → ℕ)
    (hx : ∀ r ∈ V, (bv r ≠ 0 ∨ bw r ≠ 0) → x r ≤ N r) (hxv : x v ≤ N v)
    (R1 : B V bv N < N v → B V bw N + k + 1 ≤ N w)
    (R2 : ∀ r ∈ V, 1 ≤ bv r → k + 1 ≤ bw r)
    (R3 : B V bw N ≤ N w)
    (trig : B V bv x < x v) :
    B V bw x + k + 1 ≤ N w := by
  have hmono : B V bw x ≤ B V bw N := B_mono V bw fun r hr h0 => hx r hr (Or.inr h0)
  by_cases h : B V bv N < N v
  · have := R1 h
    linarith
  · push_neg at h
    have hlt : B V bv x < B V bv N := by linarith
    unfold B at hlt
    obtain ⟨r0, hr0, hlt0⟩ := Finset.exists_lt_of_sum_lt hlt
    have hb : 1 ≤ bv r0 := by
      rcases Nat.eq_zero_or_pos (bv r0) with h0 | h0
      · rw [h0] at hlt0
        simp at hlt0
      · exact h0
    have hxr : x r0 < N r0 := by
      by_contra hc
      push_neg at hc
      have := Nat.mul_le_mul_left (bv r0) hc
      linarith
    have hk := R2 r0 hr0 hb
    have key : B V bw x + bw r0 ≤ B V bw N := by
      unfold B
      have hsplit : ∑ r ∈ V, bw r * x r + bw r0
          = ∑ r ∈ V, (bw r * x r + if r = r0 then bw r else 0) := by
        rw [Finset.sum_add_distrib, Finset.sum_ite_eq' V r0, if_pos hr0]
      rw [hsplit]
      apply Finset.sum_le_sum
      intro r hr
      by_cases hre : r = r0
      · subst hre
        have h1 : bw r * (x r + 1) ≤ bw r * N r := Nat.mul_le_mul_left _ hxr
        have h2 : bw r * (x r + 1) = bw r * x r + bw r := by ring
        rw [if_pos rfl]
        linarith
      · simp only [hre, if_false, add_zero]
        by_cases h0 : bw r = 0
        · simp [h0]
        · exact Nat.mul_le_mul_left _ (hx r hr (Or.inr h0))
    linarith

/-- A bound for word `s` is sound if it never exceeds `N s` when evaluated on any
vector that is at most `N` on all strictly longer words of `V`. -/
def SoundFor (V : Finset W) (len : W → ℕ) (N : W → ℕ) (s : W) (g : (W → ℕ) → ℕ) : Prop :=
  ∀ x : W → ℕ, (∀ r ∈ V, len s < len r → x r ≤ N r) → g x ≤ N s

/-- Extension bounds (2.5) are sound, given the extension inequality for `N`. -/
theorem ext_sound (V : Finset W) (len N : W → ℕ) (s : W) (S : Finset W)
    (hS : ∀ r ∈ S, r ∈ V ∧ len s < len r) (hN : ∑ r ∈ S, N r ≤ N s) :
    SoundFor V len N s (fun x => ∑ r ∈ S, x r) := by
  intro x hx
  exact le_trans (Finset.sum_le_sum fun r hr => hx r (hS r hr).1 (hS r hr).2) hN

/-- The required-string bound is sound. -/
theorem req_sound (V : Finset W) (len N : W → ℕ) (s : W) (h : 1 ≤ N s) :
    SoundFor V len N s (fun _ => 1) := fun _ _ => h

/-- A periodic rule (2.6) as a bound for its word `w`: `B(w; x) + k + 1` when the
trigger holds for `x`, and `0` otherwise. It is sound whenever every word it reads
is strictly longer than `w`. -/
theorem rule_bound_sound (V : Finset W) (len N : W → ℕ) (v w : W) (k : ℕ) (bv bw : W → ℕ)
    (hv : v ∈ V) (hvlen : len w < len v)
    (hsupp : ∀ r ∈ V, (bv r ≠ 0 ∨ bw r ≠ 0) → len w < len r)
    (R1 : B V bv N < N v → B V bw N + k + 1 ≤ N w)
    (R2 : ∀ r ∈ V, 1 ≤ bv r → k + 1 ≤ bw r)
    (R3 : B V bw N ≤ N w) :
    SoundFor V len N w (fun x => if B V bv x < x v then B V bw x + k + 1 else 0) := by
  intro x hx
  dsimp only
  split_ifs with trig
  · exact rule_sound V N x v w k bv bw (fun r hr h => hx r hr (hsupp r hr h))
      (hx v hv hvlen) R1 R2 R3 trig
  · exact Nat.zero_le _

/-- **Certificate theorem.** Suppose every word `s ∈ V` has a list of sound bounds and
its count `m s` is at most one of them evaluated at `m`. Then `m s ≤ N s` for every
`s ∈ V`. The proof is induction on decreasing length; it uses the final counts
only, so it applies to any procedure (such as a fixpoint iteration) whose output
passes the check. -/
theorem certificate (V : Finset W) (len N m : W → ℕ)
    (bounds : W → List ((W → ℕ) → ℕ))
    (hsound : ∀ s ∈ V, ∀ g ∈ bounds s, SoundFor V len N s g)
    (hjust : ∀ s ∈ V, ∃ g ∈ bounds s, m s ≤ g m) :
    ∀ s ∈ V, m s ≤ N s := by
  have main : ∀ d : ℕ, ∀ s ∈ V, V.sup len - len s ≤ d → m s ≤ N s := by
    intro d
    induction d with
    | zero =>
      intro s hs hd
      obtain ⟨g, hg, hmg⟩ := hjust s hs
      refine le_trans hmg (hsound s hs g hg m ?_)
      intro r hr hlt
      have := Finset.le_sup (f := len) hr
      omega
    | succ d ih =>
      intro s hs hd
      obtain ⟨g, hg, hmg⟩ := hjust s hs
      refine le_trans hmg (hsound s hs g hg m ?_)
      intro r hr hlt
      apply ih r hr
      have := Finset.le_sup (f := len) hr
      omega
  intro s hs
  exact main _ s hs le_rfl

/-- Consequence for the lower bound: summing over one-letter words, `W' ≤ ∑ N c`,
which is the length of the superstring when `N` counts its occurrences. -/
theorem lower_bound (V : Finset W) (len N m : W → ℕ)
    (bounds : W → List ((W → ℕ) → ℕ))
    (hsound : ∀ s ∈ V, ∀ g ∈ bounds s, SoundFor V len N s g)
    (hjust : ∀ s ∈ V, ∃ g ∈ bounds s, m s ≤ g m)
    (letters : Finset W) (hl : letters ⊆ V) :
    ∑ c ∈ letters, m c ≤ ∑ c ∈ letters, N c :=
  Finset.sum_le_sum fun c hc => certificate V len N m bounds hsound hjust c (hl hc)

end LazySCS

#print axioms LazySCS.rule_sound
#print axioms LazySCS.certificate
#print axioms LazySCS.lower_bound
