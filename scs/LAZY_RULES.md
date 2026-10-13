# Lazy periodic rules keep the factor-2 guarantee

This note justifies the one place where `scs-rs` / `two_approx_c_lazy` departs from
the paper ("A Polynomial-Time 2-Approximation for Shortest Common Superstring",
openai/math result 128). Numbers like (2.6) and Lemma 3.3 refer to that paper.

## The change

The paper fixes the counts `m` by the extension inequalities (2.5), the
required-string bound, and **every** periodic rule (2.6). Enumerating all rules
is the expensive part (one per word with a border, i.e. almost every word).

The lazy variant computes counts `m'` from:

1. the extension inequalities (2.5) and the required-string bound, exactly;
2. periodic rules (2.6) only for an explicit set `E` of enabled periodic texts
   `(p, least rotation)`, enforced by iterating to a fixpoint (each pass
   evaluates every enabled rule on the current counts, raises `m'(w)` where a
   triggered rule is violated, then recomputes the extension closure).

It then runs Sections 2.3-7 unchanged. If Section 5.3 finds no record, it
enables the text of the group being processed, recomputes `m'`, and restarts.

## Claim

On every input the lazy variant terminates in polynomial time and outputs a
common superstring of length at most `2W' <= 2 OPT`, where `W' = sum_c m'(c)`.

## Argument

**Step 1: every value is a valid lower bound.** The paper's proof of Lemma 2.1
needs, at word `w`, only that every count it reads (of longer words) is at
most the true occurrence count `N_T` in an arbitrary common superstring `T`.
That holds for any vector `m' <= N_T`, not only for the paper's `m`:

- The extension bounds and the required-string bound are immediate.
- For a rule with trigger `m'(v) > B_A(v; m')`, the two cases of the proof go
  through with `m'` in place of `m`. Both cases use only that `B_A(.; .)` has
  nonnegative coefficients and that `N_T - m' >= 0` termwise. If
  `N(v) > B_A(v; N)`, then `N(w) >= B_A(w; N) + k + 1 >= B_A(w; m') + k + 1`.
  Otherwise some bracketing word `r` has `N(r) - m'(r) >= 1` and
  `b_{A,w}(r) >= k + 1`, so `B_A(w; N) - B_A(w; m') >= k + 1`.

Induct over the sequence of assignments, not over word length. Every value ever
assigned is the maximum of bounds computed from values that are already
`<= N_T`, so it is itself `<= N_T`. Hence `m' <= N_T` and `W' <= OPT`. Values
only increase and are bounded by `N_T <= L`, so the fixpoint terminates.

**Step 2: the construction never needs the other rules.** Sections 2.3-7 use
the counts in three ways:

- **Balance and coverage (Lemma 2.2)** need only (2.5) and the required-string
  bound. The final pass of the fixpoint recomputes the extension closure, so
  `m'` satisfies them exactly. Multiplicities `u, d >= 0`, the base graph is
  balanced with `W'` up edges, and every required string is a used vertex.
- **Lemmas 3.1-3.5** are consequences of the edge multiplicities and of the
  definition of `B_A`. They hold for whatever counts defined the graph. In
  particular `m'(s) - B_A(s; m')` counts unblocked recorded occurrences (3.4),
  and the top layer in Section 5.3 records an unblocked occurrence of
  `v = A[t:H)`, so the trigger of the rule for `(A, v, w = A[t:H-K), k)` holds
  automatically.
- **Section 5.3 is the only place a periodic rule is invoked.** It is used to
  show that a record exists: an actual window on a layer with a different
  aligned text, agreeing with `A` and containing `[t, H-K)`. Every later step
  (5.4, Lemmas 4.3-4.4, Section 6) uses only the properties of the record that
  was found: `a <= t`, `b >= H - K`, a different alignment, and agreement with
  `A`. The bound `b - a < p + q` follows from Lemma 3.4 alone.

So the lazy variant needs nothing beyond what it checks at run time. If the
search succeeds, the remaining steps are justified exactly as in the paper.
If it fails, the rule for this text is missing. After enabling the text, the
fixpoint makes every triggered rule on that text hold. If that same text is ever
needed again, Lemma 5.1's argument applies and the search cannot fail (the
implementation asserts this; the assertion has never fired).

**Step 3: the budget.** Layer budgets are periods, and by Lemma 3.2 they sum to
the up cost of the base graph, which is `W'`. The connection procedure spends at
most the budgets (Theorem 6.3), so the Euler tour has at most `2W'` letters.

**Termination and running time.** Each restart enables a new text, and the
implementation asserts it is new. Every text that can carry a rule is a rotation
class of a length-`p` substring of an input, so there are at most `O(L^2)` of
them. Each round is polynomial, as in the paper.

## Consequences and caveats

- `W' <= W <= OPT`. The output bound `2W'` is never weaker than the paper's
  `2W`. The per-instance optimality certificate `len / W'` can be looser than
  `len / W` when rules would have raised `W`. On read data they almost never
  do: 0 lazy rules were needed on all read sets tested, including real E. coli.
- This is a proof sketch by the implementer, not a peer-reviewed proof. It
  rests on re-reading the paper's Lemma 2.1 and Section 5.3; anyone relying on
  it should check those two places.
- Empirical support: 4,503 fuzz instances (Python and Rust) all satisfy
  `len <= 2W'` and `W' <= OPT` (exact optimum via subset DP). They needed at
  most 7 rounds, and no rule-enabled text was ever requested twice.
