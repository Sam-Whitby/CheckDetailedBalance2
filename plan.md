# Implementation Plan: Continuous-Limit VMMC with Gaussian Displacement Proposals

*Date: 25 May 2026*

This document describes the planned changes to `vmmc_continuous.wl` and the
checker infrastructure (`dbc_core.wl`, `check.wls`) to support physically
meaningful continuous-limit VMMC parameterised by a single length scale
`physLen`, with Gaussian displacement proposals that can be checked
symbolically by the existing detailed-balance checker.

---

## 1. Motivation and Summary of Problems with the Current Code

`vmmc_continuous.wl` currently has two interconnected problems.

**Problem 1 — The potential is sub-lattice.**  
`sigLJ = 2^(-1/6) ≈ 0.891` places the LJ zero-crossing at `d² = 0.794`, which
is below the minimum possible particle separation on the lattice (`d² = 1`,
adjacent sites). Every particle pair in the simulation is already in the
attractive tail of the LJ well. There is no repulsive region, no hard-core
exclusion from the energy function, and no resolved LJ diameter. The model is a
nearest-neighbour attractive lattice gas, not an LJ fluid.

**Problem 2 — The direction proposal is an independent free parameter K.**  
Directions are chosen from a fixed list of K unit vectors at equally-spaced
angles, rounded to grid coordinates. K and the lattice size are two independent
parameters that both need to grow to reach the continuous limit. The rounding
also makes most of these directions identical for small step sizes, and
introduces a non-uniform angular distribution that is not fixed by any
physically motivated rule.

The goal of this plan is to replace both with a single physical length scale
`physLen` (the particle diameter in lattice units) that determines the LJ
potential shape, the interaction cutoff, and the width of a Gaussian
displacement proposal. K ceases to be a parameter of the algorithm: it emerges
automatically from the Gaussian's support and the lattice geometry.

---

## 2. Physical Design

### 2.1 The physLen parameter

`physLen` is the particle diameter in lattice units. It is the single physical
input from which all other length scales are derived.

```
sigLJ   = physLen                          (* LJ zero-crossing = particle diameter *)
d_min   = 2^(1/6) * physLen               (* LJ minimum, in lattice units *)
d_min²  = 2^(1/3) * physLen²              (* LJ minimum, in squared lattice units *)
$maxD2  = Ceiling[2 * physLen^2]          (* interaction cutoff: captures the well *)
sigStep = physLen / alpha                  (* Gaussian step std dev; alpha ~ 5 *)
```

For `physLen = 5`, `alpha = 5`: the LJ minimum is at `d ≈ 5.6` lattice units,
`$maxD2 = 50` captures the full well, and `sigStep = 1` lattice unit per step —
the cluster translates by a Brownian-scale increment relative to the particle
diameter.

At adjacent sites (`d² = 1`), the LJ energy with `sigLJ = 5` is approximately
`4 * epsLJ * (5^12 - 5^6) ~ 10^9 epsLJ` — strongly repulsive. This correctly
enforces a soft hard-core at roughly one particle diameter. The lattice
hard-sphere exclusion (no double occupancy) acts as an additional hard wall at
`d = 0`.

As `physLen → ∞` with the lattice size scaled proportionally, the particle
diameter spans many lattice sites and the discrete lattice approaches the
continuum. The dimensionless ratio `physLen / nGrid` (diameter relative to
system size) controls the thermodynamic state point.

### 2.2 The Gaussian displacement proposal

In overdamped Brownian dynamics, a particle's displacement per unit time is
drawn from a Gaussian with standard deviation `√(2 D₀ Δt)`. The natural VMMC
analogue is to draw the cluster translation vector `(dx, dy)` from a 2D
isotropic Gaussian:

```
dx ~ N(0, sigStep),    dy ~ N(0, sigStep)   independently
dir = {Round[dx], Round[dy]}
```

This satisfies the symmetry requirement for superdetailed balance: the proposal
distribution satisfies `p(dx, dy) = p(-dx, -dy)` exactly (Gaussian is symmetric
under sign flip), so the proposal factors cancel in the detailed-balance
condition, and the Whitelam-Geissler link mechanism handles the rest
unchanged.

When `dir = {0, 0}` (the zero-displacement outcome), the move is a no-op and
the state is returned unchanged. This contributes to the self-loop weight
`T(i→i)` and is correctly handled by the checker.

### 2.3 The role of K

K is not set by the user. It emerges from the discretised Gaussian:

- The support of the proposal over distinct integer displacements is
  `{-N_max, ..., N_max}` per dimension, where `N_max = Ceiling[truncSigmas * sigStep]`.
- `truncSigmas` is a computational parameter (default 4; see §3.2 for the
  `Infinity` case).
- The number of distinct `(dx, dy)` pairs with non-zero Gaussian weight is
  `(2 * N_max + 1)²` before applying periodic-boundary collapsing on a finite
  lattice.
- For the checker, on a small lattice of side `nGrid`, the effective distinct
  displacements per dimension are capped at `floor(nGrid / 2)` by the torus
  geometry, so many large displacements collapse to the same final site.

This means K is always determined by `physLen` (via `sigStep`) and the lattice
size being tested. It is no longer an independent design choice.

---

## 3. Changes to `dbc_core.wl`

This is the most significant new work. Two additions are required.

### 3.1 Intercepting `RandomVariate[NormalDistribution[...]]`

Inside the checker's BFS context, `RandomVariate` for named distributions must
be intercepted and converted to the existing `seqBernoulli` mechanism (which
the checker already uses for `RandomChoice[weights → elements]`).

The interception rule:

```mathematica
(* Set inside the checker BFS context, not globally *)
RandomVariate[NormalDistribution[mu_, sigma_]] :=
  With[{
    nMax = If[TrueQ[$dbcNMax === Infinity],
              (* lattice-size cap: see §3.2 *)
              $dbcCurrentNGrid,
              Ceiling[$dbcTruncSigmas * N[sigma]]],
    vals = Range[-nMax, nMax]
  },
  With[{
    rawW = Table[
      CDF[NormalDistribution[mu, sigma], k + 0.5] -
      CDF[NormalDistribution[mu, sigma], k - 0.5],
      {k, vals}]
  },
  seqBernoulli[rawW / Total[rawW], vals]]]
```

Two global variables (`$dbcNMax`, `$dbcTruncSigmas`) control truncation and
are set by `check.wls` from command-line arguments.

The existing `seqBernoulli` BFS machinery handles the branching: for `n` values
in `vals`, the BFS tree generated has `2n − 1` nodes (n leaves, one per
outcome; n − 1 internal branch nodes). This is linear in n, not exponential.

A `Round` pass-through UpValue is also needed so that the natural pattern
`Round[RandomVariate[NormalDistribution[0, s]]]` silently absorbs the `Round`
without modifying the already-integer-valued seqBernoulli:

```mathematica
seqBernoulli /: Round[seqBernoulli[w_, v_]] :=
  seqBernoulli[w, Round /@ v]   (* no-op when v is already integers *)
```

### 3.2 The N_max = Infinity case

When `$dbcNMax = Infinity`, the Gaussian's support is not truncated by a fixed
number of sigmas. Instead, for each lattice being tested, the effective cap is
`floor(nGrid / 2)` — the maximum distinct displacement on a periodic lattice of
side `nGrid`. On a 2×2 lattice this is 1; on a 3×3 lattice this is 1; on a 5×5
lattice this is 2; and so on.

Implementation requires threading the current lattice side length through the
BFS. A global variable `$dbcCurrentNGrid` is set at the start of each
state's BFS traversal:

```mathematica
$dbcCurrentNGrid = Round[Sqrt[Length[currentState]]];
```

This variable is then used inside the `RandomVariate` interception when
`$dbcNMax === Infinity`.

The `Infinity` setting is physically principled: the Gaussian is never
truncated below the lattice's actual resolution. On small lattices (which are
all the checker tests), the effective N_max is tiny (1 or 2), making the BFS
tree small regardless of `physLen`. This is the correct behaviour — the checker
is verifying algorithm structure on simple systems, not simulating large ones.

### 3.3 Extending to other named distributions

The same interception pattern generalises to any named distribution whose CDF
is available in Mathematica. Priority additions:

| Distribution | Use case |
|---|---|
| `UniformDistribution[{a,b}]` | Already handled by contToken; add here for uniform API |
| `ExponentialDistribution[λ]` | Alternative move-size distribution (heavier tail than Gaussian) |
| `CauchyDistribution[x₀, γ]` | Very heavy-tailed; potentially interesting for cluster moves |

The general template:

```mathematica
RandomVariate[dist_] :=
  With[{vals = ..., rawW = Table[CDF[dist, k+0.5] - CDF[dist, k-0.5], {k, vals}]},
  seqBernoulli[rawW / Total[rawW], vals]]
```

The CDF-difference weighting is correct for any unimodal distribution that is
discretised by rounding to the nearest integer.

---

## 4. Changes to `vmmc_continuous.wl`

Section 0 (the user-facing parameter block) is redesigned around `physLen`:

```mathematica
(* ---- SECTION 0: PHYSICAL PARAMETERS ---- *)

physLen = 5          (* particle diameter in lattice units — the single physical scale *)
alpha   = 5          (* step-size divisor: sigStep = physLen / alpha *)
epsLJ   = 1          (* LJ well depth in kT units (exact) *)

(* Derived — do not edit below this line *)
sigLJ   = physLen                     (* LJ diameter *)
sigStep = physLen / alpha             (* Gaussian step std dev in lattice units *)
$maxD2  = Ceiling[2 * physLen^2]     (* interaction cutoff: covers the LJ well *)
```

The concrete coupling function:

```mathematica
$couplingJConcrete[a_Integer, b_Integer, d2_Integer] :=
  If[d2 == 0 || d2 > $maxD2, 0,
     4 * epsLJ * ((sigLJ^2 / d2)^6 - (sigLJ^2 / d2)^3)]
```

With `sigLJ = physLen` (exact integer or rational), all values returned by
`$couplingJConcrete` are exact rational numbers, avoiding the `0.0 ≠ 0` issue
encountered when floating-point coupling values are passed through the symmetry
checker.

The Algorithm function is replaced:

```mathematica
Algorithm[state_List] :=
  Module[{nGrid, occupied, seed, dx, dy, dir, cluster, newState, dest},
    nGrid    = Round[Sqrt[Length[state]]];
    occupied = Flatten[Position[state, _?(# > 0 &)]];
    If[occupied === {}, Return[state]];

    seed = RandomChoice[occupied];

    (* Gaussian displacement proposal — both dx and dy independently *)
    dx = Round[RandomVariate[NormalDistribution[0, sigStep]]];
    dy = Round[RandomVariate[NormalDistribution[0, sigStep]]];
    dir = {dx, dy};

    (* Zero displacement is a no-op *)
    If[dir === {0, 0}, Return[state]];

    cluster = $vmmcBuildCluster[state, nGrid, seed, dir];
    If[cluster === None, Return[state]];

    newState = state;
    Do[newState[[cluster[[i]]]] = 0, {i, Length[cluster]}];
    Do[
      dest = $applyDir[cluster[[i]], dir, nGrid];
      If[newState[[dest]] =!= 0, Return[state, Module]];
      newState[[dest]] = state[[cluster[[i]]]],
      {i, Length[cluster]}];
    newState]
```

The `$dirVectors`, `numDirections`, and `stepSize` variables are removed
entirely. The `$vmmcBuildCluster` and neighbour functions are unchanged — they
already work for arbitrary integer displacement vectors.

---

## 5. Changes to `check.wls`

Two new command-line options:

| Option | Default | Description |
|---|---|---|
| `TruncSigmas=N` | `4` | Gaussian truncation in units of sigStep: N_max = ceil(N * sigStep) |
| `NMax=Infinity` | off | Lattice-size cap instead of sigma-based cap (see §3.2) |

These are read from `kvArgs` and stored as `$dbcTruncSigmas` and `$dbcNMax`
before the BFS begins. Both are made available to `dbc_core.wl` as globals
prior to any Algorithm execution.

The `NMax=Infinity` flag is passed as the string `"Infinity"` and converted via
`ToExpression` to the Mathematica symbol `Infinity`.

No changes to the BFS logic, the symbolic DB check, or the ergodicity check
are required.

---

## 6. What the Checker Now Does

The checker's behaviour after these changes, stated cleanly:

1. **Loads the algorithm file** and reads `physLen`, `sigStep`, `$maxD2`,
   `epsLJ`, `sigLJ` from it. These are concrete numbers. The abstract
   `couplingJ` function has no DownValues during the symbolic check.

2. **Intercepts `RandomVariate[NormalDistribution[mu, sigma]]`** within
   Algorithm execution and replaces each call with a `seqBernoulli` over
   integer values weighted by Gaussian CDF differences. The support is
   `{-N_max, ..., N_max}` where N_max is set by `TruncSigmas` or capped by
   the current lattice geometry if `NMax=Infinity`.

3. **Runs BFS** over all possible combinations of random outcomes: seed
   particle choice, (dx, dy) from the Gaussian seqBernoulli, and binary link
   decisions in the cluster builder.

4. **Builds the symbolic transition matrix** T(i→j) as before. Each entry is a
   rational combination of Gaussian weights (concrete numbers) and abstract
   `couplingJ[a, b, d2]` atoms.

5. **Checks detailed balance** symbolically for each state pair: verifies that
   `T(i→j) * exp(-β E(i)) - T(j→i) * exp(-β E(j)) = 0` after FullSimplify
   with `couplingJ` atoms treated as free reals and β > 0.

6. **Checks ergodicity** by comparing the number of BFS-reachable states to the
   falling factorial (S)_N.

The checker never evaluates `couplingJ` at specific distances. It does not know
or care what `physLen` is numerically — the energy algebra is purely symbolic.
The only things that change numerically with `physLen` are the Gaussian weights
in the seqBernoulli and the interaction cutoff `$maxD2` (which is irrelevant on
small lattices where all pairs are within any reasonable cutoff).

---

## 7. What Algorithms Must Provide

The checker contract does not substantially change. Algorithm files must define:

```mathematica
energy[state_]       (* pair energy as before, using couplingJ[a, b, d2] *)
Algorithm[state_]    (* MCMC move; may now use RandomVariate[NormalDistribution[...]] *)
BitsToState[bits_]   (* unchanged *)
numBeta              (* unchanged *)
```

For continuous VMMC, the displacement proposal must satisfy:
- `p(dx, dy) = p(-dx, -dy)` — the Gaussian satisfies this automatically
- The reverse of any proposed direction `{dx, dy}` must also be in the support —
  trivially true for a symmetric Gaussian
- The cluster builder must be able to handle any integer displacement vector `{dx, dy}`

Algorithms are not required to use `RandomVariate`. The existing mechanisms
(`RandomInteger`, `RandomChoice`, `RandomReal` in comparisons) remain valid and
unchanged. `RandomVariate` support is additive, not a replacement.

---

## 8. Uncertainties and Difficulties

### 8.1 Symbolic simplification of Gaussian weight ratios

The most significant uncertainty. The DB condition for a VMMC move with
direction `(dx, dy)` and its reverse `(-dx, -dy)` requires:

```
T(i→j; dir=(dx,dy)) * gaussW[dx] * gaussW[dy]  =
T(j→i; dir=(-dx,-dy)) * gaussW[-dx] * gaussW[-dy]
```

Since `gaussW[k] = gaussW[-k]` (Gaussian symmetry), the weight factors cancel.
But in the symbolic checker, `gaussW[dx]` and `gaussW[-dx]` are **concrete
numbers** (evaluated from the CDF formula with specific `sigStep`). Their
equality is numerical, not symbolic.

If `sigStep` is an exact rational (e.g., `physLen = 5`, `alpha = 5`, so
`sigStep = 1`), the Gaussian CDF values are exact symbolic expressions involving
`Erf`. Whether `FullSimplify` recognises that `Erf[x]` terms cancel across
forward and reverse paths is not guaranteed without careful testing.

**Possible mitigation:** leave `sigStep` as a symbolic atom (add it to
`DynamicSymParams` as an abstract parameter assumed positive). Then
`gaussW[k, sigStep]` and `gaussW[-k, sigStep]` would be equal by the symmetry
property `gaussW[-k, s] = gaussW[k, s]` which could be added as an assumption
to FullSimplify. This requires extending the `symParams`/`DynamicSymParams`
mechanism.

**Alternative mitigation:** provide a dedicated FastChecker path for the
Gaussian case that explicitly cancels the symmetric weight pairs before passing
expressions to Simplify.

This needs careful prototyping before committing to the design.

### 8.2 Threading $dbcCurrentNGrid through the BFS

When `NMax = Infinity`, the `RandomVariate` interception needs to know the
current lattice side length to cap the seqBernoulli support. This requires a
global variable `$dbcCurrentNGrid` that is set before each state's BFS
traversal. The BFS in `dbc_core.wl` iterates over states; adding the global
assignment is a small change but needs to be done carefully to avoid
contaminating parallel kernel contexts (if `ParallelMap` is active for the
symbolic check).

In practice, since the BFS itself is sequential (only the subsequent symbolic
simplification is parallelised), this is straightforward.

### 8.3 BFS tree size growth with physLen

For large `physLen` and `alpha = 5`, `sigStep = physLen/5`, `N_max ≈ 4 * physLen/5`.
The seqBernoulli for each of dx and dy has `2*N_max+1` leaves. The BFS tree for
the displacement alone has approximately `(2*N_max+1)²` leaf combinations:

| physLen | N_max | Displacement BFS leaves |
|---|---|---|
| 5 | 4 | 81 |
| 10 | 8 | 289 |
| 20 | 16 | 1,089 |
| 50 | 40 | 6,561 |

This scales as O(physLen²), remaining linear (not exponential) and tractable
for the small lattices the checker tests. However, at large `physLen` (>50),
the checker becomes noticeably slower even on 2×2 lattices. A warning should be
printed if `N_max > 20`.

For `NMax = Infinity`, the effective N_max on a 2×2 lattice is 1 (3 values per
dimension, 9 total displacement combinations), which is trivially fast regardless
of physLen. This makes `NMax = Infinity` the recommended default for checker
runs.

### 8.4 The zero-displacement outcome

When `dx = 0` and `dy = 0` (which occurs with probability `gaussW[0]² > 0`),
the algorithm returns the current state unchanged. This is a self-loop
contribution to `T(i→i)`. The checker handles self-loops correctly (they appear
in the diagonal of the transition matrix). However, if `gaussW[0]` is large
(which happens when `sigStep << 1`), almost all BFS paths lead to self-loops,
and the off-diagonal entries of T become very small. The DB check still passes,
but the checker exercises very little of the VMMC logic. A warning should be
emitted if the (0,0) probability exceeds 50%.

### 8.5 $maxD2 vs small-lattice geometry

For `physLen = 5`, `$maxD2 = 50`. On a 3×3 lattice (maximum d² = 4), this is
irrelevant — all bonds are within the cutoff. The `$uniqueBondsExt` function
already handles this correctly (it selects bonds by distance, not by cutoff
index). However, computing `$uniqueBondsExt[nGrid]` for large `$maxD2` and
small `nGrid` is slightly wasteful. The existing `d2 ≤ $maxD2` filter in
`$uniqueBondsExt` means the result is the same as `$maxD2 = max_d2_on_lattice`,
just with an unnecessary comparison. No correctness issue, minor efficiency issue.

### 8.6 sigLJ as an exact integer

Setting `sigLJ = physLen` as an exact integer (not a float) is necessary to
avoid the `0.0 ≠ 0` bug in the symmetry checker that was encountered previously.
With `physLen = 5` (exact), `$couplingJConcrete[1, 2, 1] = 4*(25^6 - 25^3)` is
an exact integer. The symmetry check passes cleanly. If a user sets `physLen =
5.0` (a float), this breaks. The plan should include a guard:

```mathematica
If[!IntegerQ[physLen], physLen = Round[physLen];
   Print["Warning: physLen rounded to nearest integer: ", physLen]]
```

### 8.7 alpha as a free parameter

The divisor `alpha` (determining `sigStep = physLen / alpha`) is left as a user
parameter in Section 0. There is no universally correct value: `alpha = 5`
(steps of 1/5 particle diameter) is appropriate for moderate temperatures and
densities, but optimal alpha depends on acceptance rates and the intended
physical timescale. This is documented but not automated.

### 8.8 Backward compatibility

Existing algorithm files (`kawasaki_1d.wl`, `vmmc_2d.wl`, etc.) use no
`RandomVariate` calls. The new interception rules are additive and will not
affect them. The `vmmc_continuous.wl` rewrite replaces the current version; no
other examples are affected.

---

## 9. Implementation Order

The changes should be made in the following order to allow incremental testing:

1. **`dbc_core.wl`**: Add `RandomVariate[NormalDistribution[...]]` interception
   with a fixed `$dbcNMax` global. Test that the seqBernoulli is generated
   correctly for a simple test case (e.g., `sigStep = 1`, compare output weights
   to known Gaussian CDF differences).

2. **`check.wls`**: Add `TruncSigmas` and `NMax` argument parsing. Set globals
   before BFS. Test that a minimal algorithm using `RandomVariate` can be loaded
   and checked.

3. **`vmmc_continuous.wl`**: Replace Section 0 with the `physLen`/`alpha` design.
   Replace the `Algorithm` function with the Gaussian displacement version.
   Verify symmetry checker passes with `physLen = 1` (the near-previous behaviour)
   and `physLen = 3` (a genuinely physical case visible on a 3×3 lattice).

4. **Symbolic cancellation testing**: Run the full symbolic check on the 3×3
   seed with `physLen = 3` (LJ minimum at `d ≈ 3.4` lattice units, visible on
   the lattice). Verify that the Gaussian weight ratios cancel in the DB
   expressions. If they do not, implement the symbolic `sigStep` atom fallback
   described in §8.1.

5. **`NMax = Infinity` support**: Add `$dbcCurrentNGrid` threading through the
   BFS. Test that on a 2×2 lattice, the effective N_max is 1 regardless of
   physLen.

6. **Extend to other distributions**: Add `ExponentialDistribution` interception
   once Normal is confirmed working. Test with a simple 1D jump algorithm using
   exponential step sizes.

7. **README and documentation**: Update examples table, add `physLen` to the
   Animation section, document `TruncSigmas` and `NMax` options.
