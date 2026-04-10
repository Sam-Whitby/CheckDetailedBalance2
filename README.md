# DetailedBalanceChecker

A symbolic exhaustive checker for detailed balance in Monte Carlo algorithms.
Given a Markov-chain algorithm written in Mathematica, the checker proves — or disproves — that it satisfies detailed balance exactly, using computer algebra (FullSimplify with β > 0) rather than statistical tests alone.

---

## How the checker works

### 1. Reducing randomness to fair coin flips

Every source of randomness in the algorithm is reduced to reads from a binary tape:

- **`RandomInteger[{lo, hi}]`** and **`RandomChoice[list]`**: use rejection sampling over `k = ⌈log₂(n)⌉` bits to select an integer uniformly from `{0, …, n−1}`. Paths that fall outside the range are silently discarded (the missing probability fraction is uniform across all starting states, so the unnormalised transition matrix still satisfies detailed balance).

- **`RandomReal[]`**: uses *interval tracking* (Approach 3, see below).

### 2. Interval tracking for `RandomReal[]`

Each call to `RandomReal[]` creates a fresh latent variable `U ~ Uniform[0,1]`, represented by a token `$dbc$irand[j, acceptTestI]` that records the current belief interval `[lo, hi]` for that variable (initially `{0, 1}`).

Whenever the token is compared to a threshold `p`:

```
P(U < p | U ∈ [lo, hi]) = (p − lo) / (hi − lo)
```

One bit is read from the tape. If the bit is 1 (accept), the weight is multiplied by `(p − lo)/(hi − lo)` and the interval narrows to `[lo, p]`. If the bit is 0 (reject), the weight is multiplied by `(hi − p)/(hi − lo)` and the interval narrows to `[p, hi]`. Subsequent comparisons of the **same** variable correctly condition on all previous narrowings.

This is equivalent to inverse-CDF sampling: the bit tape encodes the quantile of `U` at sufficient resolution for each comparison that the algorithm makes.

**Deterministic short-circuit (no bit consumed):** Before branching, the checker attempts a lightweight `PiecewiseExpand` on the threshold `p`. If the result is ≤ `lo` the comparison is certainly False; if ≥ `hi` it is certainly True. For a fresh variable with interval `{0, 1}` this means thresholds that evaluate to exactly 0 or exactly 1 are resolved without reading any bit and without creating a BFS branch. This handles common physical cases automatically — for example, a zero link-formation probability when the move does not change a bond energy (`wFwd = 0`), or a mandatory link when two particles would collide (`wFwd = 1` from hard-sphere exclusion) — even if the algorithm writes these using `Max`, `Piecewise`, or any other form, without requiring the algorithm author to add special-case guards.

**Why this matters**: In the previous (Approach 2) formulation, each `RandomReal[]` call was modelled as a single independent Bernoulli trial with probability `p` — correct only when the token is compared exactly once. Interval tracking is correct for any number of comparisons on the same variable and enables the next feature.

### 3. `RandomChoice[weights -> elements]`

Weighted choice is decomposed into `n − 1` independent sequential Bernoulli trials, each using a fresh `RandomReal[]` interval-tracking variable. For weights `w₁, …, wₙ` with sum `W`:

- Trial 1: accept element 1 with probability `w₁/W`.
- Trial 2 (if rejected): accept element 2 with probability `w₂/(W − w₁)`.
- …
- Trial `n−1` (if all rejected): default to element `n`.

This telescopes to give `P(element k) = wₖ/W` exactly, and all probabilities are symbolic expressions that flow through to `FullSimplify`.

### 4. BFS over bit sequences

Starting from a seed state, the checker performs BFS over all possible bit sequences. Each leaf of the BFS tree corresponds to one complete algorithm execution, characterised by:

- the bit sequence consumed,
- the resulting next state,
- the path weight (product of all per-bit factors).

All leaves that terminate in the same `(start state, end state)` pair are summed to give the symbolic transition matrix `T[i → j]`.

### 5. Symbolic detailed-balance check

For each pair of distinct states `(i, j)`, the checker verifies:

```
T(i→j) · exp(−β E(i)) − T(j→i) · exp(−β E(j)) = 0
```

A **tiered simplification** strategy is used, escalating only when cheaper methods fail:

| Tier | Method | Catches |
|------|--------|---------|
| 1 | Literal `=== 0` test | Both T entries are zero (disconnected state pairs) |
| 2 | `Simplify[..., β > 0]` | Algebraic cancellation without expanding Piecewise cases |
| 3 | `FullSimplify[PiecewiseExpand[...], β > 0]` | Piecewise Metropolis terms requiring sign-case analysis |

All three tiers are exact — no numerical approximations are made. The Boltzmann factor `β` is kept as a free symbol throughout. Only state pairs that survive Tier 1 and Tier 2 reach the expensive `FullSimplify + PiecewiseExpand` call, which significantly reduces runtime for systems with many disconnected or simply-balanced state pairs.

### 6. Numerical MCMC validation

As an independent check, the algorithm is run as a genuine Markov chain for `NSteps` steps. The empirical state distribution is compared to the analytical Boltzmann distribution via the KL divergence; `KL < 0.02` is treated as a pass.

---

## Usage

```
wolframscript -file check.wls <algorithm_file.wl> [options]
```

Options (no spaces around `=`):

| Option | Default | Description |
|--------|---------|-------------|
| `MaxBitString=XXXX` | `11111111` | Largest bit string to test. All bit strings from length 1 through `len(XXXX)` are generated, stopping at `XXXX` within the final length. |
| `Mode=Symbolic\|Numerical\|Both` | `Both` | Run symbolic check only, numerical MCMC only, or both. |
| `NSteps=N` | `50000` | MCMC steps per connected component. |
| `MaxBitDepth=N` | `20` | BFS depth cap per state. |
| `Verbose=True\|False` | `False` | Print per-state BFS progress. |

### Algorithm file format

Each `.wl` file must define:

```mathematica
energy[state_]       (* bare energy, no β factor, may contain symbolic J symbols *)
Algorithm[state_]    (* MCMC move; use RandomReal[], RandomInteger[], RandomChoice[] *)
BitsToState[bits_]   (* {0,1,...} list → seed state, or None to skip *)
numBeta              (* numeric inverse temperature for MCMC check *)
```

Optional:

```mathematica
symParams = <|"eps" -> {...}, "couplings" -> {...}|>   (* static symbolic parameters *)
DynamicSymParams[states_List] := ...                   (* per-component symbolic parameters *)
DisplayState[state_]  := ...                           (* custom state display string *)
ValidStateIDs[maxId_] := ...                           (* fast enumeration of valid seeds *)
```

For Metropolis acceptance, use exactly:
```mathematica
If[RandomReal[] < MetropolisProb[dE], newState, state]
```
where `dE = energy[newState] - energy[state]`. `MetropolisProb` is intercepted symbolically.

See `template.wl` for a fully annotated starting point.

---

## Supported random calls

| Call | How it is intercepted |
|------|----------------------|
| `RandomReal[]` | Interval-tracking token; each comparison reads one bit and narrows the interval |
| `Random[]` | Same as `RandomReal[]` |
| `RandomInteger[{lo, hi}]` | Rejection sampling over `⌈log₂(hi−lo+1)⌉` bits |
| `RandomInteger[n]` | Uniform over `{0, …, n}` via rejection sampling |
| `RandomChoice[list]` | Uniform choice via rejection sampling over `⌈log₂(n)⌉` bits |
| `RandomChoice[weights -> elements]` | Sequential Bernoulli decomposition using independent `RandomReal[]` variables |

**Not supported**: `RandomVariate`, `RandomSample`, `RandomPermutation`, `AbsoluteTime`, and similar calls will cause the analysis to abort with an error message.

---

## Example files (`examples3/`)

### Expected to PASS

| File | Description |
|------|-------------|
| `kawasaki_1d.wl` | 1D Kawasaki (nearest-neighbour swap) on a periodic ring. Labeled particles, symmetric proposal via `RandomInteger[{0, L-1}]` for bond selection, Metropolis acceptance. |
| `kawasaki_2d.wl` | 2D Kawasaki on a periodic square lattice. Uses `RandomInteger` for bond selection and `RandomReal[]` for Metropolis acceptance. |
| `vmmc_2d.wl` | Virtual Move Monte Carlo (Whitelam–Geissler) cluster moves on a 2D lattice. Uses `RandomChoice[occupied]`, `RandomChoice[dirs]`, and multiple symbolic `RandomReal[]` comparisons against link-formation probabilities. |
| `vmmc_2d_thesis.wl` | **Chapter 2 column VMMC**: extends `vmmc_2d.wl` to the condensate-column geometry from Sam Whitby's PhD thesis (Chapter2/src/VMMC.cpp + NucleolusModel.cpp). Hard wall in x, periodic in y; interactions at d=1,√2,2,√5 all sharing a single coupling symbol `J<a><b>` per type pair; optional chemical gradient γ(x_mid)=x_mid/L; optional cluster-size cutoff, SL moves, and Stokes drag; internal-bond Metropolis filter for DB correctness with gradient. **DB check result: PASS** — symbolic + numerical (KL < 0.01) for all states up to `MaxBitString=111111` with default settings (`$tvHardBarrier=True`, `$tvGradient=False`). |
| `jump_1d_edit.wl` | 1D jump dynamics: a particle jumps by a displacement `d` drawn via **rejection sampling** (uniform `RandomInteger[{0,L-1}]`, then thinned by `RandomReal[] >= w(d)/w(0)`) to implement a discrete normal distribution. Symmetric proposal, Metropolis acceptance. |
| `jump_1d_weighted.wl` | Same jump dynamics as `jump_1d_edit.wl`, but the displacement is drawn **directly** via `RandomChoice[normalWeights -> displacements]`. This is the canonical test of the new weighted-`RandomChoice` support (Approach 3). Both implementations are mathematically equivalent and both pass detailed balance. |

### Expected to FAIL

| File | Bug | Why it fails |
|------|-----|-------------|
| `kawasaki_1d_fail.wl` | `dE = energy[state] - energy[newState]` (sign reversed) | Metropolis acceptance penalises downhill moves and always accepts uphill moves, inverting the Boltzmann distribution. Violations appear at `MaxBitString ≥ 1111101` (4-site ring, 2 labeled particles). |
| `kawasaki_2d_fail.wl` | Same wrong-sign `dE` in 2D | Same inversion of Boltzmann distribution. |
| `cluster_1d_fail.wl` | Cluster slides only rightward | Forward and reverse proposals are not symmetric: sliding right has probability 1 but sliding left has probability 0 for the same bond, breaking detailed balance. |
| `cluster_2d.wl` | Cluster merging changes the cluster count | The number of clusters in the system is not conserved, so the forward and reverse move probabilities do not match. |

---

---

## `vmmc_2d_thesis.wl` — design notes and how to run

### What was implemented

`vmmc_2d_thesis.wl` is a Mathematica translation of the Chapter 2 C++ VMMC
(`Chapter2/src/VMMC.cpp` + `NucleolusModel.cpp`).  It shares the same
bijective integer encoding, `BitsToState`, `DisplayState`, and
`ValidStateIDs` as `vmmc_2d.wl`, but replaces the Algorithm and energy
function to match the column physics:

| Feature | `vmmc_2d.wl` | `vmmc_2d_thesis.wl` |
|---------|-------------|---------------------|
| Boundary conditions | fully periodic (torus) | periodic in y (rows), hard wall in x (cols) |
| Interaction range | d=1 (nearest neighbour) | d=1, √2, 2, √5 |
| Coupling symbols | `J_ab` (one per type pair) | `J_ab` per type pair, applied at all four distance ranges |
| Chemical gradient | none | optional `γ(x_mid)=x_mid/L` |
| Cluster-size cutoff | none | optional `⌊1/r⌋` (C++ proposeMove) |
| Saturated-link moves | none | optional, probability `$tvProbSL` |
| Internal-bond Metropolis | no-op | applied when gradient active |
| Neighbour shell for link tests | union of shells of p, pPost, pRev at d=1 | same, at d=√5 |

The hard-barrier mode (`$tvHardBarrier = True`) closes the system so that
detailed balance can be checked — no exit/replacement dynamics occur.

The internal-bond Metropolis filter (C++ lines 758–817) corrects for the fact
that VMMC link weights only Boltzmann-weight *boundary* bonds (one particle
moving, one static).  When interactions are spatially varying (gradient on),
*internal* bonds (both particles moving together) also change energy; the
filter applies `min(1, exp(−β ΔE_int))` to compensate.  For uniform coupling
(`$tvGradient=False`) every internal ΔE is zero and the filter is a no-op.

### How to run the detailed-balance check

Default (uniform coupling, no gradient, no cutoff, hard barrier) — **PASS** confirmed:

```bash
cd CheckDetailedBalance2
wolframscript -file check.wls examples3/vmmc_2d_thesis.wl \
  MaxBitString=111111 Mode=Both NSteps=20000
```

To view the full decision tree and transition matrix for a specific state (e.g. the 3-particle, 4-site system):

```bash
wolframscript -file report.wls examples3/vmmc_2d_thesis.wl BitString=101001
```

With gradient (tests the internal-bond Metropolis filter):

```bash
wolframscript -file check.wls examples3/vmmc_2d_thesis.wl \
  MaxBitString=111111 Mode=Both NSteps=20000
# Then edit the file to set $tvGradient = True before running.
```

With cluster-size cutoff enabled (adds BFS branches; slower):

```bash
# Edit: $tvCutoff = True
wolframscript -file check.wls examples3/vmmc_2d_thesis.wl \
  MaxBitString=11111 Mode=Symbolic
```

Animation (square grid, e.g. 3×3 with 4 particles):

```bash
wolframscript -file animate.wls examples3/vmmc_2d_thesis.wl \
  Sites=9 N=4 Steps=500 FPS=10 Jd112=−1.0 Jdsq212=−0.5
```

### Known issues / things to verify

1. **DB check runtime**: The checker's BFS over bit strings grows quickly with
   extended interactions (four coupling symbols per type pair instead of one).
   Start with `MaxBitString=1111111` (7 bits, states up to 2-site grids) and
   increase gradually.  For 3-site grids (`MaxBitString=11111111`) the check
   may take tens of minutes.

2. **Cutoff with DB checker**: The cluster-size cutoff (`$tvCutoff=True`)
   adds a `RandomReal[]` draw whose comparison `r ≤ 1/k` creates many BFS
   branches.  Leave `$tvCutoff=False` (default) for DB checking; the cutoff
   does not affect detailed balance theoretically.

3. **Gradient check**: With `$tvGradient=True` and a fixed L, the γ factors
   become rational numbers (e.g. γ=3/4 for the mid-column bond in a 4-column
   grid).  The DB checker should still handle these symbolically.  The internal-
   bond Metropolis filter is required for the check to pass.

4. **C++ neighbour finding**: The C++ `recursiveMoveAssignment` queries
   `interactionsCallback` at the particle's *current* position only, whereas
   `vmmc_2d_thesis.wl` (like `vmmc_2d.wl`) uses the union of neighbour shells
   of p, pPost, and pRev.  The union is required for correctness at extended
   range: a particle within √5 of pPost but outside √5 of p changes energy
   upon the move and must be link-tested.  If the C++ interaction range does
   not cover pPost-neighbours, the C++ code may miss some link tests.

---

## Other scripts

| Script | Usage | Description |
|--------|-------|-------------|
| `report.wls` | `wolframscript -file report.wls <alg_file.wl> BitString=<bits>` | Runs the full check on one specific bit string and opens an interactive graphical report window (decision trees, transition matrix, detailed-balance table, frequency bar chart). |
| `animate.wls` | `wolframscript -file animate.wls <alg_file.wl> [options]` | Runs the algorithm as a live Markov chain and produces an animation of the evolving state. |
| `show_report.py` | Called by `check.wls` when no Mathematica frontend is available | Renders the JSON report produced by the checker as a static PNG image. |
| `animate_plot.py` | Called by `animate.wls` | Renders animation frames as a video or GIF. |
| `template.wl` | — | Annotated starting template for new algorithm files. |
