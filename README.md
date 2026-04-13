# DetailedBalanceChecker

Symbolically proves or disproves detailed balance for Monte Carlo algorithms written in Mathematica. Uses computer algebra rather than statistical tests alone.

---

## How it works

1. **Randomness → bits.** Every random call is intercepted and reduced to reads from a binary tape:
   - `RandomInteger[{lo,hi}]` / `RandomChoice[list]` — rejection sampling over ⌈log₂(n)⌉ bits.
   - `RandomReal[]` — *interval tracking*: a latent variable U ∈ [lo,hi]; each comparison `U < p` reads one bit and narrows the interval to [lo,p] or [p,hi], with weight (p−lo)/(hi−lo) or (hi−p)/(hi−lo) respectively.
   - `RandomChoice[weights → elements]` — decomposed into n−1 sequential Bernoulli trials.

2. **BFS over bit sequences.** Starting from a seed state, all possible bit sequences are enumerated. Each complete path gives a (start, end, weight) triple. Weights for the same (start, end) pair are summed to form the symbolic transition matrix T.

3. **Symbolic detailed-balance check.** For each state pair (i, j):
   ```
   T(i→j)·exp(−β E(i)) − T(j→i)·exp(−β E(j))  =?  0
   ```
   β is kept as a free symbol so Boltzmann factors cancel algebraically. Two checkers are available (see below).

4. **Numerical MCMC check.** The algorithm is run as a genuine Markov chain; the empirical distribution is compared to the Boltzmann distribution via KL divergence. KL < 0.02 is a pass.

---

## Symbolic checkers

Two symbolic checking methods are available, selectable per run:

### Standard checker (default)
Uses `FullSimplify[PiecewiseExpand[expr], {β > 0, ...}]` for each state pair. Correct for all algorithm types. Slow for complex algorithms (tens of seconds per component).

### FastChecker (`FastChecker=1`)
An Exp-polynomial algebraic checker. After `PiecewiseExpand`, detailed balance expressions for Boltzmann algorithms reduce — within each feasible sign case — to sums of terms `c·exp(−β·L)` where c is a rational coefficient and L is a polynomial in coupling constants. Detailed balance holds iff all coefficient sums vanish, which is verified by `Expand[...] === 0` (microseconds).

**Algorithm:**
1. `PiecewiseExpand` splits the expression into sign cases (e.g. `Jpair12 > 0`, `Jpair13 == 0`).
2. Each case is tested for feasibility (quick `Reduce` on linear real arithmetic); infeasible cases are skipped.
3. Equality constraints in the condition (e.g. `Jpair12 == 0`) are substituted into the value.
4. The value is normalised: `(E^a)^n → E^(n·a)`, products of `Exp` are merged, and the result is expanded to a sum of `c·E^(poly)` terms.
5. Terms are grouped by their exponent polynomial (compared via `Expand[p1−p2] === 0`). If all coefficient groups sum to zero, the expression is identically zero.
6. **Fallback:** if any step is inconclusive (unusual expression structure), the pair falls back to `FullSimplify`. The result is always exact — never approximate.

The fast checker is parallelised with the same subkernel pool as the standard checker.

**Measured speedups** (on a 4-core machine, `MaxBitString=1111111111`, `Mode=Symbolic`):

| Algorithm | Standard | FastChecker | Speedup |
|-----------|----------|-------------|---------|
| `vmmc_2d.wl` | ~66s | ~49s | 26% |
| `vmmc_2d_field.wl` | ~130s | ~61s | 53% |

The speedup is larger for `vmmc_2d_field.wl` because its field and coupling terms produce complex expressions that `FullSimplify` struggles with, but which the polynomial coefficient check resolves trivially.

**Applicability:** Works directly for all algorithms using Metropolis acceptance. Falls back gracefully for other acceptance functions (Barker, heat-bath) or unusual expression structures.

---

## Usage

### Checker
```bash
wolframscript -file check.wls <algorithm.wl> [options]
```

| Option | Default | Description |
|--------|---------|-------------|
| `MaxBitString=XXXX` | `11111111` | Largest bit string tested |
| `Mode=Symbolic\|Numerical\|Both` | `Both` | Which checks to run |
| `NSteps=N` | `50000` | MCMC steps for numerical check |
| `MaxBitDepth=N` | `20` | BFS depth cap per state |
| `FastChecker=1` | off | Use Exp-polynomial fast checker (see above) |
| `Verbose=True` | `False` | Print per-state BFS progress |

### Report (single seed state)
```bash
wolframscript -file report.wls <algorithm.wl> BitString=XXXXX [options]
```

| Option | Default | Description |
|--------|---------|-------------|
| `BitString=XXXXX` | *(required)* | Seed bit string, e.g. `11001` |
| `NSteps=N` | `50000` | MCMC steps |
| `MaxBitDepth=N` | `20` | BFS depth cap |
| `ShowTrees=True\|False` | `True` | Include decision trees in report |
| `ShowT=True\|False` | `True` | Include transition matrix in report |
| `DoNumerical=True\|False` | `True` | Run numerical MCMC and include scatter plot |
| `PrintMatrix=True\|False` | `False` | Print transition matrix entries as LaTeX |

### Animation
```bash
wolframscript -file animate.wls <algorithm.wl> Sites=<n> N=<n> [options]
```

| Option | Default | Description |
|--------|---------|-------------|
| `Sites=N` | *(required)* | Total lattice sites |
| `N=N` | *(required)* | Number of labeled particle types |
| `Steps=N` | `200` | MCMC steps to run |
| `Beta=f` | from `.wl` file | Inverse temperature |
| `FPS=f` | `10` | Animation frame rate |
| `Simple=1` | off | Fast 2-colour mode (holes vs particles) |
| `RecordEvery=N` | `1` | Record state every N steps |
| `Jpair<a><b>=f` | random | Coupling constant, e.g. `Jpair12=-1.0` |

---

## Algorithm file format

Each `.wl` file must define:

```mathematica
energy[state_]       (* bare energy, no β factor *)
Algorithm[state_]    (* MCMC move; use RandomReal[], RandomInteger[], RandomChoice[] *)
BitsToState[bits_]   (* {0,1,...} list → seed state, or None if invalid *)
numBeta              (* numeric β for the numerical check *)
```

For Metropolis acceptance use exactly:
```mathematica
If[RandomReal[] < MetropolisProb[dE], newState, state]
```
where `dE = energy[newState] - energy[state]`.

### Optional definitions

```mathematica
symParams = <|"eps" -> {...}, "couplings" -> {...}|>
DynamicSymParams[states_List] := ...   (* per-component symbolic parameters *)
DisplayState[state_] := ...            (* human-readable state string *)
ValidStateIDs[maxId_] := ...           (* restrict enumeration to valid IDs *)
```

### Abstract field/coupling functions (`vmmc_2d_field.wl` style)

`vmmc_2d_field.wl` separates symbolic and numerical modes:

**Symbolic check:** `fieldF` and `couplingJ` have no DownValues. Each evaluation `fieldF[x,y,L]` and `couplingJ[a,b,d2]` is treated as an independent free real atom, proving detailed balance for *all possible* functions simultaneously.

**Numerical MCMC / animation:** Concrete implementations (`$fieldFConcrete`, `$couplingJConcrete`) are activated inside a `Block`.

To customise, edit **Section 0 only** in `vmmc_2d_field.wl`:

```mathematica
$fieldFConcrete[x_Integer, y_Integer, L_Integer] :=
  fieldAmp * Sin[Pi * (x + 1) / L]

$couplingJConcrete[a_Integer, b_Integer, d2_Integer] :=
  $jPairSym[a, b] * Exp[-lambdaJ * d2]    (* must satisfy [a,b,d2] = [b,a,d2] *)

$concreteParams = <|fieldAmp -> 1.0, lambdaJ -> 0.5|>
$maxD2 = 1        (* max squared bond distance *)
$abstractFunctions = True
```

---

## Supported random calls

| Call | Behaviour |
|------|-----------|
| `RandomReal[]` / `Random[]` | Interval-tracking token |
| `RandomInteger[{lo, hi}]` | Rejection sampling over ⌈log₂(hi−lo+1)⌉ bits |
| `RandomInteger[n]` | Uniform over {0,…,n} |
| `RandomChoice[list]` | Uniform choice via rejection sampling |
| `RandomChoice[weights → elements]` | Sequential Bernoulli decomposition |

`RandomVariate`, `RandomSample`, `RandomPermutation` are **not** supported.

---

## Example files (`examples3/`)

| File | Description | DB |
|------|-------------|----|
| `kawasaki_1d.wl` | 1D Kawasaki (nearest-neighbour swap) on a periodic ring | PASS |
| `kawasaki_2d.wl` | 2D Kawasaki on a periodic square lattice | PASS |
| `vmmc_2d.wl` | Virtual Move Monte Carlo (Whitelam–Geissler) on a 2D torus | PASS |
| `vmmc_2d_field.wl` | VMMC with user-defined field and coupling functions | PASS |
| `jump_1d_weighted.wl` | 1D jump dynamics using `RandomChoice[weights → ...]` | PASS |
| `kawasaki_1d_fail.wl` | Sign-reversed dE in Metropolis | FAIL |
| `cluster_1d_fail.wl` | Cluster slides only rightward (asymmetric proposal) | FAIL |

---

## Quick examples

```bash
# Check 2D Kawasaki symbolically + numerically
wolframscript -file check.wls examples3/kawasaki_2d.wl

# Check VMMC with fast polynomial checker (faster for complex algorithms)
wolframscript -file check.wls examples3/vmmc_2d_field.wl \
  MaxBitString=1111111111 Mode=Symbolic FastChecker=1

# Report for a specific seed state
wolframscript -file report.wls examples3/vmmc_2d.wl BitString=101001 \
  ShowTrees=False ShowT=False

# Animate VMMC with field on a 6×6 grid
wolframscript -file animate.wls examples3/vmmc_2d_field.wl \
  Sites=36 N=10 Steps=500 Beta=1 FPS=15 Simple=1
```
