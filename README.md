# DetailedBalanceChecker

Symbolically proves or disproves detailed balance for Monte Carlo algorithms written in Mathematica. Uses computer algebra (FullSimplify with β > 0) rather than statistical tests alone.

---

## How it works

1. **Randomness → bits.** Every random call is intercepted and reduced to reads from a binary tape:
   - `RandomInteger[{lo,hi}]` / `RandomChoice[list]` — rejection sampling over ⌈log₂(n)⌉ bits.
   - `RandomReal[]` — *interval tracking*: a latent variable U ∈ [lo,hi]; each comparison `U < p` reads one bit and narrows the interval to [lo,p] or [p,hi], with weight (p−lo)/(hi−lo) or (hi−p)/(hi−lo) respectively.
   - `RandomChoice[weights → elements]` — decomposed into n−1 sequential Bernoulli trials using independent interval-tracking variables.

2. **BFS over bit sequences.** Starting from a seed state, all possible bit sequences are enumerated. Each complete path gives a (start, end, weight) triple. Weights for the same (start, end) pair are summed to form the symbolic transition matrix T.

3. **Symbolic detailed-balance check.** For each pair (i, j):
   ```
   FullSimplify[ T(i→j)·exp(−β E(i)) − T(j→i)·exp(−β E(j)), {β > 0} ]  =?  0
   ```
   β is kept as a free symbol so Boltzmann factors cancel algebraically.

4. **Numerical MCMC check.** The algorithm is run as a genuine Markov chain; the empirical distribution is compared to the Boltzmann distribution via KL divergence. KL < 0.02 is a pass.

---

## Usage

### Checker
```bash
wolframscript -file check.wls <algorithm.wl> [options]
```

| Option | Default | Description |
|--------|---------|-------------|
| `MaxBitString=XXXX` | `11111111` | Largest bit string tested; all strings up to this length are checked |
| `Mode=Symbolic\|Numerical\|Both` | `Both` | Which checks to run |
| `NSteps=N` | `50000` | MCMC steps for numerical check |
| `MaxBitDepth=N` | `20` | BFS depth cap per state |
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
| `PrintMatrix=True\|False` | `False` | Print transition matrix entries as LaTeX to terminal |

### Animation
```bash
wolframscript -file animate.wls <algorithm.wl> Sites=<n> N=<n> [options]
```

| Option | Default | Description |
|--------|---------|-------------|
| `Sites=N` | *(required)* | Total lattice sites |
| `N=N` | *(required)* | Number of labeled particle types |
| `Steps=N` | `200` | MCMC steps to run |
| `Beta=f` | from `.wl` file | Inverse temperature (overrides `numBeta`) |
| `FPS=f` | `10` | Animation frame rate |
| `Simple=1` | off | Fast 2-colour mode (holes vs particles) |
| `RecordEvery=N` | `1` | Record state every N steps |
| `Jpair<a><b>=f` | random | Per-pair coupling amplitude for `vmmc_2d_field` (e.g. `Jpair12=-1.0`) |

For `vmmc_2d_field.wl` all unspecified parameters (`fieldAmp`, `lambdaJ`, `Jpair12`, …) are assigned random values automatically. Pass `Beta=1` for a physically meaningful temperature.

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
where `dE = energy[newState] - energy[state]` and `MetropolisProb[dE]` is intercepted symbolically.

### Defining field and coupling functions (`vmmc_2d_field.wl` style)

`vmmc_2d_field.wl` supports user-defined external field and pair coupling functions. The design separates two modes:

**Symbolic check (default):** `fieldF` and `couplingJ` have **no DownValues**. Every site-value `fieldF[x,y,L]` and every pair-value `couplingJ[a,b,d2]` is treated as an independent free real atom by FullSimplify. This proves detailed balance for *all possible* field and coupling functions simultaneously, not just the specific forms you supply.

**Numerical MCMC / animation:** The concrete implementations (`$fieldFConcrete`, `$couplingJConcrete`) are activated inside a `Block`, with the parameters in `$concreteParams` assigned random values (or values you pass on the command line).

To customise the model, edit **Section 0 only** in `vmmc_2d_field.wl`:

```mathematica
(* Concrete field — used for numerical runs and animation only *)
$fieldFConcrete[x_Integer, y_Integer, L_Integer] :=
  fieldAmp * Sin[Pi * (x + 1) / L]

(* Concrete coupling — MUST satisfy $couplingJConcrete[a,b,d2] = $couplingJConcrete[b,a,d2] *)
$couplingJConcrete[a_Integer, b_Integer, d2_Integer] :=
  $jPairSym[a, b] * Exp[-lambdaJ * d2]

(* Numeric parameter values for MCMC and animation runs *)
$concreteParams = <|fieldAmp -> 1.0, lambdaJ -> 0.5|>

(* Maximum squared distance for bonds; 1 = nearest-neighbour only *)
$maxD2 = 1

(* Required — do not remove *)
$abstractFunctions = True
```

The checker verifies `$couplingJConcrete` symmetry before running. The `_Integer` patterns ensure matching is always exact (all coordinates and distances are integers).

### Other optional definitions

```mathematica
symParams = <|"eps" -> {...}, "couplings" -> {...}|>
DynamicSymParams[states_List] := ...   (* returns <|"couplings" -> {...}|> *)
DisplayState[state_] := ...
ValidStateIDs[maxId_] := ...
```

See `template.wl` for a fully annotated starting point.

---

## Supported random calls

| Call | Behaviour |
|------|-----------|
| `RandomReal[]` / `Random[]` | Interval-tracking token |
| `RandomInteger[{lo, hi}]` | Rejection sampling over ⌈log₂(hi−lo+1)⌉ bits |
| `RandomInteger[n]` | Uniform over {0,…,n} |
| `RandomChoice[list]` | Uniform choice via rejection sampling |
| `RandomChoice[weights → elements]` | Sequential Bernoulli decomposition |

`RandomVariate`, `RandomSample`, `RandomPermutation` are **not** supported and will abort the analysis.

---

## Example files (`examples3/`)

### Pass

| File | Description |
|------|-------------|
| `kawasaki_1d.wl` | 1D Kawasaki (nearest-neighbour swap) on a periodic ring |
| `kawasaki_2d.wl` | 2D Kawasaki on a periodic square lattice |
| `vmmc_2d.wl` | Virtual Move Monte Carlo (Whitelam–Geissler) on a 2D torus |
| `vmmc_2d_field.wl` | VMMC with user-defined field (`fieldF`) and coupling (`couplingJ`) functions. Example: sinusoidal field in x, exponential-decay coupling. Post-cluster Metropolis handles field energy changes. |
| `jump_1d_edit.wl` | 1D jump dynamics with rejection-sampled displacement |
| `jump_1d_weighted.wl` | Same jump dynamics using `RandomChoice[weights → ...]` |

### Fail (deliberate bugs for testing)

| File | Bug |
|------|-----|
| `kawasaki_1d_fail.wl` | Sign-reversed dE in Metropolis |
| `kawasaki_2d_fail.wl` | Same wrong-sign dE in 2D |
| `cluster_1d_fail.wl` | Cluster slides only rightward (asymmetric proposal) |
| `cluster_2d.wl` | Cluster merging changes cluster count |

---

## Quick examples

```bash
# Check kawasaki_1d symbolically + numerically
wolframscript -file check.wls examples3/kawasaki_1d.wl MaxBitString=1111111

# Report for a specific state (no slow notebook sections)
wolframscript -file report.wls examples3/vmmc_2d.wl BitString=1001 \
  ShowTrees=False ShowT=False PrintMatrix=True

# Animate vmmc_2d_field on a 2×2 grid
wolframscript -file animate.wls examples3/vmmc_2d_field.wl \
  Sites=4 N=2 Steps=500 Beta=1 FPS=15 Jpair12=-1.0
```
